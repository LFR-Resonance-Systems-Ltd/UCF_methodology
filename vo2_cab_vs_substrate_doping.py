from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def _read_xy(path: Path) -> tuple[list[float], list[float]]:
    xs: list[float] = []
    ys: list[float] = []
    with path.open("r", encoding="utf-8", errors="replace", newline="") as f:
        for row in csv.reader(f):
            if not row or len(row) < 2:
                continue
            try:
                x = float(str(row[0]).strip())
                y = float(str(row[1]).strip())
            except Exception:
                continue
            xs.append(x)
            ys.append(y)

    pts = sorted(zip(xs, ys), key=lambda t: t[0])
    return [p[0] for p in pts], [p[1] for p in pts]


def _clip01(v: float) -> float:
    return 0.0 if v < 0 else (1.0 if v > 1 else v)


def _interp1(xq: float, xs: list[float], ys: list[float]) -> float:
    if xq <= xs[0]:
        return ys[0]
    if xq >= xs[-1]:
        return ys[-1]

    lo, hi = 0, len(xs) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if xs[mid] <= xq:
            lo = mid
        else:
            hi = mid

    x0, x1 = xs[lo], xs[hi]
    y0, y1 = ys[lo], ys[hi]
    if x1 == x0:
        return y0
    t = (xq - x0) / (x1 - x0)
    return y0 + t * (y1 - y0)


def _pearson_r(a: list[float], b: list[float]) -> float:
    n = len(a)
    if n != len(b) or n < 2:
        return float("nan")
    ma = sum(a) / n
    mb = sum(b) / n
    va = sum((x - ma) ** 2 for x in a)
    vb = sum((y - mb) ** 2 for y in b)
    if va <= 0 or vb <= 0:
        return float("nan")
    cov = sum((a[i] - ma) * (b[i] - mb) for i in range(n))
    return cov / math.sqrt(va * vb)


def _t_at_fraction(xs: list[float], ys: list[float], frac: float = 0.5) -> float:
    target = frac
    for i in range(1, len(xs)):
        y0, y1 = ys[i - 1], ys[i]
        if (y0 - target) * (y1 - target) <= 0:
            x0, x1 = xs[i - 1], xs[i]
            if y1 == y0:
                return (x0 + x1) / 2
            t = (target - y0) / (y1 - y0)
            return x0 + t * (x1 - x0)
    return float("nan")


@dataclass(frozen=True)
class Condition:
    doping: str  # 'undoped' | 'Cr-doped'
    substrate: str  # 'c-plane' | 'r-plane'
    branch: str  # 'up' | 'down'


def compute_metrics(
    rr_path: Path,
    mf_path: Path,
    t_window: tuple[float, float] = (25.0, 100.0),
) -> dict[str, float]:
    t_lo, t_hi = t_window
    x_rr, rr = _read_xy(rr_path)
    x_mf, mf0 = _read_xy(mf_path)
    rr = [_clip01(v) for v in rr]
    mf0 = [_clip01(v) for v in mf0]

    # Sample metallic fraction onto rr temperature grid, then restrict to window.
    xs = [x for x in x_rr if t_lo <= x <= t_hi]
    rr_s = [_interp1(x, x_rr, rr) for x in xs]
    mf_s = [_interp1(x, x_mf, mf0) for x in xs]

    c_raw = _pearson_r(rr_s, mf_s)

    t50_rr = _t_at_fraction(x_rr, [_interp1(x, x_rr, rr) for x in x_rr], 0.5)
    t50_mf = _t_at_fraction(x_rr, [_interp1(x, x_mf, mf0) for x in x_rr], 0.5)
    dt50 = t50_mf - t50_rr

    # Window-aligned coherence: compare rr(T) vs mf(T + dt50) on overlap.
    x0 = max(min(xs), min(x_mf) - dt50)
    x1 = min(max(xs), max(x_mf) - dt50)
    xs_ov = [x for x in xs if x0 <= x <= x1]
    rr_ov = [_interp1(x, x_rr, rr) for x in xs_ov]
    mf_ov = [_interp1(x + dt50, x_mf, mf0) for x in xs_ov]
    c_aligned = _pearson_r(rr_ov, mf_ov)

    return {
        "n_points": float(len(xs)),
        "C_raw": float(c_raw),
        "C_aligned": float(c_aligned),
        "T50_r_ratio": float(t50_rr),
        "T50_metallic": float(t50_mf),
        "dT50_metal_minus_r": float(dt50),
        "t_window_lo": float(t_lo),
        "t_window_hi": float(t_hi),
    }


def build_table(t_window: tuple[float, float] = (25.0, 100.0)) -> pd.DataFrame:
    base = Path(__file__).resolve().parent.parent / "data" / "Figure7"
    panel_map: dict[Condition, tuple[Path, Path]] = {
        # undoped
        Condition("undoped", "c-plane", "up"): (base / "a_r_ratio_up.csv", base / "a_c_metallic_fraction_up.csv"),
        Condition("undoped", "c-plane", "down"): (base / "a_r_ratio_down.csv", base / "a_c_metallic_fraction_down.csv"),
        Condition("undoped", "r-plane", "up"): (base / "b_r_ratio_up.csv", base / "b_c_metallic_fraction_up.csv"),
        Condition("undoped", "r-plane", "down"): (base / "b_r_ratio_down.csv", base / "b_c_metallic_fraction_down.csv"),
        # Cr-doped (filenames as exported)
        Condition("Cr-doped", "c-plane", "up"): (base / "c-(r_ratio_up.csv", base / "c-c metallic fraction_up.csv"),
        Condition("Cr-doped", "c-plane", "down"): (base / "c-(r_ratio_down).csv", base / "c-(c metallic fraction_down).csv"),
        Condition("Cr-doped", "r-plane", "up"): (base / "d-(r_ratio_up).csv", base / "d-c metallic fraction_up.csv"),
        Condition("Cr-doped", "r-plane", "down"): (base / "d-(r_ratio_down).csv", base / "d-(c metallic fraction_down).csv"),
    }

    rows: list[dict[str, object]] = []
    for cond, (rr, mf) in panel_map.items():
        m = compute_metrics(rr, mf, t_window=t_window)
        rows.append(
            {
                "doping": cond.doping,
                "substrate": cond.substrate,
                "branch": cond.branch,
                "r_ratio_csv": rr.as_posix(),
                "metallic_fraction_csv": mf.as_posix(),
                **m,
            }
        )

    df = pd.DataFrame(rows)
    df["n_points"] = df["n_points"].astype(int)
    return df.sort_values(["doping", "substrate", "branch"], kind="stable")


def render_figure(df: pd.DataFrame, out_png: Path) -> None:
    # Keep colors conservative: use default matplotlib cycle.
    df = df.copy()

    def _label(row: pd.Series) -> str:
        return f"{row['doping']} | {row['substrate']} | {row['branch']}"

    df["label"] = df.apply(_label, axis=1)

    order = [
        "undoped | c-plane | up",
        "undoped | c-plane | down",
        "undoped | r-plane | up",
        "undoped | r-plane | down",
        "Cr-doped | c-plane | up",
        "Cr-doped | c-plane | down",
        "Cr-doped | r-plane | up",
        "Cr-doped | r-plane | down",
    ]
    df["label"] = pd.Categorical(df["label"], categories=order, ordered=True)
    df = df.sort_values("label")

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.2), constrained_layout=True)

    x = list(range(len(df)))
    labels = [str(v) for v in df["label"].tolist()]

    # Panel 1: dT50
    axes[0].bar(x, df["dT50_metal_minus_r"], linewidth=0)
    axes[0].axhline(0.0, color="#888888", linewidth=1)
    axes[0].set_title(r"$\Delta T_{50} = T_{50}^{\mathrm{metal}} - T_{50}^{\mathrm{rutile}}$")
    axes[0].set_ylabel("°C")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(labels, rotation=55, ha="right", fontsize=8)

    # Panel 2: C_aligned (primary), with faint markers for C_raw
    axes[1].bar(x, df["C_aligned"], linewidth=0)
    axes[1].scatter(x, df["C_raw"], s=18, color="#333333", alpha=0.65, label=r"$\mathcal{C}_{\mathrm{raw}}$")
    axes[1].axhline(0.8, color="#cc4444", linestyle="--", linewidth=1.5, label="threshold 0.80")
    axes[1].set_ylim(0.0, 1.02)
    axes[1].set_title(r"Cross-channel coherence ($\mathcal{C}$)")
    axes[1].set_ylabel("Pearson r")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(labels, rotation=55, ha="right", fontsize=8)
    axes[1].legend(frameon=False, fontsize=8, loc="lower left")

    fig.suptitle(r"VO₂: $\mathcal{C}_{AB}$ vs substrate/strain proxy and Cr doping (digitised from Suleiman et al., 2021 Fig.7)")

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


def main() -> None:
    # Default window uses the full digitised axis. If you want strict comparability with a narrower
    # transition region, change t_window below.
    t_window = (25.0, 100.0)

    df = build_table(t_window=t_window)

    out_dir = Path(__file__).resolve().parent.parent / "results"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_csv = out_dir / "VO2_CAB_vs_substrate_doping.csv"
    out_png = out_dir / "VO2_CAB_vs_substrate_doping.png"

    df.to_csv(out_csv, index=False)
    render_figure(df, out_png)

    print(f"Wrote: {out_csv}")
    print(f"Wrote: {out_png}")


if __name__ == "__main__":
    main()
