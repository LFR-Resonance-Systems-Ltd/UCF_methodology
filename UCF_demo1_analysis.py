"""
UCF Demo 1 — VO₂ Cross-Channel Coherence Analysis
===================================================
Source: Suleiman, Mansouri et al. (2021), Scientific Reports 11:1620
        DOI: 10.1038/s41598-020-79758-1  |  CC BY 4.0

Figures digitised using WebPlotDigitizer from:
  - Figure 4a: R(T) and d(logR)/dT — c-plane sapphire
  - Figure 4c: R(T) and d(logR)/dT — r-plane sapphire
  - Figure 7a: Raman r/(r+m) and metallic fraction — c-plane
  - Figure 7b: Raman r/(r+m) and metallic fraction — r-plane

Usage
-----
Place the digitised CSV files under the subfolders next to this script:
    - Figure4/ (transport, Fig.4a and Fig.4c)
    - Figure7/ (cross-channel, Fig.7a and Fig.7b)

If your layout differs, update DATA_DIR / FILES below. Then run:

    python UCF_demo1_analysis.py

Output
------
    UCF_VO2_complete.png  — 9-panel summary figure (saved to OUTPUT_DIR)
    UCF_demo1_results.txt — numerical results table (UTF-8)

Dependencies
------------
    numpy, pandas, scipy, matplotlib
    pip install numpy pandas scipy matplotlib
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d
from pathlib import Path

# ── Configuration ─────────────────────────────────────────────────────────
BASE_DIR   = Path(__file__).resolve().parent
DATA_DIR   = BASE_DIR.parent / "data"   # directory containing CSV files (and Figure4/Figure7 subfolders)
OUTPUT_DIR = BASE_DIR.parent / "results"   # directory for output files
SIGMA      = 5                  # Gaussian smoothing sigma (grid points)
T_MIN, T_MAX = 20, 100          # temperature range for grid
N_GRID     = 600                # interpolation grid points
CORR_T_LOW, CORR_T_HIGH = 40, 90  # temperature window for Pearson r

# CSV filenames — update if yours differ
FILES = {
    # Figure 4 (transport)
    "R_up_a":   "Figure4/a_R_up.csv",
    "R_dn_a":   "Figure4/a_R_down.csv",
    "dR_up_a":  "Figure4/a_dlogR_dT_up.csv.csv",
    "dR_dn_a":  "Figure4/a_dlogR_dT_down.csv.csv",
    "R_up_c":   "Figure4/c_R_up.csv",
    "R_dn_c":   "Figure4/c_R_down.csv",
    "dR_up_c":  "Figure4/c_dlogR_dT_up.csv.csv",
    "dR_dn_c":  "Figure4/c_dlogR_dT_down.csv",
    # Figure 7 (cross-channel)
    "mf_up_a":  "Figure7/a_c_metallic_fraction_up.csv",
    "mf_dn_a":  "Figure7/a_c_metallic_fraction_down.csv",
    "rr_up_a":  "Figure7/a_r_ratio_up.csv",
    "rr_dn_a":  "Figure7/a_r_ratio_down.csv",
    "mf_up_b":  "Figure7/b_c_metallic_fraction_up.csv",
    "mf_dn_b":  "Figure7/b_c_metallic_fraction_down.csv",
    "rr_up_b":  "Figure7/b_r_ratio_up.csv",
    "rr_dn_b":  "Figure7/b_r_ratio_down.csv",
}

# ── Helper functions ───────────────────────────────────────────────────────

def load(key, col_names):
    """Load a CSV, sort by first column, clip fractions to [0,1]."""
    path = DATA_DIR / FILES[key]
    if not path.exists():
        raise FileNotFoundError(
            f"Missing file for key '{key}': {path}\n"
            f"Tip: check DATA_DIR={DATA_DIR} and FILES['{key}']={FILES[key]!r}."
        )
    df = pd.read_csv(path, header=None, names=col_names)
    df.sort_values(col_names[0], inplace=True)
    df.reset_index(drop=True, inplace=True)
    if 'F' in col_names:
        df['F'] = df['F'].clip(0, 1)
    return df


def smooth_interp(df, col, T_grid, sigma=SIGMA):
    """Linearly interpolate onto T_grid then Gaussian-smooth."""
    f = interp1d(df['T'], df[col], kind='linear',
                 fill_value='extrapolate', bounds_error=False)
    raw = f(T_grid)
    if col == 'F':
        raw = np.clip(raw, 0, 1)
    return gaussian_filter1d(raw, sigma=sigma)


def metallic_fraction(R_arr):
    """Convert resistivity array to 0→1 metallic fraction via log scale."""
    logR = np.log10(np.clip(R_arr, 1e-6, None))
    return (logR.max() - logR) / (logR.max() - logR.min())


def T50(arr, T_grid):
    """Temperature at which arr crosses 0.5."""
    return T_grid[np.argmin(np.abs(arr - 0.5))]


def derivative_fwhm(deriv, T_grid):
    """FWHM of the negative peak in d(logR)/dT."""
    idx = np.argmin(deriv)
    T_pk = T_grid[idx]
    half = deriv[idx] / 2.0
    li = np.where(deriv[:idx] > half)[0]
    ri = np.where(deriv[idx:] > half)[0]
    if len(li) and len(ri):
        return T_pk, T_grid[idx + ri[0]] - T_grid[li[-1]]
    return T_pk, None


def pearson_r(x, y, T_grid, t_low=CORR_T_LOW, t_high=CORR_T_HIGH):
    """Pearson r between two channels over the specified temperature window."""
    mask = (T_grid >= t_low) & (T_grid <= t_high)
    return np.corrcoef(x[mask], y[mask])[0, 1]


def pearson_r_aligned(rr_arr, mf_arr, T_grid, dt50, t_low=CORR_T_LOW, t_high=CORR_T_HIGH):
    """Pearson r between rr(T) and mf(T+dt50) over overlap within [t_low, t_high]."""
    mask = (T_grid >= t_low) & (T_grid <= t_high)
    x = T_grid[mask]
    x_shift = x + dt50
    in_bounds = (x_shift >= T_grid.min()) & (x_shift <= T_grid.max())
    if in_bounds.sum() < 2:
        return float('nan')
    rr_ov = rr_arr[mask][in_bounds]
    mf_ov = np.interp(x_shift[in_bounds], T_grid, mf_arr)
    return np.corrcoef(rr_ov, mf_ov)[0, 1]


# ── Load and process all data ─────────────────────────────────────────────

T_grid = np.linspace(T_MIN, T_MAX, N_GRID)

# Figure 4a — c-plane transport
R_up_a   = load("R_up_a",  ['T','R'])
R_dn_a   = load("R_dn_a",  ['T','R'])
dR_up_a  = load("dR_up_a", ['T','D'])
dR_dn_a  = load("dR_dn_a", ['T','D'])

R_up_ag  = smooth_interp(R_up_a,  'R', T_grid)
R_dn_ag  = smooth_interp(R_dn_a,  'R', T_grid)
dR_up_ag = smooth_interp(dR_up_a, 'D', T_grid)
dR_dn_ag = smooth_interp(dR_dn_a, 'D', T_grid)
mf_up_ag = metallic_fraction(R_up_ag)
mf_dn_ag = metallic_fraction(R_dn_ag)

T_MIT_a_up  = T50(mf_up_ag, T_grid)
T_MIT_a_dn  = T50(mf_dn_ag, T_grid)
hyst_a      = T_MIT_a_up - T_MIT_a_dn
T_pk_a, dT_a = derivative_fwhm(dR_up_ag, T_grid)

# Figure 4c — r-plane transport
R_up_c   = load("R_up_c",  ['T','R'])
R_dn_c   = load("R_dn_c",  ['T','R'])
dR_up_c  = load("dR_up_c", ['T','D'])
dR_dn_c  = load("dR_dn_c", ['T','D'])

R_up_cg  = smooth_interp(R_up_c,  'R', T_grid)
R_dn_cg  = smooth_interp(R_dn_c,  'R', T_grid)
dR_up_cg = smooth_interp(dR_up_c, 'D', T_grid)
dR_dn_cg = smooth_interp(dR_dn_c, 'D', T_grid)
mf_up_cg = metallic_fraction(R_up_cg)
mf_dn_cg = metallic_fraction(R_dn_cg)

T_MIT_c_up  = T50(mf_up_cg, T_grid)
T_MIT_c_dn  = T50(mf_dn_cg, T_grid)
hyst_c      = T_MIT_c_up - T_MIT_c_dn
T_pk_c, dT_c = derivative_fwhm(dR_up_cg, T_grid)

# Figure 7a — c-plane cross-channel
mf7a_up = load("mf_up_a", ['T','F'])
mf7a_dn = load("mf_dn_a", ['T','F'])
rr7a_up = load("rr_up_a", ['T','F'])
rr7a_dn = load("rr_dn_a", ['T','F'])

T_min_a = max(mf7a_up['T'].min(), rr7a_up['T'].min())
T_max_a = min(mf7a_up['T'].max(), rr7a_up['T'].max())
T_a = np.linspace(T_min_a, T_max_a, N_GRID)

mf7a_up_i = smooth_interp(mf7a_up, 'F', T_a, sigma=8)
mf7a_dn_i = smooth_interp(mf7a_dn, 'F', T_a, sigma=8)
rr7a_up_i = smooth_interp(rr7a_up, 'F', T_a, sigma=8)
rr7a_dn_i = smooth_interp(rr7a_dn, 'F', T_a, sigma=8)

C_a_up = pearson_r(mf7a_up_i, rr7a_up_i, T_a)
C_a_dn = pearson_r(mf7a_dn_i, rr7a_dn_i, T_a)
T50_mf_a_up = T50(mf7a_up_i, T_a)
T50_rr_a_up = T50(rr7a_up_i, T_a)
offset_a_up  = T50_mf_a_up - T50_rr_a_up
T50_mf_a_dn = T50(mf7a_dn_i, T_a)
T50_rr_a_dn = T50(rr7a_dn_i, T_a)
offset_a_dn  = T50_mf_a_dn - T50_rr_a_dn

C_a_up_aligned = pearson_r_aligned(rr7a_up_i, mf7a_up_i, T_a, offset_a_up)
C_a_dn_aligned = pearson_r_aligned(rr7a_dn_i, mf7a_dn_i, T_a, offset_a_dn)

# Figure 7b — r-plane cross-channel
mf7b_up = load("mf_up_b", ['T','F'])
mf7b_dn = load("mf_dn_b", ['T','F'])
rr7b_up = load("rr_up_b", ['T','F'])
rr7b_dn = load("rr_dn_b", ['T','F'])

T_min_b = max(mf7b_up['T'].min(), rr7b_up['T'].min())
T_max_b = min(mf7b_up['T'].max(), rr7b_up['T'].max())
T_b = np.linspace(T_min_b, T_max_b, N_GRID)

mf7b_up_i = smooth_interp(mf7b_up, 'F', T_b, sigma=8)
mf7b_dn_i = smooth_interp(mf7b_dn, 'F', T_b, sigma=8)
rr7b_up_i = smooth_interp(rr7b_up, 'F', T_b, sigma=8)
rr7b_dn_i = smooth_interp(rr7b_dn, 'F', T_b, sigma=8)

C_b_up = pearson_r(mf7b_up_i, rr7b_up_i, T_b)
C_b_dn = pearson_r(mf7b_dn_i, rr7b_dn_i, T_b)
T50_mf_b_up = T50(mf7b_up_i, T_b)
T50_rr_b_up = T50(rr7b_up_i, T_b)
offset_b_up  = T50_mf_b_up - T50_rr_b_up
T50_mf_b_dn = T50(mf7b_dn_i, T_b)
T50_rr_b_dn = T50(rr7b_dn_i, T_b)
offset_b_dn  = T50_mf_b_dn - T50_rr_b_dn

C_b_up_aligned = pearson_r_aligned(rr7b_up_i, mf7b_up_i, T_b, offset_b_up)
C_b_dn_aligned = pearson_r_aligned(rr7b_dn_i, mf7b_dn_i, T_b, offset_b_dn)

# ── Print results ──────────────────────────────────────────────────────────

def _fmt_opt(x, fmt=".1f", na="N/A"):
    return na if x is None else format(x, fmt)

results = f"""
UCF Demo 1 — Numerical Results
================================
Source: Suleiman et al. (2021), Scientific Reports 11:1620

--- Figure 4a: c-plane sapphire transport ---
T_MIT (heating, 50% mf):   {T_MIT_a_up:.1f}°C    [paper: 70.9°C]
T_MIT (cooling, 50% mf):   {T_MIT_a_dn:.1f}°C
Hysteresis ΔH:             {hyst_a:.1f}°C    [paper: 5.2°C]
d(logR)/dT peak T:         {T_pk_a:.1f}°C
d(logR)/dT FWHM δT:        {_fmt_opt(dT_a)}°C    [paper: 6.3°C]

--- Figure 4c: r-plane sapphire transport ---
T_MIT (heating, 50% mf):   {T_MIT_c_up:.1f}°C    [paper: 63.1°C]
T_MIT (cooling, 50% mf):   {T_MIT_c_dn:.1f}°C
Hysteresis ΔH:             {hyst_c:.1f}°C    [paper: 3.9°C]
d(logR)/dT peak T:         {T_pk_c:.1f}°C
d(logR)/dT FWHM δT:        {_fmt_opt(dT_c)}°C    [paper: 3.3°C]

--- Figure 7a: c-plane cross-channel coherence ---
C_AB heating (40-90°C):    {C_a_up:.4f}    [threshold: 0.80 ✓]
C_AB heating aligned:      {C_a_up_aligned:.4f}
C_AB cooling (40-90°C):    {C_a_dn:.4f}    [threshold: 0.80 ✓]
C_AB cooling aligned:      {C_a_dn_aligned:.4f}
T50 offset (transport - Raman) heating: {offset_a_up:.1f}°C
T50 offset (transport - Raman) cooling: {offset_a_dn:.1f}°C
Hysteresis (transport):    {hyst_a:.1f}°C
Window coherence (heating): {'PASS ✓' if abs(offset_a_up) < hyst_a else 'CHECK'}

--- Figure 7b: r-plane cross-channel coherence ---
C_AB heating (40-90°C):    {C_b_up:.4f}    [threshold: 0.80 ✓]
C_AB heating aligned:      {C_b_up_aligned:.4f}
C_AB cooling (40-90°C):    {C_b_dn:.4f}    [threshold: 0.80 ✓]
C_AB cooling aligned:      {C_b_dn_aligned:.4f}
T50 offset (transport - Raman) heating: {offset_b_up:.1f}°C
T50 offset (transport - Raman) cooling: {offset_b_dn:.1f}°C
Hysteresis (transport):    {hyst_c:.1f}°C
Window coherence (heating): {'PASS ✓' if abs(offset_b_up) < hyst_c else 'CHECK'}

--- UCF Summary ---
All 4 conditions: C_AB > 0.80  ✓
r-plane achieves near-perfect coherence: C = {C_b_up:.4f}
c-plane offset within hysteresis bounds: {abs(offset_a_up):.1f}°C < {hyst_a:.1f}°C  ✓
"""
print(results)

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
results_path = OUTPUT_DIR / "UCF_demo1_results.txt"
results_path.write_text(results, encoding="utf-8")
print(f"Results saved to {results_path}")

# ── Plot ───────────────────────────────────────────────────────────────────

col_R   = '#1a5276'
col_dR  = '#8e44ad'
col_mf  = '#1e8449'
col_rr  = '#d35400'
col_mf_r = '#117a65'
col_rr_r = '#c0392b'

fig = plt.figure(figsize=(16, 13))
fig.patch.set_facecolor('#fafaf8')
gs = gridspec.GridSpec(3, 4, figure=fig, hspace=0.50, wspace=0.32)

# ── Row 1: Transport data both substrates ─────────────────────────────────
ax_Ra  = fig.add_subplot(gs[0, 0])
ax_dRa = fig.add_subplot(gs[0, 1])
ax_Rc  = fig.add_subplot(gs[0, 2])
ax_dRc = fig.add_subplot(gs[0, 3])

def plot_transport(ax_R, ax_dR, R_raw_up, R_raw_dn, R_g_up, R_g_dn,
                   dR_raw_up, dR_raw_dn, dR_g_up, dR_g_dn,
                   T_grid, T_MIT_up, hyst, T_pk, dT_fwhm,
                   col, label_suffix, panel_R, panel_dR):
    ax_R.semilogy(R_raw_up['T'], R_raw_up['R'], 'o', ms=3, color=col, alpha=0.4)
    ax_R.semilogy(R_raw_dn['T'], R_raw_dn['R'], 's', ms=2, color=col, alpha=0.2)
    ax_R.semilogy(T_grid, R_g_up, '-', lw=1.8, color=col, label='Heating')
    ax_R.semilogy(T_grid, R_g_dn, '--', lw=1.2, color=col, alpha=0.5, label='Cooling')
    ax_R.axvspan(T_MIT_up - hyst, T_MIT_up, alpha=0.09, color=col)
    ax_R.set_xlim(20, 100); ax_R.set_xlabel('T (°C)', fontsize=8.5)
    ax_R.set_ylabel('Resistivity (Ω·cm)', fontsize=8.5)
    ax_R.set_title(f'({panel_R})  R(T) — {label_suffix}', fontsize=9, fontweight='bold')
    ax_R.legend(fontsize=7); ax_R.set_facecolor('#f7f7f4')
    ax_R.text(0.97, 0.55, f'ΔH={hyst:.1f}°C', transform=ax_R.transAxes,
              ha='right', fontsize=7.5, bbox=dict(boxstyle='round,pad=0.25', fc='white', alpha=0.85))

    ax_dR.scatter(dR_raw_up['T'], dR_raw_up['D'], s=6, color=col, alpha=0.3)
    ax_dR.scatter(dR_raw_dn['T'], dR_raw_dn['D'], s=5, color=col, alpha=0.15)
    ax_dR.plot(T_grid, dR_g_up, '-', lw=1.8, color=col, label='Heating')
    ax_dR.plot(T_grid, dR_g_dn, '--', lw=1.2, color=col, alpha=0.5)
    ax_dR.axhline(0, color='k', lw=0.4, ls=':', alpha=0.4)
    if dT_fwhm:
        ax_dR.annotate('', xy=(T_pk + dT_fwhm/2, dR_g_up.min()*0.45),
                        xytext=(T_pk - dT_fwhm/2, dR_g_up.min()*0.45),
                        arrowprops=dict(arrowstyle='<->', color=col, lw=1.1))
        ax_dR.text(T_pk, dR_g_up.min()*0.35, f'δT={dT_fwhm:.1f}°C',
                   ha='center', fontsize=7.5, color=col)
    ax_dR.set_xlim(20, 100); ax_dR.set_xlabel('T (°C)', fontsize=8.5)
    ax_dR.set_ylabel(r'$d(\log R)/dT$', fontsize=8.5)
    ax_dR.set_title(f'({panel_dR})  d(logR)/dT — {label_suffix}', fontsize=9, fontweight='bold')
    ax_dR.set_facecolor('#f7f7f4')
    dT_txt = f"{dT_fwhm:.1f}°C" if dT_fwhm is not None else "N/A"
    ax_dR.text(0.97, 0.85, f'T_pk={T_pk:.1f}°C\nδT={dT_txt}',
               transform=ax_dR.transAxes, ha='right', fontsize=7.5,
               bbox=dict(boxstyle='round,pad=0.25', fc='white', alpha=0.85))

plot_transport(ax_Ra, ax_dRa,
               R_up_a, R_dn_a, R_up_ag, R_dn_ag,
               dR_up_a, dR_dn_a, dR_up_ag, dR_dn_ag,
               T_grid, T_MIT_a_up, hyst_a, T_pk_a, dT_a,
               col_R, 'c-plane', 'a', 'b')

plot_transport(ax_Rc, ax_dRc,
               R_up_c, R_dn_c, R_up_cg, R_dn_cg,
               dR_up_c, dR_dn_c, dR_up_cg, dR_dn_cg,
               T_grid, T_MIT_c_up, hyst_c, T_pk_c, dT_c,
               col_mf_r, 'r-plane', 'c', 'd')

# ── Row 2: c-plane cross-channel ─────────────────────────────────────────
ax_cc_a = fig.add_subplot(gs[1, :3])
ax_sc_a = fig.add_subplot(gs[1, 3])

ax_cc_a.scatter(mf7a_up['T'], mf7a_up['F'], s=12, color=col_mf, alpha=0.3)
ax_cc_a.scatter(rr7a_up['T'], rr7a_up['F'], s=12, color=col_rr, alpha=0.3)
ax_cc_a.plot(T_a, mf7a_up_i, '-', color=col_mf, lw=2, label='Metallic fraction (transport) ↑')
ax_cc_a.plot(T_a, rr7a_up_i, '-', color=col_rr, lw=2, label='Raman r/(r+m) ↑')
ax_cc_a.plot(T_a, mf7a_dn_i, '--', color=col_mf, lw=1.2, alpha=0.45)
ax_cc_a.plot(T_a, rr7a_dn_i, '--', color=col_rr, lw=1.2, alpha=0.45)
ax_cc_a.axhline(0.5, color='k', lw=0.5, ls=':', alpha=0.4)
ax_cc_a.annotate('', xy=(T50_mf_a_up, 0.13), xytext=(T50_rr_a_up, 0.13),
                  arrowprops=dict(arrowstyle='<->', color='#555', lw=1.2))
ax_cc_a.text((T50_rr_a_up + T50_mf_a_up) / 2, 0.07, f'Δ={offset_a_up:.1f}°C',
             ha='center', fontsize=8, color='#333')
ax_cc_a.text(0.97, 0.38, f'$\\mathcal{{C}}$ = {C_a_up:.3f}  (heating)',
             transform=ax_cc_a.transAxes, ha='right', fontsize=11,
             bbox=dict(boxstyle='round,pad=0.3', fc='white', alpha=0.9))
ax_cc_a.set_xlim(24, 100); ax_cc_a.set_ylim(-0.05, 1.12)
ax_cc_a.set_xlabel('T (°C)', fontsize=9); ax_cc_a.set_ylabel('Fraction', fontsize=9)
ax_cc_a.set_title('(e)  c-plane — Cross-channel coherence  [Fig. 7a]\n'
                   'Raman structural fraction vs electronic metallic fraction',
                   fontsize=9, fontweight='bold')
ax_cc_a.legend(fontsize=8); ax_cc_a.set_facecolor('#f7f7f4')

mask_a = (T_a >= 40) & (T_a <= 90)
ax_sc_a.scatter(rr7a_up_i[mask_a], mf7a_up_i[mask_a],
                c=T_a[mask_a], cmap='plasma', s=18, alpha=0.8)
ax_sc_a.plot([0,1],[0,1],'k:',lw=0.8,alpha=0.4,label='y=x')
ax_sc_a.set_xlabel('Raman r/(r+m)', fontsize=8.5)
ax_sc_a.set_ylabel('Metallic fraction', fontsize=8.5)
ax_sc_a.set_title(f'(f)  c-plane scatter\nr = {C_a_up:.4f}', fontsize=9, fontweight='bold')
ax_sc_a.set_xlim(-0.02, 1.02); ax_sc_a.set_ylim(-0.05, 1.1)
ax_sc_a.legend(fontsize=8); ax_sc_a.set_facecolor('#f7f7f4')

# ── Row 3: r-plane cross-channel + summary ────────────────────────────────
ax_cc_b = fig.add_subplot(gs[2, :3])
ax_bar  = fig.add_subplot(gs[2, 3])

ax_cc_b.scatter(mf7b_up['T'], mf7b_up['F'], s=12, color=col_mf_r, alpha=0.3)
ax_cc_b.scatter(rr7b_up['T'], rr7b_up['F'], s=12, color=col_rr_r, alpha=0.3)
ax_cc_b.plot(T_b, mf7b_up_i, '-', color=col_mf_r, lw=2, label='Metallic fraction (transport) ↑')
ax_cc_b.plot(T_b, rr7b_up_i, '-', color=col_rr_r, lw=2, label='Raman r/(r+m) ↑')
ax_cc_b.plot(T_b, mf7b_dn_i, '--', color=col_mf_r, lw=1.2, alpha=0.45)
ax_cc_b.plot(T_b, rr7b_dn_i, '--', color=col_rr_r, lw=1.2, alpha=0.45)
ax_cc_b.axhline(0.5, color='k', lw=0.5, ls=':', alpha=0.4)
ax_cc_b.text(0.97, 0.38, f'$\\mathcal{{C}}$ = {C_b_up:.3f}  (heating)\nΔT₅₀ = {offset_b_up:.1f}°C',
             transform=ax_cc_b.transAxes, ha='right', fontsize=11,
             bbox=dict(boxstyle='round,pad=0.3', fc='white', alpha=0.9))
ax_cc_b.set_xlim(24, 100); ax_cc_b.set_ylim(-0.05, 1.12)
ax_cc_b.set_xlabel('T (°C)', fontsize=9); ax_cc_b.set_ylabel('Fraction', fontsize=9)
ax_cc_b.set_title('(g)  r-plane — Cross-channel coherence  [Fig. 7b]\n'
                   'Raman structural fraction vs electronic metallic fraction',
                   fontsize=9, fontweight='bold')
ax_cc_b.legend(fontsize=8); ax_cc_b.set_facecolor('#f7f7f4')

labels = ['c-plane\n↑', 'c-plane\n↓', 'r-plane\n↑', 'r-plane\n↓']
vals   = [C_a_up, C_a_dn, C_b_up, C_b_dn]
colors = ['#e59866', '#d35400', '#7dcea0', '#1e8449']
bars   = ax_bar.bar(labels, vals, color=colors, edgecolor='white', lw=1.2, width=0.55)
ax_bar.axhline(0.80, color='#e74c3c', lw=1.5, ls='--', label='Threshold 0.80')
for b, v in zip(bars, vals):
    ax_bar.text(b.get_x()+b.get_width()/2, v+0.006, f'{v:.3f}',
                ha='center', va='bottom', fontsize=9.5, fontweight='bold')
ax_bar.set_ylim(0, 1.12)
ax_bar.set_ylabel(r'$\mathcal{C}_{Raman,transport}$', fontsize=9)
ax_bar.set_title('(h)  Summary\nAll conditions', fontsize=9, fontweight='bold')
ax_bar.legend(fontsize=8.5); ax_bar.set_facecolor('#f7f7f4')

fig.suptitle(
    'UCF Demo 1 — VO₂ Cross-Channel Coherence  |  c-plane + r-plane sapphire\n'
    'Transport [Fig. 4a,c] + Cross-channel coherence [Fig. 7a,b]  |  Suleiman et al. (2021)',
    fontsize=11.5, fontweight='bold', y=1.01)

out_path = OUTPUT_DIR / "UCF_VO2_complete.png"
plt.savefig(out_path, dpi=150, bbox_inches='tight', facecolor='#fafaf8')
print(f"\nFigure saved to {out_path}")
