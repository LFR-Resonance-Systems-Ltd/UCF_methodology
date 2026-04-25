# Analysis Results

## Files in this directory

### UCF_VO2_complete.png
Main 9-panel figure summarizing VO₂ cross-channel coherence analysis:
- **Panels a–b**: c-plane transport (R vs T, d(logR)/dT vs T)
- **Panels c–d**: r-plane transport
- **Panels e–f**: c-plane cross-channel coherence (metallic fraction vs Raman r-ratio)
- **Panels g–h**: r-plane cross-channel coherence
- **Panel h**: Summary bar chart of coherence $\mathcal{C}_{AB}$ across all conditions

**How to read**: All four conditions (c-plane ↑, c-plane ↓, r-plane ↑, r-plane ↓) show $\mathcal{C}_{AB} > 0.80$, with r-plane achieving near-perfect coherence (0.9950).

### UCF_demo1_results.txt
Numerical results table in UTF-8 text format. Contains:
- Transport metrics: $T_{50}$, hysteresis $\Delta H$, peak temperature, FWHM
- Cross-channel coherence: $\mathcal{C}_{AB}$ (raw and aligned), $\Delta T_{50}$ offsets
- Summary verdict for window coherence (offset within hysteresis)

### VO2_CAB_vs_substrate_doping.png
Comparative summary of substrate/strain effects and Cr doping:
- **Left panel**: $\Delta T_{50} = T_{50}^{\text{metal}} - T_{50}^{\text{rutile}}$ across 8 conditions
- **Right panel**: Coherence $\mathcal{C}_{AB}$ (bars) and raw $\mathcal{C}_{\text{raw}}$ (scatter) with 0.80 threshold line

**Key observation**: Cr-doped samples show reversed sign of $\Delta T_{50}$ compared to undoped, suggesting a strain-driven effect.

### VO2_CAB_vs_substrate_doping.csv
Tabular data (CSV, UTF-8) with one row per condition:
- Columns: doping, substrate, branch (up/down)
- Input CSV files used for each condition
- Metrics: n_points, C_raw, C_aligned, T50_r_ratio, T50_metallic, dT50_metal_minus_r
- Temperature window: [25°C, 100°C]

**Usage**: Load in spreadsheet or Python to inspect raw numbers or create alternative visualizations.

---

## How these results were generated

All outputs are **fully reproducible** from the digitized data:

```bash
# Generate VO₂ transport + cross-channel figure
cd ../scripts
python UCF_demo1_analysis.py
# → ../results/UCF_VO2_complete.png, UCF_demo1_results.txt

# Generate substrate/doping comparison
python vo2_cab_vs_substrate_doping.py
# → ../results/VO2_CAB_vs_substrate_doping.png, .csv
```

All scripts run in-place; no command-line arguments required.

---

## Notes on digitization uncertainty

These results are based on **digitized figures**, not raw experimental data. Uncertainties:
- Coordinate extraction: ±0.5–1.5% per axis (temperature ±0.5°C, fraction ±0.01)
- Impact on $\mathcal{C}_{AB}$: ±0.01–0.02 (qualitative findings robust)
- Impact on $T_{50}$: ±0.5–1°C (but $\Delta T_{50}$ offsets are typically 2–5°C)

**Interpretation**: Qualitative results (e.g., r-plane shows higher coherence, Cr-doped sign-flip) are reliable; quantitative precision should be quoted with ±5–10% caveats.

---

*Last generated: April 25, 2026*
