# Quick Start: UCF VO₂ Demo

## 1. Prerequisites
```bash
pip install numpy pandas scipy matplotlib
```

## 2. Run the analysis

### Main VO₂ transport + cross-channel coherence analysis
```bash
cd scripts
python UCF_demo1_analysis.py
```

**Outputs**:
- `../results/UCF_VO2_complete.png` — 9-panel figure
- `../results/UCF_demo1_results.txt` — numerical results

**Expected result**: All four conditions show $\mathcal{C}_{AB} > 0.80$; r-plane achieves 0.9950 (near-perfect).

### Substrate/doping comparison
```bash
python vo2_cab_vs_substrate_doping.py
```

**Outputs**:
- `../results/VO2_CAB_vs_substrate_doping.png` — 2-panel plot
- `../results/VO2_CAB_vs_substrate_doping.csv` — metrics table

**Key finding**: Cr-doped samples show reversed $\Delta T_{50}$ sign (substrate strain effect).

---

## 3. Understand the results

See [`results/README.md`](results/README.md) for detailed explanation of output files.

See [`METHOD.md`](METHOD.md) for complete technical description of:
- The $\mathcal{C}_{AB}$ definition
- Parameter choices (smoothing, window, etc.)
- Interpretation guidelines

---

## 4. Explore the data

All data are in `data/Figure4/` and `data/Figure7/` (digitized from Suleiman et al., 2021).

See [`data/README.md`](data/README.md) for:
- File organization
- Format and accuracy
- Data extraction method
- Limitations

---

## 5. Reproduce vs modify

### To reproduce exactly:
```bash
cd scripts
python UCF_demo1_analysis.py
```
(Outputs overwrite existing results; all random seeds fixed.)

### To experiment:
Edit configuration in `scripts/UCF_demo1_analysis.py`:
```python
SIGMA      = 5           # Gaussian smoothing sigma (try 3–10)
T_MIN, T_MAX = 20, 100   # Temperature range
N_GRID     = 600         # Grid resolution
CORR_T_LOW, CORR_T_HIGH = 40, 90  # Coherence window (try different ranges)
```

Then rerun. All changes propagate to outputs.

---

## 6. Next steps

- **Add a new material system**: Digitize Raman + transport data for another MIT candidate (TiO₂, NdNiO₃, etc.)
  - Place CSVs in a new `data/Figure_new_material/` folder
  - Create a new analysis script following the same template
  
- **Improve the coherence metric**: Experiment with alternative definitions
  - Spearman rank correlation instead of Pearson
  - Adaptive window based on transition onset/end
  - Information-theoretic measures (mutual information, transfer entropy)

- **Access raw data**: If you have the original experimental raw data (not digitized), replace the CSV files with your own and rerun

---

## 7. Get help

- **Method details**: See [`METHOD.md`](METHOD.md)
- **Data source**: Suleiman et al. (2021), Scientific Reports **11**, 1620
- **Code issues**: Check the scripts in [`scripts/`](scripts/) for comments and parameter descriptions

---

*This is active research methodology — not a finished product. The goal is to make the path real and improvable.*
