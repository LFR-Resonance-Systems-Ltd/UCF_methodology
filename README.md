# UCF_methodology
A methodology for quantifying cross-channel coherence in phase transitions (e.g., Raman vs transport in VO₂), using window-aligned correlation and T50 offset.
# Unified Coherence Framework (UCF) — A Methodology for Cross-Channel Phase Transition Analysis

**Status**: Active research exploration. This is not a finished product, but a documented path of method development.

---

## Why This Framework?

In correlated electron materials undergoing metal-insulator transitions (MIT), different electronic and structural order parameters often change at different temperatures. **The central question**: how do we quantify whether these changes are genuinely coupled or merely coincidental?

Existing approaches either:
- Compare individual transition temperatures (loses information about shape/timing of the transition)
- Use generic correlation metrics (neglect the intrinsic time/temperature offsets that transitions naturally have)

**UCF's approach**: Define a window-aware coherence measure $\mathcal{C}_{AB}$ that:
- Captures correlation *across* channels (e.g., electronic vs. structural)
- Accounts for intrinsic offsets between their 50%-transition points
- Remains insensitive to arbitrarily smooth portions of the curves
- Provides a single interpretable metric: 0 = independent, 1 = perfectly tracked

---

## Current Demonstrations

### VO₂ on sapphire (Suleiman et al., 2021)
- **Data source**: Digitized from Suleiman, Mansouri et al., *Scientific Reports* **11**, 1620 (2021)
- **Finding**: c-plane and r-plane show $\mathcal{C}_{AB} > 0.95$, with r-plane achieving near-perfect coherence (0.9950)
- **Interpretation**: Transport (metallic fraction) and Raman structural order are tightly coupled across both crystal faces
- **Method**: [scripts/UCF_demo1_analysis.py](scripts/UCF_demo1_analysis.py)

### VO₂ with Cr doping
- **Data source**: Same paper, Fig.7c (Cr-doped c-plane) and Fig.7d (r-plane)
- **Preliminary findings**: Substrate-strain effects show distinct $\Delta T_{50}$ signatures; Cr doping produces sign-flip in offset trend
- **Method**: [scripts/vo2_cab_vs_substrate_doping.py](scripts/vo2_cab_vs_substrate_doping.py)

---

## Quick Start

### Prerequisites
```bash
pip install numpy pandas scipy matplotlib
```

### Run the analysis
```bash
cd scripts
python UCF_demo1_analysis.py
```

Outputs:
- `../results/UCF_VO2_complete.png` — 9-panel figure with transport and cross-channel coherence
- `../results/UCF_demo1_results.txt` — numerical results (UTF-8)

### Reproduce the substrate/doping summary
```bash
cd scripts
python vo2_cab_vs_substrate_doping.py
```

Outputs:
- `../results/VO2_CAB_vs_substrate_doping.csv` — summary table
- `../results/VO2_CAB_vs_substrate_doping.png` — comparative coherence plot

---

## Directory Structure

```
UCF_methodology/
├── README.md                          # This file
├── METHOD.md                          # Detailed method description & definitions
├── data/
│   ├── Figure4/                       # Digitized transport data (Suleiman et al., Fig.4a, 4c)
│   │   └── {a,c}_{R,dlogR_dT}_{up,down}.csv
│   └── Figure7/                       # Digitized Raman/metallic fraction data (Fig.7a-d)
│       └── {a,b,c,d}_{r_ratio,c_metallic_fraction}_{up,down}.csv
├── scripts/
│   ├── UCF_demo1_analysis.py          # Main VO₂ analysis (transport + cross-channel)
│   └── vo2_cab_vs_substrate_doping.py # Substrate/strain comparison & doping effects
└── results/
    ├── UCF_VO2_complete.png           # 9-panel summary figure
    ├── UCF_demo1_results.txt          # Numerical results
    ├── VO2_CAB_vs_substrate_doping.png # Comparative coherence plot
    └── VO2_CAB_vs_substrate_doping.csv # Summary metrics table
```

---

## Key Parameters & Choices

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| Coherence window | 40–90°C | Covers the transition region for both c- and r-plane; avoids extreme T where digitization noise dominates |
| Gaussian smoothing σ | 5–8 grid points | Reduces digitization artifacts without erasing transition sharpness |
| T_grid resolution | 600 points over [20, 100]°C | ~0.13°C per point; smooth enough for derivative calculations |
| T50 definition | Fraction = 0.5 | Standard electronic/structural transition midpoint |

For full justification, see [METHOD.md](METHOD.md).

---

## Known Limitations & Next Steps

### Current limitations
1. **Data source**: Digitized from published figures, not raw experimental files
   - Inherent uncertainty in coordinate extraction (~0.5–1% per axis)
   - No access to original error bars or measurement details
2. **Single material system**: Only VO₂ demonstrated
   - Need results from another MIT candidate (e.g., TiO₂, NdNiO₃, or HEA) to validate universality
3. **Preliminary Cr-doping results**: Available data incomplete
   - Fig.7c/d filenames contain export artifacts; not all conditions have complete up/down branches
4. **No theory**: $\mathcal{C}_{AB}$ is empirical
   - Would benefit from equilibrium statistical mechanics grounding

### Next steps
- [ ] Analyze TiO₂ anatase→rutile transition (another canonical MIT)
- [ ] Access raw (non-digitized) data for VO₂ to reduce uncertainty
- [ ] Complete Cr-doped VO₂ Raman data if available
- [ ] Develop minimal theoretical model for $\mathcal{C}_{AB}$ interpretation
- [ ] Compare against other cross-channel metrics (e.g., mutual information, transfer entropy)

---

## How to Cite This Work

Since this is ongoing research methodology, the standard citation is not yet finalized. If you use this framework, please cite the Suleiman et al. data and mention this repository:

> Methodology repository: https://github.com/[your-username]/UCF_methodology  
> Data source: Suleiman et al., *Scientific Reports* **11**, 1620 (2021), DOI: 10.1038/s41598-020-79758-1

---

## Contributing & Feedback

This is an active record of a real scientific exploration. If you:
- Find a bug in the scripts, open an issue
- Have data from another material system to test, please share
- Want to propose a different coherence definition, let's discuss it

**The goal is to make this path real and improvable, not to declare it finished.**

---

## License

This work is released under the [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) license, matching the license of the original Suleiman et al. (2021) data.
