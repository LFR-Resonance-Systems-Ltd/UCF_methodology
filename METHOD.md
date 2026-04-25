# UCF Methodology: Detailed Description

## 1. Problem Setup

### Background
In materials undergoing metal-insulator transitions (MIT), multiple order parameters evolve with temperature:
- **Electronic**: resistivity $\rho(T)$ → metallic fraction $f_{\text{metal}}(T)$
- **Structural/Magnetic**: Raman spectral weight, magnetic moment, lattice parameter, etc.

**Observation**: These rarely transition at the *same* temperature.

### Existing Approaches & Their Limitations

#### Approach A: Direct T50 Comparison
Simply compare $T_{50}^{\text{elec}}$ vs. $T_{50}^{\text{struct}}$.  
**Problem**: Ignores transition shape, width, and reversibility. Two materials with different transition slopes but same $T_{50}$ would score equally.

#### Approach B: Generic Pearson Correlation
Compute $r = \text{Pearson}(f_{\text{metal}}, f_{\text{struct}})$ over some temperature window.  
**Problems**:
- Strongly window-dependent (results fragile if window chosen poorly)
- Doesn't account for natural offsets; two perfectly anti-correlated but offset curves still give spurious high r
- No interpretation: r=0.7 doesn't tell you whether coupling is physical or sampling artifact

#### Approach C: Mutual Information / Transfer Entropy
Information-theoretic measures of channel interdependence.  
**Appeal**: Model-free, captures nonlinear coupling.  
**Problems**:
- Computationally intensive for sparse digitized data
- Sensitive to binning choices and noise
- Hard to compare across different material systems

---

## 2. UCF Definition

We propose a **window-aligned coherence** measure $\mathcal{C}_{AB}$ defined as follows:

### Step 1: Prepare the data
Given two signals $A(T) = f_{\text{metal}}(T)$ and $B(T) = f_{\text{struct}}(T)$ over temperature range $[T_{\min}, T_{\max}]$:

1. **Read** raw $(T, \text{value})$ pairs from CSV
2. **Sort** by temperature
3. **Interpolate** onto a uniform $T$-grid: $T_i = T_{\min} + i \cdot \Delta T$, $i = 0, 1, \ldots, N$
4. **Smooth** with Gaussian filter (σ in grid units) to reduce digitization noise

### Step 2: Compute raw coherence
Over a coherence window $[t_{\text{low}}, t_{\text{high}}]$ (typically [40°C, 90°C] for VO₂):

$$\mathcal{C}_{\text{raw}} = \text{Pearson}(A[t_{\text{low}}, t_{\text{high}}], B[t_{\text{low}}, t_{\text{high}}])$$

This captures whether the two channels *track* each other in shape.

### Step 3: Compute T₅₀ offset
Interpolate the 50%-transition point for each channel:

$$T_{50}^A = T \text{ where } A(T) = 0.5$$
$$T_{50}^B = T \text{ where } B(T) = 0.5$$
$$\Delta T_{50} = T_{50}^A - T_{50}^B$$

### Step 4: Compute aligned coherence
Shift $B$ by $\Delta T_{50}$ to temporally align their transitions, then recompute coherence:

$$\mathcal{C}_{\text{aligned}} = \text{Pearson}(A[t_{\text{low}}, t_{\text{high}}], B(T + \Delta T_{50})[t_{\text{low}}, t_{\text{high}}])$$

If $\mathcal{C}_{\text{aligned}} > \mathcal{C}_{\text{raw}}$, the channels are *intrinsically* coherent but offset; if equal, there is no natural offset in the physics.

---

## 3. Interpretation

| $\mathcal{C}_{AB}$ | $\Delta T_{50}$ | Interpretation |
|-------------------|-----------------|---|
| 0.95–1.00 | < hysteresis | Strongly coupled; offset is within hysteresis scatter |
| 0.80–0.95 | Variable | Moderately coupled; worth investigating further |
| < 0.80 | — | Weakly coupled or independent; coupling questionable |

---

## 4. Parameter Choices for VO₂

### Coherence window: [40°C, 90°C]
- **Why**: VO₂ on sapphire transitions sharply around 68–70°C; window covers the full transition region
- **Why not wider**: Below 40°C and above 90°C, both the digitized data and physical transition are poorly defined
- **Robustness**: Results don't change qualitatively if window shifts to [35°C, 95°C]

### Gaussian smoothing σ = 5–8 grid points
- **Grid resolution**: 600 points over [20°C, 100°C] → Δ*T* ≈ 0.13°C per point
- **σ = 5 points**: ~0.65°C characteristic width
- **Effect**: Smooths out noisy digitization spikes without erasing transition sharpness (real transitions are broader than 0.65°C at this resolution)
- **Validation**: Results stable if σ ∈ [3, 10]; above 10, transition features start to blur

### T₅₀ definition (50% transition)
- **Standard choice** in MIT literature
- **Robustness**: T₃₀ or T₇₀ give qualitatively similar ΔT₅₀ offsets; C_AB remains stable

---

## 5. Digitization & Uncertainty

All data come from **WebPlotDigitizer** extraction of published figures (Suleiman et al., 2021).

### Sources of uncertainty
1. **Coordinate extraction**: ~0.5–1.5% per axis (temperature ±0.5°C, fraction ±0.01)
2. **Figure resolution**: Original figures are ~600×800 px; local digitization error ~1–2 px
3. **No error bars in source**: Original measurements may have their own uncertainty, not accessible here

### Impact on results
- **C_AB**: Robust to ±1% coordinate noise; effect on correlation ~0.01–0.02
- **T₅₀**: Affected by ~0.5–1°C, but ΔT₅₀ offsets are typically 2–5°C, so relative error ~10–20%
- **Conclusion**: Qualitative findings (e.g., "c-plane and r-plane show different offsets") are reliable; quantitative precision should be quoted with ±5–10% caveats

---

## 6. Implementation Notes

### Language & dependencies
- **Python 3.7+** (f-strings, pathlib)
- **NumPy**: numerical arrays, linear interpolation, statistics
- **Pandas**: data I/O and table manipulation
- **SciPy**: Gaussian filtering, Pearson correlation
- **Matplotlib**: visualization

### File I/O
- Input: CSV files with 2 columns (T, value), no header
- Output: PNG figures (150–200 dpi for arXiv/figure quality), CSV summary tables (UTF-8 encoding)

### Reproducibility
- All random seeds set explicitly (no randomness in analysis)
- Scripts accept only configuration parameters (smoothing σ, window, grid resolution) at module top; no command-line arguments
- Output filenames and paths hardcoded; scripts run in-place

---

## 7. Validation & Sanity Checks

### VO₂ demo sanity checks
| Check | Value | Expected | Status |
|-------|-------|----------|--------|
| c-plane C_AB (heating) | 0.9464 | > 0.80 (coupling exists) | ✓ |
| c-plane C_AB (cooling) | 0.9737 | > 0.80 | ✓ |
| r-plane C_AB (heating) | 0.9950 | near 1.0 (perfect) | ✓ |
| r-plane C_AB (cooling) | 0.9955 | near 1.0 | ✓ |
| c-plane ΔT₅₀ (heating) | 4.2°C | < hysteresis (5.2°C) | ✓ |
| r-plane ΔT₅₀ (heating) | 0.1°C | < hysteresis (4.3°C) | ✓ |
| C_aligned > C_raw | Yes (c-plane, cooling) | Expected if offset is natural | ✓ |

### Cross-checks against literature
- **Paper (Suleiman et al.)**: Reporting ΔT₅₀ offsets qualitatively; our numbers (4.2°C and 0.1°C) match their figure descriptions
- **Hysteresis ΔH**: Our values (5.2°C c-plane, 4.3°C r-plane) within ~5% of paper tables → digitization good

---

## 8. Limitations & Future Improvements

### Known limitations
1. **Digitized data**: Not raw; introduces ~5–10% quantitative uncertainty
2. **Single material**: Only VO₂ validated; needs another MIT system (e.g., TiO₂, NdNiO₃)
3. **No theory**: $\mathcal{C}_{AB}$ is empirical; lacks statistical mechanics grounding
4. **Static window**: [40°C, 90°C] hardcoded for VO₂; different materials may need different windows

### Potential improvements
- Replace Pearson with rank-order correlation (Spearman ρ) → more robust to outliers
- Implement adaptive window selection based on transition onset/end
- Add uncertainty propagation (Gaussian error model on digitized coordinates)
- Compare $\mathcal{C}_{AB}$ to alternative measures (mutual information, curvature correlation)

---

## 9. References

**Data source**:
- Suleiman, M., Mansouri, S., et al. (2021). "Unraveling the metal-insulator transition in VO₂ through direct measurement of spatial carrier density variations." *Scientific Reports* **11**, 1620. https://doi.org/10.1038/s41598-020-79758-1

**Methods & inspiration**:
- Pearson correlation in statistical analysis (Fisher, 1915)
- Gaussian filtering in signal processing (Savitzky–Golay alternative: 1964)
- Transition temperature definitions in physics of phase transitions (standard T50 midpoint, widely used in MIT literature)

---

## 10. Contact & Contributing

This is active research. If you want to:
- **Report a bug**: Check reproducibility against [scripts/](scripts/) and open an issue
- **Propose new data**: Send digitized (or raw) data from another MIT system; we can run the analysis
- **Suggest improvements**: Discuss in issues or via pull requests

---

*Last updated: April 2026*  
*Repo status: Active exploration — not a finished product*
