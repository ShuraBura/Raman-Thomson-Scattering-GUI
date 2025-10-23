# RSTSGUI2 - Raman & Thomson Scattering Analysis User Manual

**Version 9.9**
**Date: October 2024**

---

## Table of Contents
1. [Overview](#overview)
2. [Installation](#installation)
3. [Quick Start](#quick-start)
4. [GUI Layout](#gui-layout)
5. [Theory & Equations](#theory--equations)
6. [Fitting Methods](#fitting-methods)
7. [Workflow Guide](#workflow-guide)
8. [Parameter Tuning](#parameter-tuning)
9. [Troubleshooting](#troubleshooting)
10. [References](#references)

---

## Overview

RSTSGUI2 is a MATLAB-based graphical user interface for analyzing Raman Scattering (RS) and Thomson Scattering (TS) spectroscopy data from plasma diagnostics. The tool provides:

- **Automated intensity calibration** using Raman scattering from known gases
- **Multiple fitting algorithms** for electron temperature (Te) and density (ne)
- **Realistic error estimation** using Monte Carlo sampling
- **Interactive parameter adjustment** and visualization

### Key Features (V9.9)

- ✅ Weighted fitting prioritizes central spectral region (532 nm)
- ✅ Monte Carlo error estimation (~15-20% realistic errors)
- ✅ Three independent fitting methods for cross-validation
- ✅ Real-time diagnostic output and parameter adjustment

---

## Installation

### Requirements

- **MATLAB R2018b or newer** (requires `uifigure` support)
- **Image Processing Toolbox**
- **Optimization Toolbox**

### Setup

1. **Download RSTSGUI2.m** from GitHub:
   ```
   https://github.com/ShuraBura/Raman-Thomson-Scattering-GUI
   Branch: claude/matlab-support-011CUNPs1aDJHFFbyb5mLPrh
   ```

2. **Place in MATLAB path** or navigate to the directory:
   ```matlab
   cd /path/to/Raman-Thomson-Scattering-GUI
   ```

3. **Edit default file paths** (lines 26-107) to point to your data

4. **Run the GUI**:
   ```matlab
   analyze_RSTS
   ```

---

## Quick Start

### Typical Analysis Workflow

1. **Load data** → GUI opens with initial Gaussian fit displayed
2. **Adjust ROI/Exclusion** → Set signal region and laser line exclusion
3. **Recalibrate Raman** → Optimize wavelength calibration (optional)
4. **Choose fitting method**:
   - **Int.Calib.Full** (recommended) → Best absolute ne
   - **Shape Calib.** → Verification/consistency check
   - **Int.Calib.Gauss** → Fast, simplified analysis
5. **Review results** → Check fit quality (R², SSR, visual inspection)
6. **Save outputs** → Export CSV data and plots

---

## GUI Layout

### Display Panels (Top)

```
┌─────────────────────────────────────────────────────┐
│  Raman Image       │  Thomson Image  │  File List   │
├─────────────────────────────────────────────────────┤
│  Raman Spectrum    │  Thomson Spectrum Fit           │
└─────────────────────────────────────────────────────┘
```

- **Raman Image**: Raw Raman scattering image with ROI overlay
- **Thomson Image**: Thomson scattering image with ROI overlay
- **Raman Spectrum**: Experimental vs. theoretical Raman spectrum
- **Thomson Spectrum**: Fitted spectrum with current parameters
- **File List**: Browse Thomson files in current directory

### Control Panel (Bottom)

#### Row 1: Image Parameters
- **ROI yc, wy**: Vertical center and width of signal integration region
- **Excl: start to end**: Pixel columns to exclude (laser line)
- **Refresh**: Update images and ROI overlays

#### Row 2: Acquisition Parameters
- **RS Shots, TS Shots**: Number of laser shots averaged
- **a, b**: Wavelength calibration (λ = a·pixel + b)
- **Toggle TS sym**: Switch between line and symbol plot styles

#### Row 3: Energy & Calibration
- **RS E(mJ), TS E(mJ)**: Laser energy per shot (millijoules)
- **Linewidth**: Raman spectral linewidth (nm)
- **P(Pa)**: Gas pressure (Pascals)

#### Row 4: Gas Composition
- **T(K)**: Gas temperature (Kelvin)
- **N2, O2**: Gas mixing ratios (must sum to 1.0)
- **Recalib disp. from Raman**: Re-optimize wavelength calibration

#### Row 5: Filters
- **RmOut**: Remove outliers using IQR method
- **MedF**: Apply median filter to spectrum
- **Window**: Median filter window size
- **Recalc Raman Spectrum**: Update Raman analysis

#### Row 6: Plasma Parameters
- **Te(eV)**: Electron temperature (electron volts)
- **ne**: Electron density (m⁻³)
- **α**: Plasma screening parameter (read-only)
- **γ**: Duty cycle factor (for pulsed plasmas)

#### Row 7: Fitting Controls
- **Multi**: Number of multi-start optimizations (1-20)
- **Recalc TS**: Manually update spectrum from Te/ne/α/γ
- **Shape Calib.**: Fit shape to get Te/α/ne
- **Int.Calib.Gauss**: Gaussian fit with intensity calibration
- **Int.Calib.Full**: Full coherent fit with intensity calibration
- **< / >**: Pixel shift controls for manual alignment

### Utilities Panel (Right)

- **Decouple**: Open plots in separate figure windows
- **EEDF**: Plot electron energy distribution function
- **Save CSV**: Export wavelength, experimental, and fitted data
- **Load buttons**: Select new Raman/Thomson files
- **BG buttons**: Load background images for subtraction

---

## Theory & Equations

### Raman Scattering (Intensity Calibration)

Raman scattering from known gases (N₂, O₂) provides absolute intensity calibration.

**Calibration Factor**:
```
C = I_exp / (E_total × I_theory)
```

Where:
- `I_exp` = ∫ experimental Raman spectrum dλ (measured counts)
- `E_total` = laser energy × shots (Joules)
- `I_theory` = ∫ theoretical Raman spectrum dλ (calculated from cross-sections)

**Theoretical Raman Spectrum**:
```
I_Raman(λ) = Σ_species σ_species(λ) × n_species × G(λ, FWHM)
```

- `σ_species(λ)`: Raman cross-section for each species (N₂, O₂)
- `n_species`: Number density of species
- `G(λ, FWHM)`: Instrumental broadening (Gaussian)

### Thomson Scattering (Plasma Diagnostics)

Thomson scattering measures electron temperature and density from spectral shape and intensity.

#### Electron Temperature from Spectral Width

**FWHM-Te Relation** (for Gaussian approximation):
```
Te = (m_e c²) / (32 ln2) × (FWHM / λ_laser)²
```

Where:
- `m_e c² = 0.511 MeV` (electron rest energy)
- `FWHM`: Full width at half maximum (nm)
- `λ_laser = 532 nm` (laser wavelength)

**Physical Interpretation**: Wider spectrum → higher electron thermal velocity → higher Te

#### Electron Density from Absolute Intensity

**Density Calculation**:
```
n_e = N_total / (C × E_total × σ_T)
```

Where:
- `N_total` = ∫ Thomson spectrum dλ (total scattered photons)
- `C`: Calibration factor from Raman (m⁻³ J⁻¹)
- `E_total`: Total laser energy (J)
- `σ_T = 7.94 × 10⁻³⁰ m²` (Thomson differential cross-section at 90°)

#### Plasma Screening Parameter (α)

**Definition**:
```
α = 1 / (k × λ_D)
```

Where:
- `k = (4π/λ_laser) sin(θ/2)` (scattering wave vector, θ = 90°)
- `λ_D = √(ε₀ k_B T_e / (n_e e²))` (Debye length)

**Relation to Density**:
```
n_e = (α² k² ε₀ k_B T_e) / e²
```

**Physical Interpretation**:
- `α << 1`: Incoherent scattering (individual electrons)
- `α ~ 1`: Transitional regime
- `α >> 1`: Coherent/collective scattering (plasma waves)

#### Coherent Thomson Scattering Spectrum

**Spectral Form Factor** (Salpeter approximation):
```
S_e(k, ω) = (1/√π) × exp(-x_e²) / |1 + α² W(x_e)|²
```

Where:
- `x_e = Δω / (k v_e)` (normalized frequency shift)
- `v_e = √(2 k_B T_e / m_e)` (electron thermal velocity)
- `W(x_e)` = plasma dispersion function (includes Dawson integral)
- `Δω = ω_scattered - ω_laser` (frequency shift)

**Key Features**:
- Central peak: free electron scattering
- Satellites: ion acoustic waves (if α > 1)
- Asymmetry: collective effects

---

## Fitting Methods

### Method 1: Int.Calib.Gauss (Simplest)

**Model**: Gaussian approximation to Thomson spectrum

**Fitted Parameters**:
- `A`: Amplitude
- `λ₀`: Center wavelength
- `FWHM`: Full width at half maximum
- `b`: Baseline offset

**Gaussian Model**:
```
I(λ) = A × exp(-4 ln2 × ((λ - λ₀) / FWHM)²) + b
```

**Calculations**:
1. Fit Gaussian to spectrum → get `FWHM`, `A`
2. Calculate `Te` from `FWHM` (see FWHM-Te relation)
3. Calculate total counts: `N_total = ∫ I(λ) dλ ≈ A × FWHM × √(π/(4ln2))`
4. Calculate `n_e = N_total / (C × E_total × σ_T)`

**Pros**:
- ✅ Fast computation
- ✅ Robust convergence
- ✅ Good for initial estimates

**Cons**:
- ❌ Assumes Gaussian shape (ignores collective effects)
- ❌ Less accurate for α > 0.5
- ❌ No α fitting

**When to Use**: Quick analysis, low-density plasmas (α < 0.3), initial estimates

---

### Method 2: Shape Calib. (Shape-Based)

**Model**: Full coherent Thomson scattering with plasma effects

**Fitted Parameters**:
- `Te`: Electron temperature (eV)
- `α`: Plasma screening parameter
- `A`: Amplitude (arbitrary normalization)

**Objective Function**:
```
minimize Σ [I_exp(λ) - A × S_coherent(λ, Te, α)]²
```

Where `S_coherent(λ, Te, α)` is calculated from Salpeter theory.

**Calculations**:
1. Fit (Te, α, A) to match spectral shape
2. Calculate `n_e` from `α` and `Te` using dispersion relation
3. **Does NOT use** absolute intensity calibration

**Pros**:
- ✅ Accounts for collective scattering effects
- ✅ Independent verification of ne
- ✅ Sensitive to plasma physics

**Cons**:
- ❌ ne depends on α fit (not direct intensity measurement)
- ❌ More sensitive to spectral shape errors
- ❌ Can differ from intensity-based ne

**When to Use**: Verification, consistency checks, high-α plasmas

**Diagnostic Output**:
```
NOTE - Shape Calib uses ne from alpha (shape only), not intensity:
  ne from alpha (used) = 2.45e+17 m^-3
  ne from intensity = 5.82e+17 m^-3 (ratio: 2.38 x)
```

This shows how ne from shape compares to what intensity calibration would give.

---

### Method 3: Int.Calib.Full (Recommended)

**Model**: Full coherent Thomson scattering with absolute intensity

**Fitted Parameters**:
- `Te`: Electron temperature (eV)
- `α`: Plasma screening parameter
- `A`: Amplitude (physical units)

**Two-Stage Process**:

**Stage 1: Shape Fitting**
```
minimize Σ w(λ) × [I_exp(λ) - A × S_coherent(λ, Te, α)]²
```

Where `w(λ)` = Gaussian weights centered at 532 nm:
```
w(λ) = exp(-((λ - 532) / σ)²),  σ = 5 nm
```

**Why Weighted**: Prioritizes central region where signal is strongest, prevents fits being pulled by noisy wings.

**Stage 2: Intensity Calculation**
1. Calculate total counts: `N_total = ∫ A × S_coherent(λ, Te, α) dλ`
2. Calculate `n_e = N_total / (C × E_total × σ_T)`
3. **Independent of α** - uses direct photon counting

**Pros**:
- ✅ Best absolute ne accuracy (uses Raman calibration)
- ✅ Accounts for collective effects in shape
- ✅ Weighted fitting improves central peak fit
- ✅ Two independent measurements (shape → Te, intensity → ne)

**Cons**:
- ❌ Slowest method (~10-30 seconds)
- ❌ Requires accurate Raman calibration
- ❌ Most sensitive to calibration errors

**When to Use**: Final quantitative analysis, publication-quality results

---

## Error Estimation (Monte Carlo Method)

All three fitting methods use **Monte Carlo sampling** to estimate realistic parameter uncertainties.

### Traditional Method (Jacobian - Too Small)

```
σ_param = √(diag(σ_residual² × (J^T J)^-1))
```

**Problem**: Only captures statistical noise, ignores model uncertainty → errors 10-100× too small!

### Monte Carlo Method (V9.9 - Realistic)

**Algorithm**:

1. **Perform best fit** → get best_Te, best_ne, best_SSR

2. **Define acceptance threshold**:
   ```
   SSR_threshold = 1.5 × best_SSR
   ```
   (Adjust `1.5` for more/less conservative errors)

3. **Grid search parameter space**:
   ```
   Te_range = [0.7 × best_Te, 1.3 × best_Te] (25 points)
   ne_range = [0.7 × best_ne, 1.3 × best_ne] (25 points)
   ```

4. **For each (Te_test, ne_test)**:
   - Generate synthetic spectrum
   - Calculate SSR_test
   - If `SSR_test < SSR_threshold` → accept parameters

5. **Calculate errors**:
   ```
   Te_error = (max(Te_accepted) - min(Te_accepted)) / 2
   ne_error = (max(ne_accepted) - min(ne_accepted)) / 2
   ```

**Physical Interpretation**:
- Error = range of parameters producing "visually indistinguishable" fits
- Captures both statistical noise AND model uncertainty
- Typically gives 15-20% errors (matches empirical observation)

**Console Output**:
```
Calculating realistic parameter errors via Monte Carlo sampling...
  SSR threshold for acceptable fits: 3.51e+04 (best fit SSR = 2.34e+04)
  Monte Carlo results: 189 parameter sets within acceptable SSR
  Te range: 1.21 to 1.87 eV (error = 0.33 eV, 21.7%)
  ne range: 1.45e+17 to 3.21e+17 m^-3 (error = 8.80e+16 m^-3, 19.4%)
```

### Tuning Error Estimation

**Line 1259 (Int.Calib.Full), 1407 (Shape Calib.), 1621 (Int.Calib.Gauss)**:
```matlab
acceptable_ssr_threshold = best_resnorm * 1.5;
```

- **Increase multiplier** (e.g., 2.0) → larger errors (more conservative)
- **Decrease multiplier** (e.g., 1.3) → smaller errors (tighter bounds)
- **Recommended**: 1.3-1.7 depending on data quality

---

## Workflow Guide

### Standard Analysis Procedure

#### 1. Data Preparation
- Ensure Raman and Thomson images are in same format (TIFF, PNG)
- Record acquisition parameters:
  - Laser shots
  - Laser energy per shot
  - Gas pressure and composition
- Obtain background images (laser off) if available

#### 2. Initial Setup
```matlab
% Edit RSTSGUI2.m lines 26-107
default_set = 2;  % Choose your dataset

% Update file paths
p.raman_file = 'path/to/raman.tif';
p.thomson_file = 'path/to/thomson.tif';
p.raman_bg_file = 'path/to/raman_bg.tif';  % optional
p.thomson_bg_file = 'path/to/thomson_bg.tif';  % optional

% Set acquisition parameters
p.shots_raman = 300;
p.shots_thomson = 300;
p.energy_J_raman = 52e-3;  % 52 mJ
p.energy_J_thomson = 52e-3;
```

#### 3. ROI Adjustment
- **ROI yc**: Center of signal vertically (pixel row)
- **ROI wy**: Width of signal vertically (pixels)
- **Excl start-end**: Columns to exclude (laser line, typically 455-550)

**Tips**:
- Look at raw images to identify signal region
- Make ROI wide enough to capture all signal (typically 30-60 pixels)
- Exclude enough around laser line to avoid saturation

#### 4. Wavelength Calibration

**Initial Calibration** (if known):
```matlab
p.a_initial = 0.010050;  % nm/pixel
p.b_initial = 526.7384;   % nm offset
```

**Optimization**:
- Click **"Recalib disp. from Raman"**
- Uses peak matching between experimental and theoretical Raman
- Updates `a` and `b` coefficients
- Check that Raman peaks align in spectrum plot

**Manual Fine-Tuning**:
- Adjust `a` and `b` fields directly
- Click **"Recalc Raman Spectrum"** to update
- Iterate until theoretical/experimental peaks align

#### 5. Raman Analysis Verification

Check Raman spectrum plot:
- **Good calibration**: Blue (experimental) and red (model) peaks aligned
- **Calibration factor**: Should be ~10⁻¹⁵ to 10⁻¹⁴ (order of magnitude)
- **If mismatch**: Re-run calibration or adjust gas parameters

#### 6. Thomson Fitting

**Option A: Quick Analysis**
- Use initial Gaussian fit (displayed on file load)
- Click **"Int.Calib.Gauss"** for Monte Carlo errors
- Fastest, good for surveys

**Option B: Accurate Analysis** (Recommended)
- Click **"Int.Calib.Full"**
- Wait ~10-30 seconds for weighted fit + Monte Carlo
- Check console output for diagnostics
- Best absolute accuracy

**Option C: Verification**
- Click **"Shape Calib."**
- Compare ne from shape vs. intensity
- Ratio should be 1-3× typically
- Large discrepancy (>5×) indicates calibration issues

#### 7. Result Validation

**Check Fit Quality**:
- **R²**: Should be > 0.95 for good fits
- **SSR**: Lower is better (compare between methods)
- **Visual**: Fitted curve should match data well

**Check Errors**:
- Should be 10-25% typically
- Much smaller (<5%) → suspicious, check MC sampling
- Much larger (>30%) → noisy data or poor fit

**Cross-Validation**:
- Compare Te between methods → should agree within errors
- Compare ne from Int.Calib.Full vs. Shape Calib. → typically 1-3× ratio

#### 8. Parameter Sweeps

To analyze multiple files:

```matlab
% In GUI file list (right panel)
% 1. Click on file name
% 2. Wait for automatic reanalysis
% 3. Click fitting button of choice
% 4. Save results

% Or: loop programmatically
files = dir('*.tif');
for i = 1:length(files)
    % Load file, fit, save results
end
```

#### 9. Data Export

**CSV Export**:
- Click **"Save CSV"**
- Creates file: `TS_<method>_Te<value>_ne<value>_SSR<value>.csv`
- Contains: wavelength, experimental intensity, fitted intensity

**Plot Export**:
- Click **"Save TS Spectrum"** (or Raman/Image buttons)
- Creates JPG files with embedded parameters in filename

**Manual Export**:
```matlab
% Access from workspace after GUI runs:
handles = guidata(gcf);
Te = handles.edit_Te_manual.Value;
ne = handles.edit_ne_manual.Value;
wavelength = handles.plot_data.wavelength_nm;
spectrum = handles.current_ts_spectrum;
```

---

## Parameter Tuning

### Optimization Settings

**Multi-start Fitting** (Line 1196, 1363, 1562):
```matlab
num_starts = handles.edit_multistarts.Value;  % Default: 6
```

- **Purpose**: Avoid local minima in optimization
- **Typical**: 5-10 starts
- **More starts**: More robust, but slower
- **Fewer starts**: Faster, risk missing global minimum

**GUI Control**: "Multi" field (row 7)

---

### Weighted Fitting (Int.Calib.Full)

**Central Weighting** (Line 1176):
```matlab
central_wl = 532;     % Center wavelength (nm)
weight_sigma = 5;     % Decay width (nm)
```

**Weight Function**:
```
w(λ) = exp(-((λ - 532)² / (2 × 5²)))
```

**Effect**:
- `σ = 3 nm`: Very tight focus on center
- `σ = 5 nm`: Moderate (default)
- `σ = 8 nm`: Broader, includes more wings

**When to Adjust**:
- **Narrow peaks** (cold plasma) → decrease σ
- **Broad peaks** (hot plasma) → increase σ
- **Noisy wings** → decrease σ (ignore wings)

---

### Monte Carlo Sampling

**Acceptance Threshold** (Lines 1259, 1407, 1621):
```matlab
acceptable_ssr_threshold = best_resnorm * 1.5;
```

**Effect**:
- `1.3`: Tighter errors (~10-15%)
- `1.5`: Moderate errors (~15-20%, default)
- `2.0`: Looser errors (~20-30%)

**When to Adjust**:
- **High-quality data** (low noise) → 1.3
- **Noisy data** → 2.0
- **Publication** → 1.5 (standard)

**Grid Resolution** (Lines 1254, 1405, 1619):
```matlab
Te_range = linspace(Te_fit * 0.7, Te_fit * 1.3, 25);
ne_range = linspace(ne_fit * 0.7, ne_fit * 1.3, 25);
```

- **25 points**: Default (625 combinations)
- **Increase** (e.g., 35): Finer resolution, slower
- **Decrease** (e.g., 15): Faster, coarser errors

---

### Gas Parameters

**Pressure** (Line 474):
```matlab
handles.edit_pressure.Value = 101325;  % Pa (1 atm)
```

**Temperature** (Line 479):
```matlab
handles.edit_temp_K.Value = 293.15;  % K (20°C)
```

**Mixing Ratios** (Lines 481, 483):
```matlab
N2_fraction = 0.79;  % 79% nitrogen
O2_fraction = 0.21;  % 21% oxygen
```

**Important**: N2 + O2 must equal 1.0 (automatically normalized)

---

## Troubleshooting

### Poor Fit Quality (Low R²)

**Symptoms**: R² < 0.90, large SSR, visual mismatch

**Possible Causes**:
1. **Incorrect wavelength calibration**
   - Solution: Run "Recalib disp. from Raman"
   - Verify Raman peaks align

2. **Wrong ROI/exclusion**
   - Solution: Adjust yc, wy, exclusion columns
   - Ensure signal fully captured, laser line excluded

3. **Stray light/noise**
   - Solution: Use background subtraction
   - Enable outlier removal (RmOut checkbox)
   - Apply median filter (MedF, window=5)

4. **Incorrect acquisition parameters**
   - Solution: Verify shots, energy, pressure, gas composition
   - Check that Raman calibration factor is reasonable (~10⁻¹⁵)

---

### Very Small Errors (<5%)

**Symptoms**: Monte Carlo reports <5% errors

**Possible Causes**:
1. **Acceptance threshold too tight**
   - Solution: Increase multiplier (1.5 → 2.0)

2. **Few accepted parameter sets**
   - Check console: should be 100-300 sets
   - If <50 sets: widen search range or threshold

3. **Perfect fit** (unusual)
   - May be real if data quality is exceptional
   - Verify by visual inspection

---

### Large Errors (>30%)

**Symptoms**: Monte Carlo reports >30% errors

**Possible Causes**:
1. **Noisy data**
   - Solution: Average more shots
   - Use background subtraction
   - Enable filters

2. **Poor fit convergence**
   - Increase multi-start attempts
   - Check initial parameter guesses

3. **Model mismatch**
   - Check α value (should be 0.1-2.0 typically)
   - Verify plasma assumptions valid

---

### Density Discrepancies Between Methods

**Symptoms**: Int.Calib.Full and Shape Calib. give very different ne

**Expected**: Ratio of 1-3× is normal

**If Ratio > 5×**:
1. **Check Raman calibration**
   - Verify gas pressure, temperature, composition
   - Re-run "Recalib disp. from Raman"

2. **Check Thomson acquisition parameters**
   - Verify shots, energy are correct
   - Ensure no saturation in images

3. **Check spectral shape**
   - If α is very high (>2), collective effects may be strong
   - Shape-based ne may be less reliable

**Recommendation**: Trust Int.Calib.Full for absolute ne (uses direct intensity measurement)

---

### Fit Converges to Wrong Shape

**Symptoms**: Fitted spectrum visibly misaligned with data

**Possible Causes**:
1. **Local minimum in optimization**
   - Solution: Increase multi-start attempts (6 → 12)

2. **Bad initial guess**
   - Manually set Te, ne closer to expected values
   - Click "Recalc TS" to update

3. **Wavelength calibration off**
   - Use pixel shift controls (< / >) to align
   - Re-calibrate from Raman

---

### GUI Crashes or Errors

**Common Issues**:

1. **"File not found"**
   - Check file paths in lines 26-107
   - Use absolute paths

2. **"Undefined function"**
   - Missing toolbox (Optimization, Image Processing)
   - Install required toolboxes

3. **"Out of memory"**
   - Reduce image size
   - Close other MATLAB processes

4. **"Invalid parameter"**
   - Check that N2 + O2 = 1.0
   - Verify all fields have valid numbers

---

## Best Practices

### Data Acquisition

- ✅ Use consistent laser energy (minimize shot-to-shot variation)
- ✅ Average 100-500 shots for good statistics
- ✅ Acquire background with laser off (same integration time)
- ✅ Use reference Raman scattering from known gas at known pressure
- ✅ Minimize stray light (proper baffling, alignment)

### Analysis Workflow

- ✅ Start with Int.Calib.Gauss for quick initial estimates
- ✅ Use Int.Calib.Full for final quantitative results
- ✅ Verify with Shape Calib. for consistency
- ✅ Check that R² > 0.95 for all methods
- ✅ Compare errors across methods (should be similar)

### Parameter Selection

- ✅ Use multi-start = 6-10 for robust fitting
- ✅ Use acceptance threshold = 1.5 for standard errors
- ✅ Adjust weighted fitting σ based on spectral width
- ✅ Re-calibrate wavelength if peaks don't align

### Reporting Results

- ✅ Report method used (Int.Calib.Full recommended)
- ✅ Include Monte Carlo errors (±15-20% typical)
- ✅ State acquisition conditions (shots, energy, gas)
- ✅ Report fit quality (R², SSR)
- ✅ Show fitted spectrum vs. data (visual verification)

---

## Advanced Features

### Custom Fitting Functions

To modify fitting behavior, edit these functions:

**`coherent_fitter_direct`** (Line 1487):
- Calculates coherent Thomson spectrum
- Modify to change shape model

**`calculate_coherent_ts`** (Line 1959):
- Core Salpeter formula
- Includes plasma dispersion function

**`intensity_calibration_callback`** (Line 1163):
- Main fitting routine for Int.Calib.Full
- Modify weighting, optimization

### Batch Processing

Process multiple files programmatically:

```matlab
% Define file list
files = {
    'file1.tif',
    'file2.tif',
    'file3.tif'
};

% Initialize results storage
results = struct('Te', [], 'ne', [], 'Te_err', [], 'ne_err', []);

% Loop through files
for i = 1:length(files)
    % Set Thomson file
    params.thomson_file = files{i};

    % Run analysis
    [res, data] = run_full_analysis(params);

    % Store results
    results.Te(i) = res.corrected.Te_eV;
    results.ne(i) = res.corrected.ne_m3;
    results.Te_err(i) = res.deconv_results.Te_eV_error;
    results.ne_err(i) = res.deconv_results.ne_error;
end

% Plot results
errorbar(1:length(files), results.Te, results.Te_err);
xlabel('File Index');
ylabel('Te (eV)');
```

### Exporting Custom Data

Access internal data structures:

```matlab
% After GUI is open
handles = guidata(gcf);

% Get all data
wavelength = handles.plot_data.wavelength_nm;
spectrum = handles.current_ts_spectrum;
fitted = get(handles.h_coherent_fit, 'YData');
Te = handles.edit_Te_manual.Value;
ne = handles.edit_ne_manual.Value;
alpha = handles.edit_alpha.Value;

% Export
save('my_analysis.mat', 'wavelength', 'spectrum', 'fitted', 'Te', 'ne', 'alpha');
```

---

## References

### Theory

1. **Sheffield, J.** (1975). *Plasma Scattering of Electromagnetic Radiation*. Academic Press.
   - Comprehensive treatment of Thomson scattering theory

2. **Salpeter, E.E.** (1960). "Electron Density Fluctuations in a Plasma", *Phys. Rev.* **120**, 1528.
   - Original plasma dispersion function formulation

3. **Hutchinson, I.H.** (2002). *Principles of Plasma Diagnostics* (2nd ed.). Cambridge University Press.
   - Modern plasma diagnostic techniques

### Experimental

4. **Muraoka, K. & Maeda, M.** (1991). "Laser-Aided Diagnostics of Plasmas and Gases", *IOP Publishing*.
   - Practical laser diagnostic methods

5. **Evans, D.E. & Katzenstein, J.** (1969). "Laser light scattering in laboratory plasmas", *Rep. Prog. Phys.* **32**, 207.
   - Classic review of Thomson scattering diagnostics

### Raman Scattering

6. **Long, D.A.** (2002). *The Raman Effect: A Unified Treatment of the Theory of Raman Scattering*. Wiley.
   - Comprehensive Raman theory and applications

### Error Analysis

7. **Press, W.H. et al.** (2007). *Numerical Recipes* (3rd ed.). Cambridge University Press.
   - Monte Carlo methods and parameter estimation

---

## Version History

### V9.9 (October 2024) - RSTSGUI2
- Added weighted fitting prioritizing 532 nm central region
- Implemented Monte Carlo error estimation for all methods
- Fixed initial Gaussian fit ne calculation bug
- Added diagnostic output for Shape vs. Intensity ne comparison
- Unified error methodology across all fitting callbacks

### V9.8 (Prior)
- Automatic Te, ne, alpha field updates after Gaussian fit
- Editable gamma field for duty cycle
- Multi-start optimization for Raman calibration
- Status message display

### Earlier Versions
- Basic GUI framework
- Gaussian and coherent fitting
- Raman calibration
- ROI and exclusion controls

---

## Support

For issues, questions, or contributions:

**GitHub**: https://github.com/ShuraBura/Raman-Thomson-Scattering-GUI
**Branch**: claude/matlab-support-011CUNPs1aDJHFFbyb5mLPrh

**Reporting Bugs**:
1. Provide MATLAB version
2. Describe steps to reproduce
3. Include error messages
4. Attach sample data if possible

---

## License

[Specify license here - e.g., MIT, GPL, etc.]

---

## Acknowledgments

Development supported by [Institution/Grant].

Monte Carlo error estimation and weighted fitting improvements implemented with assistance from Claude Code (Anthropic).

---

**Document Version**: 1.0
**Last Updated**: October 2024
**Software Version**: RSTSGUI2 V9.9
