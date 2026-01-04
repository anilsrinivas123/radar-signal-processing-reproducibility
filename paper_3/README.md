# FFT-Based Three-Point Frequency Estimation

This repository contains a methodological reproduction and comparative implementation of
FFT-based single-tone frequency estimation algorithms described in:

> Eric Jacobsen and Peter Kootsookos,  
> *Fast, Accurate Frequency Estimators*,  
> IEEE Signal Processing Magazine, May 2007.

The work focuses on computationally efficient three-point interpolation techniques that refine
coarse FFT peak estimates to achieve sub-bin frequency accuracy without increasing FFT size.

---

## 📌 Project Overview

Frequency estimation using the FFT is limited by bin resolution. The referenced paper shows that
using only three FFT samples around the spectral peak allows accurate estimation of the true
frequency with minimal computational overhead.

This project:
- Implements all core estimators presented in the paper
- Provides a unified comparison framework
- Includes supporting utilities for efficient coefficient computation
- Explores practical extensions beyond the original paper’s scope

---

## 🔬 Implemented Estimators (From the Paper)

All estimators operate on three FFT samples:  
\( X_{k-1}, X_k, X_{k+1} \), where \( k \) is the peak FFT bin.

### Core Estimators
- **Parabolic Interpolation (Eq. 2)**  
  Simple magnitude-based three-point estimator (baseline, biased).

- **Unbiased Complex FFT Estimator (Eq. 3)**  
  Uses complex FFT values to eliminate bias and improve accuracy for unwindowed data.

- **Window-Adaptive Estimators (Eqs. 4–5)**  
  Incorporate window-specific scaling constants for non-rectangular windowing.

- **Candan’s Three-Point Estimator**  
  A refined algebraic three-point method from related literature, included for comparison.

---

## 📂 Repository Structure

### Core Reproduction
- `FFT_coarse_main.m`  
  Coarse FFT peak detection.

- `freq_est_fft_parabolic.m`  
  Parabolic interpolation estimator (Eq. 2).

- `freq_est_fft_improved.m`  
  Unbiased complex FFT estimator (Eq. 3).

- `freq_est_methods_three_points.m`  
  Unified comparison framework for multiple estimators.

- `candan_three_point_freq.m`  
  Implementation of Candan’s estimator.

---

### Supporting Utilities
- `get_coefs2delta_real_2N.m`  
  Precompute real-valued coefficients for efficient estimation.

- `get_coefs2delta_imag_2N.m`  
  Precompute imaginary-valued coefficients.

- `inverse_mappings.m`  
  Convert fractional bin offset \( (k, \delta) \) to physical frequency.

---

### Experiment Drivers
- `Frequency_estimators_main.m`  
  Main evaluation script for estimator comparison.

- `candan_freq_est_main.m`  
  Focused evaluation of Candan’s estimator.

---

## 🚀 Extensions Beyond the Original Paper

The following scripts extend the original work and are **clearly separated** from reproduction:

- `freq_est_cosine_2N.m`  
  Cosine-model-based high-resolution estimator using doubled FFT length (2N).

- `FMCW.m`  
  Application of frequency estimators to FMCW radar signals, demonstrating the impact of
  sub-bin frequency accuracy on range estimation.

These are included as **exploratory and application-level extensions**, not as original
contributions of the referenced paper.

---

## ▶️ How to Run

1. Open MATLAB and set the repository root as the working directory.
2. Run:
```matlab
Frequency_estimators_main
