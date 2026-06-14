# FluSight: Recursive Quantile Regression Forecasting for Influenza Hospitalizations

## Overview

This repository contains a **recursive hybrid forecasting model** for predicting weekly influenza-related inpatient hospitalizations in South Carolina. The model applies **1-step quantile regression recursively across multiple horizons**, combining **direct quantile regression with phase-aligned analog correction**.


## Key Features

### 1. Recursive Forecasting Architecture 
  - Fit H0 model (1-step ahead quantile regression)
  - Append H0 prediction to feature series and recompute all features (slopes, lags, rolling stats, phase)
  - Use updated features to forecast H1 with the same H0 model
  - Repeat for H2, H3

### 2. Hybrid Forecasting Approach
- **Quantile Regression (QR)**: 1-step horizon-specific quantile regression model
- **Analog Correction**: Identifies historical periods with similar epidemiological states and uses their outcomes to adjust predictions
- **Decline-Aware Weighting**: Applies phase-specific penalties during decline periods to reduce overly optimistic rebound predictions

### 3. Feature Engineering
The model uses engineered features including:
- **Lag features**: Target values at 1, 2, 3, 4, 5, 8, and 12-week lags
- **Slope features**: 1-4 week slopes and slope changes to capture trajectory
- **Acceleration features**: Rate of change in slopes (2-week and 3-week acceleration)
- **Decline persistence**: Indicators for sustained negative trends
- **Rolling statistics**: 4-, 8-, and 12-week means and standard deviations
- **Position in cycle**: Where current levels sit in the 12-week rolling range
- **Peak dynamics**: Weeks since 12-week peak, distance from peak, and post-peak decline speed
- **EHR indicators**: Hospitalization and positive test slopes with decline confirmation
- **Seasonal features**: Week of year harmonics (sin/cos) and holiday indicators

### 4. Phase-Based Classification
The model identifies and adapts to 9 epidemic phases:
1. **FAST_POST_PEAK_DECLINE** - Rapid post-peak decline
2. **POST_PEAK_DECLINE** - Standard post-peak decline
3. **ACCELERATING_DECLINE** - Declining trend with acceleration
4. **SLOWING_DECLINE** - Declining trend with deceleration
5. **ACCELERATING_INCREASE** - Rising trend with acceleration
6. **SLOWING_INCREASE** - Rising trend with deceleration
7. **NEAR_PEAK_DECELERATION** - Near peak with decelerating trend
8. **VOLATILE_PLATEAU** - Stable levels with high variability
9. **STABLE_PLATEAU** - Stable levels with low variability

Phase detection drives adaptive interval scaling and QR floor constraints.

### 5. Data Integration
The model integrates multiple data sources:
- **CDC NHSN Data**: Official weekly flu hospitalizations (location 45 = South Carolina)
- **Electronic Health Record (EHR) Data**: 
  - Prisma Health weekly influenza surveillance data


## Repository Contents

### Files

- **`QR.r`** - Main forecasting pipeline 
  - Data loading and merging
  - Feature engineering
  - Phase labeling with origin-specific thresholds (no temporal leakage)
  - **Tuning for H0 only** (single-step quantile regression)
  - **Recursive forecast generation**:
    - Fit H0 → predict → append predicted row → recompute all features
    - Use updated features to forecast H1 (reusing H0 model) → predict → append → recompute
    - Repeat for H2, H3
  - Quantile interval construction (blended analog + residual errors)
  - Diagnostic outputs and visualization
  - CDC submission file generation
  - Software implementation and evaluation file outputs

- **`DMAPRIME-QR.yml`** - Model metadata and documentation
  - Team and contact information
  - Model description and methods
  - Data source specification
  - Funding and licensing information

- **`README.md`** - This file

## Workflow

### 1. **Data Preparation**
- Load CDC and Prisma weekly hospitalization data
- Align datasets by week
- Standardize column names

### 2. **Feature Engineering**
- Compute predictive features from historical time series
- Apply log transformation with safe handling of zero/missing values
- Calculate rolling statistics and seasonal components
- Determine epidemic phase using origin-specific thresholds

### 3. **Rolling Validation (H0 Only)**
- For each week in the validation period (Sept 2022 – Oct 2024):
  - Compute origin-specific thresholds from data ≤ origin_week (no leakage)
  - Train H0 quantile regression model on all prior data using complete cases
  - Identify analogs from historical pool
  - Evaluate tuning parameters:
    - `k`: Number of analogs (candidate values: 8, 12, 16, 20, 25, 30)
    - `bandwidth`: Analog weight decay (0.50–2.00)
    - `blend_qr`: Mixing proportion QR vs. analog (0–1.0 in 0.1 increments)
  - Compute Mean Absolute Error (MAE) in log space
  - Select best hyperparameters for H0

### 4. **Hyperparameter Propagation (H0 → H1, H2, H3)**
- H1, H2, H3 use the same H0 quantile regression model
- Tuning parameters (k, bandwidth) from H0 are reused across all horizons
- `blend_qr` is slightly reduced for longer horizons to encourage analog influence (penalty: h × 0.05)

### 5. **Model Training**
- Fit **single quantile regression model (H0)** using all training data up to label cutoff
- Capture residuals by phase for uncertainty quantification
- This model is reused recursively for all horizons

### 6. **Recursive Forecast Generation**
- **H0 (nowcast)**: Fit H0 QR on current features → predict 1-step ahead
- **Append H0 prediction**: Create synthetic row with predicted log_target, preserve EHR state
- **Recompute Features**: Recalculate all lag, slope, rolling stats, phase features for the new synthetic row using extended time series
- **H1 (1-wk ahead)**: Apply H0 QR model to the updated feature row → generates H1 prediction
- **Repeat for H2, H3**: Append predictions, recompute features iteratively
- **Blend predictions**: Combine QR median and analog quantiles using phase-specific weights
- **Generate quantiles**: Predictions at 23 quantile levels (0.01–0.99)
- **Apply phase-conditioned interval scaling**: Narrow intervals during decline to reduce uncertainty
- **Ensure monotonicity**: via isotonic regression across quantiles

### 7. **Output Generation**
- **CDC Submission File** (Section 12):
  - Format: CDC FluSight Hub submission format
  - Contains all quantile predictions (0–4 weeks ahead)
  
- **Software Implementation File** (Section 13):
  - Combines historical data (RFA + EHR sources) with forecasts
  - Ready for deployment in forecasting infrastructure
  
- **Evaluation File** (Section 16):
  - Training data and median-only forecasts
  - Used for model performance assessment

### 8. **Visualization**
- Plots observed time series (20 weeks back)
- Overlays median forecast and 50% / 95% prediction intervals
- Marks forecast origin line
- Indicates detected phase in title


## Key Hyperparameters

| Parameter | Default Range | Purpose |
|-----------|---------------|---------|
| `candidate_k` | 8, 12, 16, 20, 25, 30 | Number of analog neighbors (H0 tuning only) |
| `candidate_bandwidth` | 0.50–2.00 | Decay rate for analog weights (H0 tuning only) |
| `candidate_blend_qr` | 0.0–1.0 (step 0.1) | QR-analog mixing proportion (H0 tuning only) |
| `quantiles` | 0.01–0.99 (23 levels) | Forecast quantile levels (all horizons) |
| `horizons` | 0, 1, 2, 3 | Weeks ahead to forecast (recursive application of H0 model) |

## Tuning Metrics

- **Metric**: Mean Absolute Error (MAE) in log-transformed space
- **Tuning Period**: Sept 1, 2022 – Oct 1, 2024
- **Candidates Evaluated**: 6 × 6 = **36 configurations** per origin (H0 only; hyperparameters propagated to H1–H3)
- **Selection**: Best MAE for H0; same (k, bandwidth) reused for H1–H3 with horizon-indexed blend_qr adjustment

## Decline-Aware Adjustments

During decline phases:
- **QR Ceiling**: Maximum blend weight for QR (e.g., 30% for FAST_POST_PEAK_DECLINE), allowing analogs more influence
- **Rebound Penalty**: Reduce weights of analogs predicting increases by 85%
- **Weak Decline Penalty**: Reduce weights of analogs with marginal declines
- **Interval Scaling**: Narrow intervals during decline to reduce uncertainty
- **Late Peak Penalty**: During H2+ forecasts in peak phases, penalize analogs predicting peaks >2 weeks away

## Requirements

### R Packages
```r
library(dplyr)          # Data manipulation
library(lubridate)      # Date handling
library(ggplot2)        # Visualization
library(purrr)          # Functional programming
library(readr)          # CSV I/O
library(quantreg)       # Quantile regression
library(tidyr)          # Tidy data
library(slider)         # Rolling windows
library(scales)         # Plot scaling
```

### Data Files
- RFA weekly influenza data (2017–2025)
- CDC hospital admissions (NHSN)
- Prisma Health weekly influenza surveillance

**Note**: Data paths are currently hardcoded to local Box Cloud Storage directories. Modify paths before running.

## Workflow Diagram

```
Data Loading
    ↓
Feature Engineering 
    ↓
Phase Classification (origin-specific thresholds)
    ↓
Rolling Validation (H0 only)
    ↓
Hyperparameter Selection (best MAE for H0)
    ↓
Final Model Training (H0 QR)
    ↓
Recursive Forecast Generation
    ├→ H0: QR Prediction + Analog Blending
    ├→ Append H0 → Recompute Features
    ├→ H1: QR Prediction (updated features) + Analog Blending
    ├→ Append H1 → Recompute Features
    ├→ H2: QR Prediction (updated features) + Analog Blending
    ├→ Append H2 → Recompute Features
    ├→ H3: QR Prediction (updated features) + Analog Blending
    ├→ Quantile Construction (23 levels)
    └→ Monotonicity Enforcement (isotonic regression)
    ↓
Output File Generation
    ├→ CDC Submission
    ├→ Software Implementation
    └→ Evaluation
    ↓
Visualization & Diagnostics
```

## Model Validation

### Metrics Produced
- **MAE (log scale)** - Tuning metric
- **Residual distributions** - By phase and horizon
- **Analog diagnostics** - Mean distance, weight distribution, rebound/decline shares
- **Blend composition** - QR vs. analog contribution per forecast
- **AI Score** - Alignment Index: min(actual, forecast) / max(actual, forecast)

### Diagnostics Reported
- Latest observed week and value
- Detected phase and phase-specific state variables
- Forecast medians by horizon
- Analog pool statistics (n, mean change, weight concentration)
- QR vs. analog blend weights
- Peak probability forecasts (prob_peak_by_1/2/3/4)


## License

Creative Commons Attribution 4.0 International (CC-BY-4.0)

---

**Last Updated**: June 14, 2026  
**Model Version**: 1.9
