# FluSight: Recursive Quantile Regression Forecasting for Influenza Hospitalizations

## Overview

This repository contains the model implementation for forecasting weekly influenza hospital admissions for South Carolina and for the U.S. national total. The workflow is a recursive quantile-regression model that blends trend-aware quantile regression with historical analogs, with EHR and cross-state signals available as additional predictors in the v2.0 implementation.



## Key Features

### 1. Recursive Forecasting Architecture
- Fit the H0 model as a one-step quantile regression on the current feature set
- Append the predicted value to the series and recompute lag, slope, rolling-statistic, and phase features
- Reuse the same QR model for recursive forecasting at H1, H2, and H3
- Forecast horizons supported by the implementation are 0, 1, 2, and 3 weeks ahead

### 2. Hybrid Forecasting Approach
- Quantile Regression (QR): recursive one-step model using covariates engineered from the observed and synthetic series
- Analog Correction: historical periods with similar epidemiologic state are used to adjust predictions
- Phase-aware blending: decline and peak conditions can modify the QR/analog balance
- External signals: the v2.0 workflow supports a cross-state signal (`pct_states_rising`) and optional pooled-state training

### 3. Feature Engineering
The model uses engineered features including:
- Lag features: 1-, 2-, 3-, 4-, 5-, 8-, and 12-week lags
- Slope features: 1-, 2-, 3-, and 4-week slopes and prior-slope terms
- Acceleration features: 2- and 3-week acceleration terms
- Decline persistence indicators: sustained negative change conditions
- Rolling statistics: 4-, 8-, and 12-week means, standard deviations, maxima, and minima
- Position in cycle: values relative to recent 8- and 12-week ranges
- Peak dynamics: weeks since recent peak, drop from recent peak, and recent peak ratio
- EHR indicators: hospital, positive-test, and testing series slopes with decline confirmation
- Seasonal features: week-of-year sine/cosine terms and holiday indicators

### 4. Phase-Based Classification
The model identifies and adapts to epidemic phases including:
1. FAST_POST_PEAK_DECLINE
2. POST_PEAK_DECLINE
3. ACCELERATING_DECLINE
4. SLOWING_DECLINE
5. ACCELERATING_INCREASE
6. SLOWING_INCREASE
7. NEAR_PEAK_DECELERATION
8. VOLATILE_PLATEAU
9. STABLE_PLATEAU

Phase detection drives adaptive interval scaling and QR/analog balance adjustments.

### 5. Data Integration
The v2.0 implementation integrates multiple data sources:
- CDC weekly influenza hospital admission data for state-level panels and the U.S. total
- South Carolina target data extracted from the CDC panel using the SC FIPS/location mapping
- EHR surveillance data 
- External state-level signal for the fraction of other states rising (`pct_states_rising`)
- Optional pooling of other states for additional training rows, rescaled to South Carolina population

## Repository Contents

### Files

- `QR.r` - Main forecasting pipeline
  - Data loading and merging
  - Feature engineering and origin-specific phase assignment
  - Leakage-protected time-series logic
  - Validation and tuning for the recursive model
  - Recursive forecast generation for horizons 0, 1, 2, and 3
  - Quantile interval construction and forecast outputs
  - Diagnostics and visualization
  - CDC submission and evaluation outputs

- `DMAPRIME-QR.yml` - Model metadata and documentation
  - Team and model metadata
  - Model version and methods statement
  - Data input description
  - Funding and licensing information

- `README.md` - Project overview and model documentation

## Workflow

### 1. Data Preparation
- Load CDC panel data for states and U.S. totals
- Load the South Carolina target series and national-total series as needed
- Load weekly EHR series
- Align datasets by week and standardize column names
- Create the weekly panel used for forecasting and pooled-state training

### 2. Feature Engineering
- Compute target and predictor time-series features from observed and recursive forecast rows
- Apply log-transformed target features and safe handling for zeros and missing values
- Calculate rolling statistics, lags, slopes, acceleration, and seasonal terms
- Assign phase and regime flags using origin-specific thresholds, without exposing future observations to the model

### 3. Validation and Tuning
- Use validation origins within the configured tuning window
- Train the QR model on historical rows with labels known by the origin
- Evaluate candidate analog settings and blend weights
- Select tuning parameters by minimizing MAE in log space for the recursive model configuration

### 4. Model Training
- Fit the QR model using the available historical training data for each forecast horizon
- Reuse the same recursive forecasting recipe for H0 through H3
- Build per-horizon models and residual diagnostics from training rows where labels are known at the forecast origin

### 5. Recursive Forecast Generation
- H0 (nowcast): predict from the current feature row
- Append the H0 prediction to create a synthetic next row and recompute features
- Forecast H1 from the updated feature row
- Repeat recursively for H2 and H3
- Blend QR predictions with analog-based adjustments under phase-aware settings
- Construct a set of quantile forecasts for the final output

### 6. Output Generation
- Forecast targets include weekly influenza hospital admissions for South Carolina and the U.S. national total
- Generated outputs include quantile forecasts at horizons 0, 1, 2, and 3
- Export files are produced for CDC submission, software implementation, and evaluation workflows

### 7. Visualization and Diagnostics
- Plot recent observed series and forecast medians across the forecast window
- Display interval bands and phase labels
- Summarize analog diagnostics, peak probabilities, and model components

## Key Hyperparameters

| Parameter | Value | Purpose |
|-----------|-------|---------|
| `quantiles` | 0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99 | Forecast quantile levels |
| `horizons` | 0, 1, 2, 3 | Forecast horizon values supported by the implementation |
| `candidate_k` | 8, 12, 16, 20, 25, 30 | Candidate analog counts |
| `candidate_bandwidth` | 0.50, 0.75, 1.00, 1.25, 1.50, 2.00 | Candidate analog decay widths |
| `candidate_blend_qr` | 0.0 to 1.0 by 0.10 | Candidate QR/analog blending levels |
| `train_cut` | 2022-09-01 | Training cutoff for validation work |
| `valid_end` | 2025-10-01 | Validation-end date used in the code |
| `eval_start` | 2025-11-06 | Evaluation-start date used in the code |
| `eval_end` | 2026-03-07 | Evaluation-end date used in the code |

## Requirements

### R Packages
```r
suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(quantreg)
  library(tidyr)
  library(slider)
  library(scales)
})
```

### Data Files
The implementation uses the following data inputs as represented in `QR.r`:
- CDC admissions panel data for states and U.S. national total
- Prisma Health weekly influenza EHR and testing data
- Output directories for CDC submission, software implementation, and evaluation files

**Note**: Data paths are hardcoded in the implementation and are not repository-managed file paths.

## Workflow Diagram

```text
Data Loading
    ↓
CDC + EHR + external signal preparation
    ↓
Feature Engineering
    ↓
Origin-specific phase and regime assignment
    ↓
Validation and tuning
    ↓
Recursive QR forecast generation
    ├─ H0: current-week prediction
    ├─ append prediction and recompute features
    ├─ H1: 1-week-ahead prediction
    ├─ append prediction and recompute features
    ├─ H2: 2-week-ahead prediction
    ├─ append prediction and recompute features
    ├─ H3: 3-week-ahead prediction
    └─ blend QR + analog outputs
    ↓
Quantile forecast construction
    ↓
Output generation
    ├─ CDC submission format
    ├─ software implementation file
    └─ evaluation file
    ↓
Visualization and diagnostics
```

## Model Validation

### Metrics Produced
- MAE in log space for tuning
- Residual diagnostics by phase and horizon
- Analog diagnostics including distance, weighting, and phase composition
- Blend composition summaries for QR vs. analog contributions
- Forecast coverage and retrospective evaluation summaries

### Diagnostics Reported
- Latest observed week and value
- Detected phase and support variables
- Forecast medians by horizon
- Summary of analog pool composition
- Peak probability outputs for the forecast path

## License

Creative Commons Attribution 4.0 International (CC-BY-4.0)

---

**Last Updated**: 2026-10-10  
**Model Version**: 2.0
