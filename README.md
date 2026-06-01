# FluSight: Hybrid Quantile Regression Forecasting for Influenza Hospitalizations

## Overview

This repository contains a **hybrid forecasting model** for predicting weekly influenza-related inpatient hospitalizations in South Carolina. The model combines **direct horizon-specific quantile regression** with **decline-aware analog correction** to produce probabilistic forecasts with phase-conditioned prediction intervals.

## Model Name

- **Full Name**: DMAPRIME Quantile Regression (DMAPRIME-QR)
- **Team**: CPHMR (Center for Public Health and Medical Research), Clemson University
- **Model Abbr**: QR
- **Version**: 1.8
- **License**: CC-BY-4.0

## Key Features

### 1. Hybrid Forecasting Approach
- **Quantile Regression (QR)**: Horizon-specific direct quantile regression models for each forecast week (0-3 weeks ahead)
- **Analog Correction**: Identifies historical periods with similar epidemiological states and uses their outcomes to adjust predictions
- **Decline-Aware Weighting**: Applies phase-specific penalties during decline periods to reduce overly optimistic rebound predictions

### 2. Rich Feature Engineering
The model uses 38+ engineered features including:
- **Lag features**: Target values at 1, 2, 3, 4, 5, 8, and 12-week lags
- **Slope features**: 1-4 week slopes and slope changes to capture trajectory
- **Acceleration features**: Rate of change in slopes (2-week and 3-week acceleration)
- **Decline persistence**: Indicators for sustained negative trends
- **Rolling statistics**: 4-, 8-, and 12-week means and standard deviations
- **Position in cycle**: Where current levels sit in the 12-week rolling range
- **Peak dynamics**: Weeks since 12-week peak, distance from peak, and post-peak decline speed
- **EHR indicators**: Hospitalization and positive test slopes with decline confirmation
- **Seasonal features**: Week of year harmonics (sin/cos) and holiday indicators

### 3. Phase-Based Classification
The model identifies and adapts to 9 epidemic phases:
1. **FAST_POST_PEAK_DECLINE** - Rapid post-peak decline
2. **POST_PEAK_DECLINE** - Standard post-peak decline
3. **ACCELERATING_DECLINE** - Declining trend with acceleration
4. **SLOWING_DECLINE** - Declining trend with deceleration
5. **ACCELERATING_INCREASE** - Rising trend with acceleration
6. **SLOWING_INCREASE** - Rising trend with deceleration
7. **VOLATILE_PLATEAU** - Stable levels with high variability
8. **STABLE_PLATEAU** - Stable levels with low variability
9. **UNCERTAIN** - Ambiguous dynamics

Phase detection drives adaptive interval scaling and QR floor constraints.

### 4. Data Integration
The model integrates multiple data sources:
- **CDC NHSN Data**: Official weekly flu hospitalizations (location 45 = South Carolina)
- **Electronic Health Record (EHR) Data**: 
  - MUSC (Medical University of South Carolina) weekly influenza data
  - Prisma Health weekly influenza surveillance data
- **Viral Testing Data**: Rapid assessment of circulation intensity via positive test rates

## Repository Contents

### Files

- **`QR.r`** - Main forecasting pipeline (993 lines)
  - Data loading and merging
  - Feature engineering
  - Phase labeling
  - Model training and tuning via rolling validation
  - Forecast generation with quantile intervals
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
- Load CDC, MUSC, and Prisma weekly hospitalization data
- Align datasets by week
- Standardize column names

### 2. **Feature Engineering**
- Compute 38+ predictive features from historical time series
- Apply log transformation with safe handling of zero/missing values
- Calculate rolling statistics and seasonal components
- Determine epidemic phase

### 3. **Rolling Validation** (Lines 440-587)
- For each week in the validation period (Sept 2022 – Oct 2025):
  - Train quantile regression models on all prior data
  - Identify analogs from historical pool
  - Evaluate multiple tuning parameters:
    - `k`: Number of analogs (candidate values: 8, 12, 16, 20, 25, 30)
    - `bandwidth`: Analog weight decay (0.50–2.00)
    - `blend_qr`: Mixing proportion (0–1.0 in 0.1 increments)
  - Compute Mean Absolute Error (MAE) in log space
  - Select best hyperparameters per horizon

### 4. **Model Training** (Lines 594-618)
- Fit final quantile regression models for horizons 0–3 using all training data
- Capture residuals by phase for uncertainty quantification

### 5. **Forecast Generation** (Lines 622-772)
- Generate predictions at 23 quantile levels (0.01–0.99)
- Blend QR and analog predictions using phase-specific weights
- Apply phase-conditioned interval scaling
- Use residual distributions to compute prediction intervals
- Ensure monotonicity via isotonic regression

### 6. **Output Generation**
- **CDC Submission File** (Section 9):
  - Format: CDC FluSight Hub submission format
  - Contains all quantile predictions (0–4 weeks ahead)
  
- **Software Implementation File** (Section 9):
  - Combines historical data (RFA + EHR sources) with forecasts
  - Ready for deployment in forecasting infrastructure
  
- **Evaluation File** (Section 10):
  - Training data and median-only test predictions
  - Used for model performance assessment

### 7. **Visualization**
- Plots observed time series (12 weeks back)
- Overlays median forecast and 50% / 95% prediction intervals
- Marks forecast origin line
- Indicates detected phase in title

## Key Hyperparameters

| Parameter | Default Range | Purpose |
|-----------|---------------|---------|
| `candidate_k` | 8, 12, 16, 20, 25, 30 | Number of analog neighbors |
| `candidate_bandwidth` | 0.50–2.00 | Decay rate for analog weights |
| `candidate_blend_qr` | 0.0–1.0 (step 0.1) | QR-analog mixing proportion |
| `quantiles` | 0.01–0.99 (23 levels) | Forecast quantile levels |
| `horizons` | 0, 1, 2, 3 | Weeks ahead to forecast |

## Decline-Aware Adjustments

During decline phases:
- **QR Floor**: Minimum blend weight for QR (e.g., 65% for FAST_POST_PEAK_DECLINE)
- **Rebound Penalty**: Reduce weights of analogs predicting increases by 80%
- **Weak Decline Penalty**: Reduce weights of analogs with marginal declines
- **Interval Scaling**: Narrow intervals during decline to reduce uncertainty

## Tuning Metrics

- **Metric**: Mean Absolute Error (MAE) in log-transformed space
- **Validation Period**: Sept 1, 2022 – Oct 1, 2025
- **Candidates Evaluated**: 6 × 6 × 11 = **396 configurations** per horizon
- **Selection**: Best MAE per horizon

## Output Format

### CDC Submission File
```
reference_date, target, horizon, location, target_end_date, output_type, output_type_id, value
2025-06-07, wk inc flu hosp, 0, 45, 2025-06-07, quantile, 0.025, 125
2025-06-07, wk inc flu hosp, 0, 45, 2025-06-07, quantile, 0.500, 150
...
```

### Software Implementation / Evaluation Files
Contains columns:
- `reference_date` - Week of submission
- `target` - Target specification (e.g., "wk inc flu hosp")
- `target_end_date` - Week ending date
- `location_general` - "state"
- `location` - "SC"
- `value` - Point estimate or historical value
- `disease`, `population`, `outcome_measure` - Metadata
- `output_type`, `output_type_id` - Quantile specification

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
```

### Data Files
- RFA weekly influenza data (2017–2025)
- CDC hospital admissions (NHSN)
- MUSC weekly influenza surveillance
- Prisma Health weekly influenza surveillance

**Note**: Data paths are currently hardcoded to local Box Cloud Storage directories. Modify paths before running.

## Workflow Diagram

```
Data Loading
    ↓
Feature Engineering (38+ features)
    ↓
Phase Classification
    ↓
Rolling Validation (396 configs × 4 horizons)
    ↓
Hyperparameter Selection (best MAE)
    ↓
Final Model Training
    ↓
Current Week Forecast Generation
    ├→ QR Predictions
    ├→ Analog Search & Weighting
    ├→ Blended Point Forecast
    ├→ Interval Construction
    └→ Quantile Levels (23 points)
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

### Diagnostics Reported
- Latest observed week and value
- Detected phase
- Phase-specific state variables (slope, acceleration, decline persistence)
- Forecast medians
- Analog pool statistics (n, mean change, weight concentration)

## Future Extensions

Potential improvements:
- Regional/substate level forecasts
- Multi-target modeling (tests, outpatient visits)
- Ensemble with other methods
- Real-time data integration
- Dynamic phase transitions
- Batch forecast generation

## Contact

**Principal Investigator**: Tanvir Ahammed  
**Affiliation**: Clemson University, South Carolina  
**Email**: tahamme@g.clemson.edu

## Funding

Supported by the Center for Forecasting and Outbreak Analytics of the Centers for Disease Control and Prevention (CDC) under award number NU38FT000011 and the National Library of Medicine of the National Institutes of Health.

## License

Creative Commons Attribution 4.0 International (CC-BY-4.0)

---

**Last Updated**: June 1, 2026  
**Model Version**: 1.8
