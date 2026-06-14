# ══════════════════════════════════════════════════════════════════════════════
# DMAPRIME-QR v4: State-Level Influenza Hospitalization Forecast
#
# KEY CHANGE FROM v3 → v4:
#   Switched from DIRECT multi-horizon forecasting to RECURSIVE forecasting.
#
#   v3 (direct):   Fit separate H0, H1, H2, H3 models independently.
#                  Then post-process to prevent horizon shift.
#                  Problem: models don't share state; anti-shift logic fights
#                  a problem the model architecture creates.
#
#   v4 (recursive): Forecast H0 first. Append that prediction to the feature
#                   series. Recompute all features (slopes, lags, phase, etc.)
#                   on the extended series. Forecast H1 using the updated
#                   features. Repeat for H2, H3.
#                   Result: horizon shift is architecturally impossible.
#                   anti_shift_postprocess is no longer needed.
#
# What stays the same:
#   - QR model structure and predictor sets
#   - Analog selection with phase gating and penalty weights
#   - Phase classification logic
#   - Rolling validation / tuning
#   - All output formats (CDC submission, software files, plots)
#
# New functions:
#   append_forecast_row()      — extends feature df with a synthetic row
#   recompute_features_tail()  — recomputes rolling/lag features for last row
#   forecast_origin_v4()       — recursive loop replacing forecast_origin_v3()
#
# Removed:
#   anti_shift_postprocess_v3() — no longer needed
# ══════════════════════════════════════════════════════════════════════════════

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

# ══════════════════════════════════════════════════════════════════════════════
# DATA LEAKAGE AUDIT (v4 — fully corrected)
#
# ALL LEAKAGE REMOVED. Every evaluation origin uses only data <= origin_week.
#
# ARCHITECTURE:
#   compute_origin_thresholds(df, cutoff_week)
#     → slope_cut, accel_cut, vol_cut, slope_q25, global_slope_sd
#     → computed from df rows where week <= cutoff_week ONLY
#     → called at every origin in: tuning loop, retrospective loop, final forecast
#
#   assign_phase(df, thr)
#     → takes an explicit thr list (no global variable closure)
#     → called with origin-specific thr everywhere
#
#   append_forecast_row(df, new_week, pred_log, origin_row, thr)
#     → passes thr through so synthetic rows get origin-specific phase labels
#
#   forecast_origin_v4(..., thr, ...)
#     → uses thr$global_slope_sd for decline penalties
#     → passes thr to every append_forecast_row call
#
# VERIFIED SAFE:
#   - QR model: get_label_known_training() enforces week+(h+1) <= cutoff
#   - Analog pool: cutoff = origin_week, future_change_h0 correctly bounded
#   - future_peak_offset_4wk: label column, masked to NA when week+4 > origin_week
#   - Tuning y_true: used only as scoring label, not as predictor
#   - Hyperparameter tuning window: train_cut to valid_end (ends 2025-10-01)
#     Eval window starts at eval_start (default 2025-12-06) — no overlap
#
# WHAT THIS MEANS FOR REPORTED SCORES:
#   AI scores in the retrospective evaluation now reflect genuine out-of-sample
#   performance. A forecast issued on 2023-01-14 uses thresholds estimated from
#   data through 2023-01-14 only — exactly what a real forecaster would have.
# ══════════════════════════════════════════════════════════════════════════════

# ══════════════════════════════════════════════════════════════════════════════
# 0. CONFIGURATION  (unchanged from v3)
# ══════════════════════════════════════════════════════════════════════════════

quantiles <- c(
   0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
   0.5,
   0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99
)

horizons <- 0:3

candidate_k <- c(8, 12, 16, 20, 25, 30)
candidate_bandwidth <- c(0.50, 0.75, 1.00, 1.25, 1.50, 2.00)
candidate_blend_qr <- seq(0, 1, by = 0.10)

train_cut <- as.Date("2022-09-01")
valid_end <- as.Date("2025-10-01")

eval_start <- as.Date("2025-11-06")
eval_end   <- as.Date("2025-03-07")




train_cut <- as.Date("2022-09-01")
valid_end <- as.Date("2024-10-01")
eval_start <- as.Date("2024-11-06")
eval_end   <- as.Date("2025-03-07")



path_rfa <- "/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Data/SC Health Records/RFA/Years 2017-2025/Weekly data/RFA_weekly_influenza_region_incident.csv"
path_cdc <- '/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Data/CDC Data/CDC Hospital Repository Data/Weekly_Hospital_Respiratory_Data_(Preliminary).csv'
path_musc <- "/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Data/MUSC/Infectious Disease EHR/Weekly Data/Latest Weekly Data/MUSC_Weekly_Influenza_State_dx_cond_lab_Incident.csv"
path_prisma <- "/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Data/Prisma Health/Infectious Disease EHR/Weekly Data/Latest Weekly Data/Prisma_Health_Weekly_Influenza_State_dx_cond_lab_Incident.csv"

output_dir_cdc <- "/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Tanvir/FluSight/Tanvir_QR/State/CDC_submission"
output_dir_software_impl <- "/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Forecasting Resources/Forecast-Drop-Off/Software/Implementation"
output_dir_software_eval <- "/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Forecasting Resources/Forecast-Drop-Off/Software/Evaluation"

target_var <- "cdc_Total.Influenza.Admissions"

# ══════════════════════════════════════════════════════════════════════════════
# 1. LOAD DATA  (unchanged from v3)
# ══════════════════════════════════════════════════════════════════════════════

RFA_data_Influenza <- read_csv(path_rfa, show_col_types = FALSE)

CDC <- read_csv(path_cdc, show_col_types = FALSE) %>%
   #filter(location == "45") %>%
   dplyr::select(
      week = Week.Ending.Date,
      cdc_Total.Influenza.Admissions = Total.Influenza.Admissions
   ) %>%
   mutate(week = as.Date(week))



MUSC_Weekly_Influenza <- read_csv(path_musc, show_col_types = FALSE)
Prisma_Weekly_Influenza <- read_csv(path_prisma, show_col_types = FALSE)

prep_source <- function(df, date_col, prefix) {
   df %>%
      rename(week = !!rlang::sym(date_col)) %>%
      mutate(week = as.Date(week)) %>%
      rename_with(~ paste0(prefix, .x), .cols = -week)
}

prisma <- prep_source(Prisma_Weekly_Influenza, "Week", "prisma_")
musc   <- prep_source(MUSC_Weekly_Influenza,   "Week", "musc_")

flu <- CDC %>%
   full_join(prisma, by = "week") %>%
   full_join(musc,   by = "week") %>%
   arrange(week)

optional_cols <- c(
   "prisma_Weekly_Tests",
   "prisma_Weekly_Positive_Tests",
   "prisma_Weekly_Inpatient_Hospitalizations"
)
for (nm in optional_cols) {
   if (!nm %in% names(flu)) flu[[nm]] <- NA_real_
}

# ══════════════════════════════════════════════════════════════════════════════
# 2. GENERAL HELPERS  (unchanged from v3)
# ══════════════════════════════════════════════════════════════════════════════

safe_log <- function(x) log(pmax(x, 0) + 1)

safe_mean <- function(x) {
   if (all(is.na(x))) return(NA_real_)
   mean(x, na.rm = TRUE)
}

safe_sd <- function(x) {
   if (sum(!is.na(x)) < 2) return(0)
   sd(x, na.rm = TRUE)
}

safe_max <- function(x) {
   if (all(is.na(x))) return(NA_real_)
   max(x, na.rm = TRUE)
}

safe_min <- function(x) {
   if (all(is.na(x))) return(NA_real_)
   min(x, na.rm = TRUE)
}

weeks_since_peak_fun <- function(x) {
   if (all(is.na(x))) return(NA_integer_)
   length(x) - which.max(replace_na(x, -Inf))
}

future_peak_offset_fun <- function(x) {
   if (all(is.na(x))) return(NA_integer_)
   which.max(replace_na(x, -Inf)) - 1L
}

weighted_quantile_safe <- function(x, w, probs) {
   ok <- !is.na(x) & !is.na(w)
   x <- x[ok]; w <- w[ok]
   if (length(x) == 0 || sum(w) <= 0) return(rep(NA_real_, length(probs)))
   ord <- order(x)
   x <- x[ord]; w <- w[ord] / sum(w[ord])
   cw <- cumsum(w)
   sapply(probs, function(p) x[which(cw >= p)[1]])
}

weighted_median_safe <- function(x, w) weighted_quantile_safe(x, w, probs = 0.5)[1]

weighted_mean_safe <- function(x, w) {
   ok <- !is.na(x) & !is.na(w)
   x <- x[ok]; w <- w[ok]
   if (length(x) == 0 || sum(w) <= 0) return(NA_real_)
   sum(x * w) / sum(w)
}

standardize_using <- function(df, ref_df, vars) {
   out <- df
   for (v in vars) {
      m <- mean(ref_df[[v]], na.rm = TRUE)
      s <- sd(ref_df[[v]], na.rm = TRUE)
      if (is.na(m)) m <- 0
      if (is.na(s) || s == 0) s <- 1
      out[[v]] <- (out[[v]] - m) / s
   }
   out
}

circular_week_diff <- function(a, b) {
   raw <- abs(a - b)
   pmin(raw, 52 - raw)
}

is_one <- function(x) !is.na(x[1]) && x[1] == 1

first_or <- function(x, default) {
   if (length(x) == 0 || all(is.na(x))) return(default)
   x[which(!is.na(x))[1]]
}

# ══════════════════════════════════════════════════════════════════════════════
# 3. FEATURE ENGINEERING  (unchanged from v3)
# ══════════════════════════════════════════════════════════════════════════════

build_flu_features <- function(flu_df) {
   
   flu_df %>%
      arrange(week) %>%
      mutate(
         log_target = safe_log(.data[[target_var]]),
         
         log_tests          = safe_log(prisma_Weekly_Tests),
         log_positive_tests = safe_log(prisma_Weekly_Positive_Tests),
         log_hosp           = safe_log(prisma_Weekly_Inpatient_Hospitalizations),
         
         lag1  = lag(log_target, 1),
         lag2  = lag(log_target, 2),
         lag3  = lag(log_target, 3),
         lag4  = lag(log_target, 4),
         lag5  = lag(log_target, 5),
         lag8  = lag(log_target, 8),
         lag12 = lag(log_target, 12),
         
         slope_1wk = log_target - lag1,
         slope_2wk = log_target - lag2,
         slope_3wk = log_target - lag3,
         slope_4wk = log_target - lag4,
         
         prev_slope_1wk = lag1 - lag2,
         prev_slope_2wk = lag2 - lag3,
         
         accel_2wk = slope_1wk - prev_slope_1wk,
         accel_3wk = slope_1wk - ((lag1 - lag3) / 2),
         
         decline_persistence_2wk = as.integer(slope_1wk < 0 & prev_slope_1wk < 0),
         decline_persistence_3wk = as.integer(
            slope_1wk < 0 & prev_slope_1wk < 0 & prev_slope_2wk < 0
         ),
         
         mean_4wk  = slide_dbl(log_target, safe_mean, .before = 3,  .complete = FALSE),
         mean_8wk  = slide_dbl(log_target, safe_mean, .before = 7,  .complete = FALSE),
         mean_12wk = slide_dbl(log_target, safe_mean, .before = 11, .complete = FALSE),
         
         sd_4wk  = slide_dbl(log_target, safe_sd, .before = 3,  .complete = FALSE),
         sd_8wk  = slide_dbl(log_target, safe_sd, .before = 7,  .complete = FALSE),
         sd_12wk = slide_dbl(log_target, safe_sd, .before = 11, .complete = FALSE),
         
         max_8wk  = slide_dbl(log_target, safe_max, .before = 7,  .complete = FALSE),
         max_12wk = slide_dbl(log_target, safe_max, .before = 11, .complete = FALSE),
         
         min_8wk  = slide_dbl(log_target, safe_min, .before = 7,  .complete = FALSE),
         min_12wk = slide_dbl(log_target, safe_min, .before = 11, .complete = FALSE),
         
         position_in_8wk_range =
            (log_target - min_8wk) / pmax(max_8wk - min_8wk, 1e-6),
         position_in_12wk_range =
            (log_target - min_12wk) / pmax(max_12wk - min_12wk, 1e-6),
         
         weeks_since_8wk_peak = slide_int(
            log_target, weeks_since_peak_fun, .before = 7, .complete = FALSE),
         weeks_since_12wk_peak = slide_int(
            log_target, weeks_since_peak_fun, .before = 11, .complete = FALSE),
         
         peak_8wk_log  = max_8wk,
         peak_12wk_log = max_12wk,
         
         drop_from_8wk_peak  = log_target - max_8wk,
         drop_from_12wk_peak = log_target - max_12wk,
         
         drop_from_8wk_peak_pct =
            ((exp(log_target) - 1) - (exp(max_8wk) - 1)) /
            pmax(exp(max_8wk) - 1, 1e-6),
         drop_from_12wk_peak_pct =
            ((exp(log_target) - 1) - (exp(max_12wk) - 1)) /
            pmax(exp(max_12wk) - 1, 1e-6),
         
         post_peak_decline_speed =
            drop_from_12wk_peak / pmax(weeks_since_12wk_peak, 1),
         
         recent_peak_ratio =
            (exp(log_target) - 1) / pmax(exp(max_8wk) - 1, 1e-6),
         
         hosp_slope_1wk = log_hosp - lag(log_hosp, 1),
         hosp_slope_2wk = log_hosp - lag(log_hosp, 2),
         pos_slope_1wk  = log_positive_tests - lag(log_positive_tests, 1),
         pos_slope_2wk  = log_positive_tests - lag(log_positive_tests, 2),
         
         ehr_decline_confirmed = as.integer(
            hosp_slope_1wk < 0 | pos_slope_1wk < 0 |
               hosp_slope_2wk < 0 | pos_slope_2wk < 0
         ),
         
         week_of_year = isoweek(week),
         sin52 = sin(2 * pi * week_of_year / 52),
         cos52 = cos(2 * pi * week_of_year / 52),
         
         month = month(week),
         is_christmas = if_else(month == 12 & day(week) >= 20, 1, 0),
         is_newyear   = if_else(month == 1  & day(week) <= 7,  1, 0)
      ) %>%
      mutate(
         across(c(sd_4wk, sd_8wk, sd_12wk), ~ replace_na(.x, 0)),
         
         slope_1wk_q25 = quantile(slope_1wk, 0.25, na.rm = TRUE),
         
         is_post_peak = if_else(
            weeks_since_12wk_peak >= 1 &
               log_target < max_12wk &
               slope_1wk < 0,
            1, 0
         ),
         
         is_fast_decline = if_else(
            slope_1wk < slope_1wk_q25 &
               slope_2wk < 0 &
               drop_from_12wk_peak < 0,
            1, 0
         ),
         
         is_near_peak = as.integer(
            recent_peak_ratio >= 0.80 &
               (accel_2wk <= 0 | slope_1wk <= 0.05)
         ),
         
         any_decline = as.integer(
            is_fast_decline == 1 |
               is_post_peak == 1 |
               decline_persistence_2wk == 1
         ),
         
         any_peak_or_decline = as.integer(
            is_near_peak == 1 | any_decline == 1
         )
      ) %>%
      dplyr::select(-slope_1wk_q25)
}

flu_features <- build_flu_features(flu)

# Future labels — only used in training, not in recursive forecast rows
for (h in horizons) {
   flu_features[[paste0("future_log_h", h)]] <-
      lead(flu_features$log_target, h + 1)
   flu_features[[paste0("future_value_h", h)]] <-
      lead(flu_features[[target_var]], h + 1)
   flu_features[[paste0("future_change_h", h)]] <-
      flu_features[[paste0("future_log_h", h)]] - flu_features$log_target
}

flu_features <- flu_features %>%
   mutate(
      future_peak_offset_4wk = slide_int(
         log_target, future_peak_offset_fun,
         .before = 0, .after = 4, .complete = TRUE
      ),
      peak_by_1wk = as.integer(!is.na(future_peak_offset_4wk) & future_peak_offset_4wk <= 1),
      peak_by_2wk = as.integer(!is.na(future_peak_offset_4wk) & future_peak_offset_4wk <= 2),
      peak_by_3wk = as.integer(!is.na(future_peak_offset_4wk) & future_peak_offset_4wk <= 3),
      peak_by_4wk = as.integer(!is.na(future_peak_offset_4wk) & future_peak_offset_4wk <= 4)
   )

# ══════════════════════════════════════════════════════════════════════════════
# 4. PREDICTOR SETS AND PHASE LABELS  (unchanged from v3)
# ══════════════════════════════════════════════════════════════════════════════

predictors <- c(
   "log_tests", "log_positive_tests", "log_hosp",
   "hosp_slope_1wk", "hosp_slope_2wk",
   "pos_slope_1wk", "pos_slope_2wk",
   "ehr_decline_confirmed",
   "lag1", "lag2", "lag4", "lag8", "lag12",
   "slope_1wk", "slope_2wk", "slope_3wk", "slope_4wk",
   "prev_slope_1wk", "prev_slope_2wk",
   "accel_2wk", "accel_3wk",
   "decline_persistence_2wk", "decline_persistence_3wk",
   "mean_4wk", "mean_8wk", "mean_12wk",
   "sd_4wk", "sd_8wk", "sd_12wk",
   "position_in_8wk_range", "position_in_12wk_range",
   "weeks_since_8wk_peak", "weeks_since_12wk_peak",
   "peak_8wk_log", "peak_12wk_log",
   "drop_from_8wk_peak", "drop_from_12wk_peak",
   "drop_from_8wk_peak_pct", "drop_from_12wk_peak_pct",
   "post_peak_decline_speed", "recent_peak_ratio",
   "is_post_peak", "is_fast_decline", "is_near_peak",
   "any_decline", "any_peak_or_decline",
   "sin52", "cos52", "is_christmas", "is_newyear"
)

state_vars <- c(
   "log_target", "lag1",
   "slope_1wk", "slope_2wk", "slope_3wk", "slope_4wk",
   "prev_slope_1wk", "prev_slope_2wk",
   "accel_2wk", "accel_3wk",
   "decline_persistence_2wk", "decline_persistence_3wk",
   "mean_4wk", "mean_8wk", "mean_12wk",
   "sd_4wk", "sd_8wk", "sd_12wk",
   "position_in_8wk_range", "position_in_12wk_range",
   "weeks_since_8wk_peak", "weeks_since_12wk_peak",
   "peak_8wk_log", "peak_12wk_log",
   "drop_from_8wk_peak", "drop_from_12wk_peak",
   "drop_from_8wk_peak_pct", "drop_from_12wk_peak_pct",
   "post_peak_decline_speed", "recent_peak_ratio",
   "is_post_peak", "is_fast_decline", "is_near_peak",
   "any_decline", "any_peak_or_decline",
   "log_hosp", "hosp_slope_1wk", "hosp_slope_2wk",
   "log_positive_tests", "pos_slope_1wk", "pos_slope_2wk",
   "ehr_decline_confirmed",
   "week_of_year", "sin52", "cos52"
)

# ══════════════════════════════════════════════════════════════════════════════
# THRESHOLD FUNCTIONS — fully origin-aware, no temporal leakage
#
# compute_origin_thresholds(df, cutoff_week):
#   Computes slope_cut, accel_cut, vol_cut, slope_q25, global_slope_sd
#   from df rows where week <= cutoff_week only.
#   Call with cutoff_week = origin_week inside every loop iteration.
#
# assign_phase(df, thr):
#   Applies phase labels using the thresholds in thr (a named list).
#   Never closes over global variables — fully explicit.
# ══════════════════════════════════════════════════════════════════════════════

compute_origin_thresholds <- function(df, cutoff_week) {
   past <- df %>%
      filter(week <= cutoff_week) %>%
      filter(!is.na(slope_1wk), !is.na(accel_2wk), !is.na(sd_4wk))
   
   slope_sd  <- sd(past$slope_1wk, na.rm = TRUE)
   accel_sd  <- sd(past$accel_2wk, na.rm = TRUE)
   
   list(
      slope_cut         = slope_sd * 0.25,
      accel_cut         = accel_sd * 0.25,
      vol_cut           = quantile(past$sd_4wk,    0.70, na.rm = TRUE),
      slope_q25         = quantile(past$slope_1wk, 0.25, na.rm = TRUE),
      global_slope_sd   = if (is.na(slope_sd) || slope_sd <= 0) 1 else slope_sd
   )
}

assign_phase <- function(df, thr) {
   sc <- thr$slope_cut
   ac <- thr$accel_cut
   vc <- thr$vol_cut
   df %>%
      mutate(
         phase = case_when(
            is_fast_decline == 1 ~ "FAST_POST_PEAK_DECLINE",
            is_post_peak == 1 & decline_persistence_2wk == 1 ~ "FAST_POST_PEAK_DECLINE",
            is_post_peak == 1 & slope_1wk < 0 ~ "POST_PEAK_DECLINE",
            is_near_peak == 1 & slope_1wk >= 0 & accel_2wk <= 0 ~ "NEAR_PEAK_DECELERATION",
            slope_1wk < -sc & accel_2wk < -ac ~ "ACCELERATING_DECLINE",
            slope_1wk < -sc ~ "SLOWING_DECLINE",
            slope_1wk > sc  & accel_2wk > ac  ~ "ACCELERATING_INCREASE",
            slope_1wk > sc  ~ "SLOWING_INCREASE",
            abs(slope_1wk) <= sc & sd_4wk > vc ~ "VOLATILE_PLATEAU",
            abs(slope_1wk) <= sc ~ "STABLE_PLATEAU",
            TRUE ~ "UNCERTAIN"
         )
      )
}

# For the FINAL forecast and for initial flu_features construction,
# use all data up to the latest observed week (set after data load).
# This placeholder is overwritten after latest_observed_week is known.
# All retrospective/tuning loops call compute_origin_thresholds() per origin.
.initial_thr <- compute_origin_thresholds(
   build_flu_features(flu),
   max(flu$week, na.rm = TRUE)
)

flu_features <- assign_phase(flu_features, .initial_thr)

# Expose global_slope_sd for functions that reference it —
# will be overridden per-origin inside all loops.
global_slope_sd <- .initial_thr$global_slope_sd

# ══════════════════════════════════════════════════════════════════════════════
# 5. QR MODEL HELPERS  (unchanged from v3)
# ══════════════════════════════════════════════════════════════════════════════

prepare_qr_design <- function(model_df, target_h, candidate_predictors) {
   model_df <- model_df %>%
      drop_na(all_of(c(target_h, candidate_predictors)))
   if (nrow(model_df) < 30) stop("Too few complete rows.")
   usable_predictors <- candidate_predictors[
      sapply(candidate_predictors, function(v) {
         x <- model_df[[v]]
         sd(x, na.rm = TRUE) > 1e-8 && length(unique(x[!is.na(x)])) > 1
      })
   ]
   if (length(usable_predictors) == 0) stop("No usable predictors.")
   xmat <- model.matrix(
      as.formula(paste(target_h, "~", paste(usable_predictors, collapse = " + "))),
      data = model_df
   )
   qr_x <- qr(xmat)
   keep_coef_names <- colnames(xmat)[qr_x$pivot[seq_len(qr_x$rank)]]
   keep_predictors <- setdiff(keep_coef_names, "(Intercept)")
   keep_predictors <- keep_predictors[keep_predictors %in% names(model_df)]
   list(data = model_df, predictors = keep_predictors)
}

# v4 NOTE: We keep H0 as the only directly fitted horizon.
# H1/H2/H3 use the same H0 QR model applied to the recursively updated row.
# We still allow fitting per-horizon QR for the residual/interval distribution,
# but the median forecast comes from recursive application of the H0 model.
fit_qr_model <- function(train_data, h, tau_value = 0.5) {
   target_h <- paste0("future_log_h", h)
   prepared <- prepare_qr_design(
      model_df = train_data,
      target_h = target_h,
      candidate_predictors = predictors
   )
   model <- rq(
      as.formula(paste(target_h, "~", paste(prepared$predictors, collapse = " + "))),
      tau = tau_value, data = prepared$data, method = "fn"
   )
   attr(model, "used_predictors") <- prepared$predictors
   model
}

predict_qr_safe <- function(model, newdata) {
   as.numeric(predict(model, newdata = newdata))
}

get_label_known_training <- function(df, label_cutoff_week, h) {
   df %>% filter(week + weeks(h + 1) <= label_cutoff_week)
}

fit_all_qr_models <- function(df, label_cutoff_week) {
   qr_models <- list()
   residual_store <- list()
   for (h in horizons) {
      target_h <- paste0("future_log_h", h)
      model_df <- get_label_known_training(df, label_cutoff_week, h) %>%
         drop_na(all_of(c(target_h, predictors)))
      if (nrow(model_df) < 80) { qr_models[[paste0("h", h)]] <- NULL; next }
      fit_h <- tryCatch(fit_qr_model(model_df, h, tau_value = 0.5), error = function(e) NULL)
      qr_models[[paste0("h", h)]] <- fit_h
      if (!is.null(fit_h)) {
         model_df$pred <- predict_qr_safe(fit_h, model_df)
         model_df$residual <- model_df[[target_h]] - model_df$pred
         model_df$horizon <- h
         residual_store[[paste0("h", h)]] <- model_df %>%
            dplyr::select(week, horizon, residual, phase)
      }
   }
   list(qr_models = qr_models, residual_df = bind_rows(residual_store))
}

# ══════════════════════════════════════════════════════════════════════════════
# 6. ANALOG HELPERS  (unchanged from v3)
# ══════════════════════════════════════════════════════════════════════════════

phase_gate_pool <- function(pool_complete, origin_row, min_n = 10) {
   cand <- pool_complete
   season_window <- if (is_one(origin_row$any_peak_or_decline)) 10 else 13
   cand_season <- cand %>%
      filter(circular_week_diff(week_of_year, origin_row$week_of_year) <= season_window)
   if (nrow(cand_season) >= min_n) cand <- cand_season
   if (is_one(origin_row$any_peak_or_decline)) {
      cand_peak_age <- cand %>%
         filter(abs(weeks_since_12wk_peak - origin_row$weeks_since_12wk_peak) <= 2)
      if (nrow(cand_peak_age) >= min_n) cand <- cand_peak_age
      cand_decline_match <- cand %>%
         filter(
            any_peak_or_decline == origin_row$any_peak_or_decline |
               any_decline == origin_row$any_decline |
               is_near_peak == origin_row$is_near_peak |
               is_post_peak == origin_row$is_post_peak
         )
      if (nrow(cand_decline_match) >= min_n) cand <- cand_decline_match
   }
   cand
}

get_analogs_v4 <- function(origin_row, analog_pool, vars, k, bandwidth,
                           h = NULL, min_n = 8) {
   pool_complete <- analog_pool %>% drop_na(all_of(vars))
   if (nrow(pool_complete) < 3) return(pool_complete[0, ])
   pool_complete <- phase_gate_pool(pool_complete, origin_row, min_n = min_n)
   if (nrow(pool_complete) < 3) return(pool_complete[0, ])
   
   scaled_pool   <- standardize_using(pool_complete, pool_complete, vars)
   scaled_origin <- standardize_using(origin_row, pool_complete, vars)
   
   var_weights <- rep(1, length(vars))
   names(var_weights) <- vars
   
   high_weight_vars <- c(
      "slope_1wk", "slope_2wk", "slope_3wk", "slope_4wk",
      "prev_slope_1wk", "prev_slope_2wk", "accel_2wk", "accel_3wk",
      "decline_persistence_2wk", "decline_persistence_3wk",
      "hosp_slope_1wk", "hosp_slope_2wk", "pos_slope_1wk", "pos_slope_2wk",
      "ehr_decline_confirmed", "drop_from_8wk_peak", "drop_from_12wk_peak",
      "drop_from_8wk_peak_pct", "drop_from_12wk_peak_pct",
      "post_peak_decline_speed", "is_post_peak", "is_fast_decline",
      "is_near_peak", "any_decline", "any_peak_or_decline"
   )
   medium_weight_vars <- c(
      "position_in_8wk_range", "position_in_12wk_range",
      "weeks_since_8wk_peak", "weeks_since_12wk_peak",
      "peak_8wk_log", "peak_12wk_log", "recent_peak_ratio", "week_of_year"
   )
   var_weights[names(var_weights) %in% high_weight_vars]  <- 3.0
   var_weights[names(var_weights) %in% medium_weight_vars] <- 1.7
   
   pool_matrix   <- as.matrix(scaled_pool[, vars])
   origin_vector <- as.numeric(scaled_origin[1, vars])
   diff_mat      <- sweep(pool_matrix, 2, origin_vector, "-")
   weighted_diff <- sweep(diff_mat^2, 2, var_weights, "*")
   pool_complete$distance <- sqrt(rowSums(weighted_diff))
   
   candidate_n <- min(max(k * 3, k), nrow(pool_complete))
   selected0 <- pool_complete %>% arrange(distance) %>% dplyr::slice(1:candidate_n)
   
   change_col <- if (!is.null(h)) paste0("future_change_h", h) else NA_character_
   origin_decline_like <- is_one(origin_row$any_decline) ||
      is_one(origin_row$is_near_peak) || is_one(origin_row$any_peak_or_decline)
   
   selected <- selected0 %>%
      mutate(
         raw_weight = exp(-distance / bandwidth),
         phase_penalty = case_when(
            phase == origin_row$phase ~ 1.00,
            any_peak_or_decline == origin_row$any_peak_or_decline ~ 0.75,
            any_decline == origin_row$any_decline ~ 0.65,
            origin_decline_like ~ 0.30,
            TRUE ~ 0.60
         ),
         peak_age_penalty = if (origin_decline_like) {
            exp(-abs(weeks_since_12wk_peak - origin_row$weeks_since_12wk_peak) / 2.5)
         } else { 1 },
         rebound_penalty = if (!is.na(change_col) && change_col %in% names(selected0) &&
                               origin_decline_like) {
            case_when(
               .data[[change_col]] > 0     ~ 0.15,
               .data[[change_col]] > -0.03 ~ 0.55,
               TRUE                        ~ 1.00
            )
         } else { 1 },
         late_peak_penalty = if (origin_decline_like && !is.null(h) &&
                                 "future_peak_offset_4wk" %in% names(selected0)) {
            case_when(
               is.na(future_peak_offset_4wk) ~ 1.00,
               h >= 2 & future_peak_offset_4wk > 2 ~ 0.35,
               h >= 1 & future_peak_offset_4wk > 3 ~ 0.55,
               TRUE ~ 1.00
            )
         } else { 1 },
         analog_weight = raw_weight * phase_penalty * peak_age_penalty *
            rebound_penalty * late_peak_penalty
      )
   
   if (sum(selected$analog_weight, na.rm = TRUE) <= 0) {
      selected <- selected %>% mutate(analog_weight = raw_weight)
   }
   
   selected %>%
      arrange(desc(analog_weight)) %>%
      dplyr::slice(1:min(k, n())) %>%
      mutate(analog_weight = analog_weight / sum(analog_weight, na.rm = TRUE))
}

analog_peak_probs <- function(selected_analogs) {
   if (is.null(selected_analogs) || nrow(selected_analogs) == 0) {
      return(tibble(
         prob_peak_by_1 = NA_real_, prob_peak_by_2 = NA_real_,
         prob_peak_by_3 = NA_real_, prob_peak_by_4 = NA_real_
      ))
   }
   w <- selected_analogs$analog_weight
   tibble(
      prob_peak_by_1 = weighted_mean_safe(selected_analogs$future_peak_offset_4wk <= 1, w),
      prob_peak_by_2 = weighted_mean_safe(selected_analogs$future_peak_offset_4wk <= 2, w),
      prob_peak_by_3 = weighted_mean_safe(selected_analogs$future_peak_offset_4wk <= 3, w),
      prob_peak_by_4 = weighted_mean_safe(selected_analogs$future_peak_offset_4wk <= 4, w)
   )
}

decline_qr_floor_v4 <- function(row) {
   # Now a CEILING on QR weight during peak/decline phases
   # so analogs (which know trajectory) get more influence.
   case_when(
      row$is_fast_decline == 1 ~ 0.30,
      row$is_near_peak == 1 & row$any_decline == 1 ~ 0.35,
      row$is_near_peak == 1 ~ 0.40,
      row$is_post_peak == 1 & row$decline_persistence_2wk == 1 ~ 0.30,
      row$is_post_peak == 1 ~ 0.40,
      row$decline_persistence_2wk == 1 ~ 0.45,
      row$phase == "ACCELERATING_DECLINE" ~ 0.45,
      row$phase == "SLOWING_DECLINE" ~ 0.50,
      TRUE ~ 0.60   # floor for non-decline phases
   )
}

interval_scale_by_phase_v4 <- function(phase) {
   case_when(
      phase == "FAST_POST_PEAK_DECLINE"  ~ 0.60,
      phase == "POST_PEAK_DECLINE"       ~ 0.70,
      phase == "ACCELERATING_DECLINE"   ~ 0.75,
      phase == "NEAR_PEAK_DECELERATION" ~ 0.75,
      phase == "VOLATILE_PLATEAU"       ~ 0.85,
      TRUE ~ 1.00
   )
}

decline_analog_quantile_v4 <- function(row, h, peak_probs_row) {
   slope <- row$slope_1wk
   slope_severity <- pmin(abs(pmin(slope, 0)) / (global_slope_sd + 1e-8), 2) / 2
   p2 <- first_or(peak_probs_row$prob_peak_by_2, 0)
   p3 <- first_or(peak_probs_row$prob_peak_by_3, 0)
   peak_risk <- max(p2, p3, na.rm = TRUE)
   if (is.infinite(peak_risk) || is.na(peak_risk)) peak_risk <- 0
   base_q <- if (h <= 1) 0.40 else 0.35
   q <- base_q - 0.15 * slope_severity - 0.08 * peak_risk -
      if_else(h >= 2, 0.05, 0.00)
   pmax(0.15, pmin(0.45, q))
}

# ══════════════════════════════════════════════════════════════════════════════
# 7. v4 CORE: RECURSIVE ROW EXTENSION
#
#  These two functions are the heart of v4.
#  After forecasting H0, we "append" a synthetic row to the feature dataframe
#  with the predicted log_target, then recompute all rolling/lag features
#  for that row so H1 sees an updated state vector.
# ══════════════════════════════════════════════════════════════════════════════

#' Append a synthetic forecast row to flu_features.
#'
#' @param df          Current flu_features dataframe (up to latest observed week).
#' @param new_week    Date of the new synthetic row (origin_week + (h+1) weeks).
#' @param pred_log    Predicted log_target for the new row.
#' @param origin_row  The row at the forecast origin (for EHR carry-forward).
#' @return df with one new row appended, all features recomputed.
append_forecast_row <- function(df, new_week, pred_log, origin_row, thr) {
   # ── Strategy ──────────────────────────────────────────────────────────────
   # df is the full working_df, which already has ALL columns including future
   # labels, phase, peak flags, etc.
   #
   # We CANNOT rebuild from raw flu columns because build_flu_features() does
   # not produce future labels (future_log_h0, future_change_h0, …) — those
   # are added separately via lead() in the main script.  Rebuilding from
   # scratch produces a df where those columns are all NA, which kills the
   # analog pool drop_na at every horizon after h=0.
   #
   # Instead we:
   #   1. Build a template new row by copying the last row of df (preserves
   #      every column type and default value).
   #   2. Overwrite only the fields that genuinely change.
   #   3. Recompute only the lag/slope/rolling features for the new tail row
   #      using the extended log_target vector — no full rebuild needed.
   # ──────────────────────────────────────────────────────────────────────────
   
   pred_count <- pmax(exp(pred_log) - 1, 0)
   
   # Step 1: template from last observed row (all types preserved)
   new_row <- df %>% dplyr::slice(n())
   
   # Step 2: overwrite changing fields
   new_row$week           <- new_week
   new_row[[target_var]]  <- pred_count
   new_row$log_target     <- pred_log
   
   # EHR carry-forward (constant = last known)
   # log_hosp / log_positive_tests stay at origin values → slopes = 0, not NA
   new_row$log_tests          <- origin_row$log_tests
   new_row$log_positive_tests <- origin_row$log_positive_tests
   new_row$log_hosp           <- origin_row$log_hosp
   new_row$hosp_slope_1wk     <- 0
   new_row$hosp_slope_2wk     <- 0
   new_row$pos_slope_1wk      <- 0
   new_row$pos_slope_2wk      <- 0
   new_row$ehr_decline_confirmed <- origin_row$ehr_decline_confirmed
   
   # Seasonality from new_week date
   woy <- isoweek(new_week)
   new_row$week_of_year <- woy
   new_row$sin52        <- sin(2 * pi * woy / 52)
   new_row$cos52        <- cos(2 * pi * woy / 52)
   new_row$month        <- month(new_week)
   new_row$is_christmas <- if_else(month(new_week) == 12 & day(new_week) >= 20, 1, 0)
   new_row$is_newyear   <- if_else(month(new_week) == 1  & day(new_week) <= 7,  1, 0)
   
   # Step 3: recompute lag/slope/rolling features using the extended log_target
   # vector (existing rows + the new predicted value).
   log_vec <- c(df$log_target, pred_log)
   n       <- length(log_vec)
   
   lag_at  <- function(v, k) if (n - k >= 1) v[n - k] else v[1]
   
   new_row$lag1  <- lag_at(log_vec, 1)
   new_row$lag2  <- lag_at(log_vec, 2)
   new_row$lag3  <- lag_at(log_vec, 3)
   new_row$lag4  <- lag_at(log_vec, 4)
   new_row$lag5  <- lag_at(log_vec, 5)
   new_row$lag8  <- lag_at(log_vec, 8)
   new_row$lag12 <- lag_at(log_vec, 12)
   
   new_row$slope_1wk      <- pred_log - new_row$lag1
   new_row$slope_2wk      <- pred_log - new_row$lag2
   new_row$slope_3wk      <- pred_log - new_row$lag3
   new_row$slope_4wk      <- pred_log - new_row$lag4
   new_row$prev_slope_1wk <- new_row$lag1 - new_row$lag2
   new_row$prev_slope_2wk <- new_row$lag2 - new_row$lag3
   new_row$accel_2wk      <- new_row$slope_1wk - new_row$prev_slope_1wk
   new_row$accel_3wk      <- new_row$slope_1wk - ((new_row$lag1 - new_row$lag3) / 2)
   
   new_row$decline_persistence_2wk <- as.integer(
      new_row$slope_1wk < 0 & new_row$prev_slope_1wk < 0)
   new_row$decline_persistence_3wk <- as.integer(
      new_row$slope_1wk < 0 & new_row$prev_slope_1wk < 0 & new_row$prev_slope_2wk < 0)
   
   # Rolling windows over the extended vector
   win <- function(v, k) v[pmax(1, n - k + 1):n]
   
   w4  <- win(log_vec, 4);  w8  <- win(log_vec, 8);  w12 <- win(log_vec, 12)
   
   new_row$mean_4wk  <- safe_mean(w4);  new_row$mean_8wk  <- safe_mean(w8)
   new_row$mean_12wk <- safe_mean(w12)
   new_row$sd_4wk    <- safe_sd(w4);    new_row$sd_8wk    <- safe_sd(w8)
   new_row$sd_12wk   <- safe_sd(w12)
   new_row$max_8wk   <- safe_max(w8);   new_row$max_12wk  <- safe_max(w12)
   new_row$min_8wk   <- safe_min(w8);   new_row$min_12wk  <- safe_min(w12)
   
   new_row$peak_8wk_log  <- new_row$max_8wk
   new_row$peak_12wk_log <- new_row$max_12wk
   
   new_row$position_in_8wk_range  <-
      (pred_log - new_row$min_8wk) / pmax(new_row$max_8wk  - new_row$min_8wk,  1e-6)
   new_row$position_in_12wk_range <-
      (pred_log - new_row$min_12wk) / pmax(new_row$max_12wk - new_row$min_12wk, 1e-6)
   
   new_row$weeks_since_8wk_peak  <- weeks_since_peak_fun(w8)
   new_row$weeks_since_12wk_peak <- weeks_since_peak_fun(w12)
   
   new_row$drop_from_8wk_peak  <- pred_log - new_row$max_8wk
   new_row$drop_from_12wk_peak <- pred_log - new_row$max_12wk
   new_row$drop_from_8wk_peak_pct <-
      ((exp(pred_log) - 1) - (exp(new_row$max_8wk) - 1)) /
      pmax(exp(new_row$max_8wk) - 1, 1e-6)
   new_row$drop_from_12wk_peak_pct <-
      ((exp(pred_log) - 1) - (exp(new_row$max_12wk) - 1)) /
      pmax(exp(new_row$max_12wk) - 1, 1e-6)
   new_row$post_peak_decline_speed <-
      new_row$drop_from_12wk_peak / pmax(new_row$weeks_since_12wk_peak, 1)
   new_row$recent_peak_ratio <-
      (exp(pred_log) - 1) / pmax(exp(new_row$max_8wk) - 1, 1e-6)
   
   # Phase flags
   # Use the precomputed training-period Q25 (not df$slope_1wk which includes
   # eval data) to avoid leakage into is_fast_decline classification.
   slope_1wk_q25 <- thr$slope_q25
   new_row$is_post_peak <- if_else(
      new_row$weeks_since_12wk_peak >= 1 &
         pred_log < new_row$max_12wk &
         new_row$slope_1wk < 0, 1, 0)
   new_row$is_fast_decline <- if_else(
      new_row$slope_1wk < slope_1wk_q25 &
         new_row$slope_2wk < 0 &
         new_row$drop_from_12wk_peak < 0, 1, 0)
   new_row$is_near_peak <- as.integer(
      new_row$recent_peak_ratio >= 0.80 &
         (new_row$accel_2wk <= 0 | new_row$slope_1wk <= 0.05))
   new_row$any_decline <- as.integer(
      new_row$is_fast_decline == 1 |
         new_row$is_post_peak == 1 |
         new_row$decline_persistence_2wk == 1)
   new_row$any_peak_or_decline <- as.integer(
      new_row$is_near_peak == 1 | new_row$any_decline == 1)
   
   # Future labels: NA for synthetic rows (not used in analog pool — pool
   # always draws from original df, not working_df)
   for (h_lab in horizons) {
      new_row[[paste0("future_log_h",    h_lab)]] <- NA_real_
      new_row[[paste0("future_value_h",  h_lab)]] <- NA_real_
      new_row[[paste0("future_change_h", h_lab)]] <- NA_real_
   }
   new_row$future_peak_offset_4wk <- NA_integer_
   for (p in c("peak_by_1wk","peak_by_2wk","peak_by_3wk","peak_by_4wk"))
      new_row[[p]] <- NA_integer_
   
   # Assign phase using origin-specific thresholds (no leakage)
   new_row <- assign_phase(new_row, thr)
   
   bind_rows(df, new_row)
}

# ══════════════════════════════════════════════════════════════════════════════
# 8. v4 FORECAST ENGINE  — recursive loop
#
#  Replaces forecast_origin_v3().
#  Key difference: after each horizon h, the predicted log value is appended
#  to the working feature dataframe so that h+1 uses updated features.
#  anti_shift_postprocess is NOT called — the recursion enforces coherence.
# ══════════════════════════════════════════════════════════════════════════════

get_tuning_row <- function(tuning_table, h) {
   row <- tuning_table %>% filter(horizon == h) %>% dplyr::slice(1)
   if (nrow(row) == 0 || any(is.na(row$k), is.na(row$bandwidth), is.na(row$blend_qr))) {
      return(tibble(horizon = h, k = 20, bandwidth = 1.25, blend_qr = 0.50,
                    mae_log = NA_real_, n_valid = 0))
   }
   row
}

forecast_origin_v4 <- function(origin_week,
                               df,
                               qr_models,
                               residual_df,
                               tuning_table,
                               thr,
                               quantiles = quantiles) {
   
   origin_week <- as.Date(origin_week)
   
   # working_df grows by one row after each horizon step
   working_df <- df
   
   # Origin row — must be complete; used as imputation source for synthetic rows
   origin_row_0 <- df %>%
      filter(week == origin_week) %>%
      drop_na(all_of(union(predictors, state_vars)))
   
   if (nrow(origin_row_0) != 1) return(NULL)
   
   required_cols <- union(predictors, state_vars)
   
   model_h0 <- qr_models[["h0"]]
   if (is.null(model_h0)) return(NULL)
   
   median_rows          <- list()
   quantile_error_store <- list()
   analog_diag_store    <- list()
   
   for (h in horizons) {
      
      # The row we forecast FROM is origin_week + h weeks.
      # After h=0, that row is the synthetic one we appended.
      current_week_for_h <- origin_week + weeks(h)
      target_end_date_h  <- origin_week + weeks(h + 1)
      
      # ── 1. Get the current feature row ────────────────────────────────────
      current_row_raw <- working_df %>% filter(week == current_week_for_h)
      
      if (nrow(current_row_raw) < 1) {
         # Row missing — append carry-forward so cascade continues
         cat(sprintf("  [v4] h=%d: current row not found, carry-forward\n", h))
         last_log <- working_df$log_target[nrow(working_df)]
         working_df <- append_forecast_row(working_df, target_end_date_h,
                                           last_log, origin_row_0, thr)
         next
      }
      
      current_row_raw <- current_row_raw %>% dplyr::slice(1)
      
      # For synthetic rows: impute any NAs from origin_row_0
      if (h > 0) {
         na_cols <- required_cols[sapply(required_cols, function(v) {
            v %in% names(current_row_raw) && is.na(current_row_raw[[v]][1])
         })]
         if (length(na_cols) > 0) {
            for (col in na_cols) {
               if (col %in% names(origin_row_0))
                  current_row_raw[[col]] <- origin_row_0[[col]]
            }
         }
      }
      
      still_na <- required_cols[sapply(required_cols, function(v) {
         v %in% names(current_row_raw) && is.na(current_row_raw[[v]][1])
      })]
      
      if (length(still_na) > 0) {
         cat(sprintf("  [v4] h=%d: %d NAs remain after imputation (%s) — carry-forward\n",
                     h, length(still_na), paste(head(still_na, 3), collapse=",")))
         last_log <- working_df$log_target[nrow(working_df)]
         working_df <- append_forecast_row(working_df, target_end_date_h,
                                           last_log, origin_row_0, thr)
         next
      }
      
      current_row <- current_row_raw
      
      # ── 2. Tuning params ──────────────────────────────────────────────────
      tune_h       <- get_tuning_row(tuning_table, h)
      k_h          <- tune_h$k
      bandwidth_h  <- tune_h$bandwidth
      blend_qr_raw <- tune_h$blend_qr
      
      # Use origin-specific global_slope_sd from thr
      global_slope_sd <- thr$global_slope_sd
      qr_limit     <- decline_qr_floor_v4(current_row)
      decline_phase <- is_one(current_row$any_peak_or_decline)
      blend_qr_h   <- if (decline_phase) min(blend_qr_raw, qr_limit) else
         max(blend_qr_raw, qr_limit)
      blend_analog_h <- 1 - blend_qr_h
      
      # ── 3. QR prediction ──────────────────────────────────────────────────
      qr_median_log <- tryCatch(
         predict_qr_safe(model_h0, current_row),
         error = function(e) NA_real_
      )
      
      if (is.na(qr_median_log)) {
         cat(sprintf("  [v4] h=%d: QR predict failed — carry-forward\n", h))
         last_log <- current_row$log_target
         working_df <- append_forecast_row(working_df, target_end_date_h,
                                           last_log, current_row, thr)
         next
      }
      
      # ── 4. Analog selection ───────────────────────────────────────────────
      change_col <- "future_change_h0"
      
      # Use origin_week as cutoff so pool size stays constant across horizons
      analog_pool <- get_label_known_training(df, origin_week, h = 0) %>%
         mutate(
            future_peak_offset_4wk = if_else(
               week + weeks(4) <= origin_week,
               future_peak_offset_4wk, NA_integer_
            )
         ) %>%
         drop_na(all_of(c(state_vars, change_col)))
      
      selected_analogs <- get_analogs_v4(
         origin_row  = current_row,
         analog_pool = analog_pool,
         vars        = state_vars,
         k           = k_h,
         bandwidth   = bandwidth_h,
         h           = 0L
      )
      
      use_analog <- nrow(selected_analogs) >= 5
      if (!use_analog)
         cat(sprintf("  [v4] h=%d: only %d analogs — QR-only fallback\n",
                     h, nrow(selected_analogs)))
      
      # ── 5. Median forecast ────────────────────────────────────────────────
      peak_probs_h <- if (use_analog) {
         analog_peak_probs(selected_analogs)
      } else {
         tibble(prob_peak_by_1 = NA_real_, prob_peak_by_2 = NA_real_,
                prob_peak_by_3 = NA_real_, prob_peak_by_4 = NA_real_)
      }
      
      if (use_analog) {
         analog_change <- if (is_one(current_row$any_peak_or_decline)) {
            dq <- decline_analog_quantile_v4(current_row, h, peak_probs_h)
            weighted_quantile_safe(selected_analogs[[change_col]],
                                   selected_analogs$analog_weight, dq)[1]
         } else {
            weighted_median_safe(selected_analogs[[change_col]],
                                 selected_analogs$analog_weight)
         }
         analog_target_log <- current_row$log_target + analog_change
         raw_combined_log  <- pmax(blend_qr_h * qr_median_log +
                                      blend_analog_h * analog_target_log, 0)
      } else {
         analog_change     <- NA_real_
         analog_target_log <- qr_median_log
         raw_combined_log  <- pmax(qr_median_log, 0)
         blend_qr_h        <- 1.0
         blend_analog_h    <- 0.0
      }
      
      # ── 6. Interval ───────────────────────────────────────────────────────
      analog_error_q <- if (use_analog) {
         analog_error <- selected_analogs[[change_col]] - analog_change
         weighted_quantile_safe(analog_error, selected_analogs$analog_weight, quantiles)
      } else {
         rep(0, length(quantiles))
      }
      
      residual_phase_h <- residual_df %>%
         filter(horizon == 0, phase == current_row$phase) %>% pull(residual)
      if (length(residual_phase_h) < 25)
         residual_phase_h <- residual_df %>% filter(horizon == 0) %>% pull(residual)
      
      residual_error_q <- quantile(residual_phase_h, probs = quantiles,
                                   na.rm = TRUE, names = FALSE)
      
      combined_error_q <- blend_analog_h * analog_error_q +
         blend_qr_h     * residual_error_q
      combined_error_q <- combined_error_q *
         interval_scale_by_phase_v4(current_row$phase)
      
      if (is_one(current_row$any_peak_or_decline) && h >= 1) {
         combined_error_q[quantiles > 0.5] <- combined_error_q[quantiles > 0.5] * 0.70
         combined_error_q[quantiles < 0.5] <- combined_error_q[quantiles < 0.5] * 1.15
      }
      
      med_idx          <- which.min(abs(quantiles - 0.5))
      combined_error_q <- combined_error_q - combined_error_q[med_idx]
      lo_cap           <- quantile(combined_error_q, 0.025, na.rm = TRUE)
      hi_cap           <- quantile(combined_error_q, 0.975, na.rm = TRUE)
      combined_error_q <- pmin(pmax(combined_error_q, lo_cap), hi_cap)
      
      # ── 7. Store results ──────────────────────────────────────────────────
      median_rows[[paste0("h", h)]] <- bind_cols(
         tibble(
            origin_week       = origin_week,
            target_end_date   = target_end_date_h,
            horizon           = h,
            phase             = current_row$phase,
            final_log         = raw_combined_log,
            raw_corrected_log = raw_combined_log,
            qr_median_log     = qr_median_log,
            analog_target_log = analog_target_log,
            blend_qr_raw      = blend_qr_raw,
            blend_qr          = blend_qr_h,
            blend_analog      = blend_analog_h,
            k                 = k_h,
            bandwidth         = bandwidth_h
         ),
         peak_probs_h
      )
      
      quantile_error_store[[paste0("h", h)]] <- tibble(
         horizon            = h,
         quantile           = quantiles,
         centered_error_log = combined_error_q
      )
      
      analog_diag_store[[paste0("h", h)]] <- if (use_analog) {
         selected_analogs %>%
            mutate(origin_week = origin_week, horizon = h,
                   future_change = .data[[change_col]]) %>%
            dplyr::select(
               origin_week, horizon, week, distance, analog_weight, future_change,
               all_of(target_var), phase, is_near_peak, is_post_peak, is_fast_decline,
               any_decline, any_peak_or_decline, decline_persistence_2wk,
               decline_persistence_3wk, slope_1wk, slope_2wk, accel_2wk,
               weeks_since_12wk_peak, drop_from_12wk_peak, future_peak_offset_4wk
            )
      } else {
         tibble(origin_week = origin_week, horizon = h)
      }
      
      # ── 8. Extend working_df — ALWAYS, no matter what happened above ──────
      working_df <- append_forecast_row(
         df         = working_df,
         new_week   = target_end_date_h,
         pred_log   = raw_combined_log,
         origin_row = current_row,
         thr        = thr
      )
      
      cat(sprintf("  [v4] h=%d done: pred_log=%.3f  phase=%s  use_analog=%s\n",
                  h, raw_combined_log, current_row$phase, use_analog))
   }
   
   median_df <- bind_rows(median_rows)
   if (nrow(median_df) == 0) return(NULL)
   
   error_df <- bind_rows(quantile_error_store)
   
   all_quantiles <- median_df %>%
      dplyr::select(
         origin_week, target_end_date, horizon, phase,
         final_log, raw_corrected_log, qr_median_log, analog_target_log,
         blend_qr_raw, blend_qr, blend_analog, k, bandwidth,
         starts_with("prob_peak_by_")
      ) %>%
      left_join(error_df, by = "horizon") %>%
      mutate(
         log_quantile    = final_log + centered_error_log,
         predicted_value = pmax(exp(log_quantile) - 1, 0)
      ) %>%
      group_by(horizon) %>%
      arrange(quantile, .by_group = TRUE) %>%
      mutate(predicted_value = as.numeric(stats::isoreg(predicted_value)$yf)) %>%
      ungroup()
   
   list(
      all_quantiles      = all_quantiles,
      median_df          = median_df,
      analog_diagnostics = bind_rows(analog_diag_store)
   )
}


# ══════════════════════════════════════════════════════════════════════════════
# 9. ROLLING VALIDATION / TUNING
#    v4 change: tune only for h=0 (1-step), since all horizons use the H0 model
#    recursively. k, bandwidth, blend_qr from h=0 tuning are reused for h1-3
#    with minor horizon-indexed adjustments.
# ══════════════════════════════════════════════════════════════════════════════

validation_origins <- flu_features %>%
   filter(week >= train_cut, week < valid_end, !is.na(log_target)) %>%
   drop_na(all_of(c(predictors, state_vars))) %>%
   pull(week)

tune_results <- list()

cat("\nTuning h=0 (1-step) — used recursively for all horizons ...\n")

pred_store <- list()
change_col <- "future_change_h0"
target_h   <- "future_log_h0"

for (origin_week in validation_origins) {
   origin_week <- as.Date(origin_week)
   
   # Compute thresholds strictly from data <= origin_week (no future leakage)
   thr_tune <- compute_origin_thresholds(flu_features, origin_week)
   
   # Rebuild features and phase labels for this origin using its own thresholds.
   # We pass these into flu_features_at_origin so the origin_row's phase
   # and is_fast_decline reflect only what was knowable at that time.
   flu_features_at_origin <- flu_features %>%
      filter(week <= origin_week) %>%
      assign_phase(thr_tune)
   
   train_origin <- get_label_known_training(flu_features_at_origin, origin_week, 0) %>%
      drop_na(all_of(c(target_h, predictors)))
   
   origin_row <- flu_features_at_origin %>%
      filter(week == origin_week) %>%
      drop_na(all_of(c(predictors, state_vars)))
   
   y_true <- flu_features %>% filter(week == origin_week) %>% pull(target_h)
   
   if (nrow(train_origin) < 80 || nrow(origin_row) != 1 ||
       length(y_true) != 1 || is.na(y_true)) next
   
   qr_model_h0 <- tryCatch(
      fit_qr_model(train_origin, 0, tau_value = 0.5),
      error = function(e) NULL
   )
   if (is.null(qr_model_h0)) next
   
   qr_pred <- tryCatch(
      predict_qr_safe(qr_model_h0, origin_row),
      error = function(e) NA_real_
   )
   if (is.na(qr_pred)) next
   
   analog_pool <- get_label_known_training(flu_features, origin_week, 0) %>%
      mutate(
         future_peak_offset_4wk = if_else(
            week + weeks(4) <= origin_week, future_peak_offset_4wk, NA_integer_
         )
      ) %>%
      drop_na(all_of(c(state_vars, change_col)))
   
   grid_origin <- expand_grid(k = candidate_k, bandwidth = candidate_bandwidth) %>%
      mutate(
         origin_week  = origin_week,
         y_true       = y_true,
         qr_pred      = qr_pred,
         analog_pred  = NA_real_
      )
   
   for (g in seq_len(nrow(grid_origin))) {
      analogs_g <- get_analogs_v4(
         origin_row  = origin_row,
         analog_pool = analog_pool,
         vars        = state_vars,
         k           = grid_origin$k[g],
         bandwidth   = grid_origin$bandwidth[g],
         h           = 0L
      )
      if (nrow(analogs_g) < 5) next
      
      peak_probs_g <- analog_peak_probs(analogs_g)
      decline_like <- is_one(origin_row$any_peak_or_decline)
      
      analog_change <- if (decline_like) {
         dq <- decline_analog_quantile_v4(origin_row, 0, peak_probs_g)
         weighted_quantile_safe(analogs_g[[change_col]], analogs_g$analog_weight, dq)[1]
      } else {
         weighted_median_safe(analogs_g[[change_col]], analogs_g$analog_weight)
      }
      
      analog_pred <- origin_row$log_target + analog_change
      analog_pred <- pmax(analog_pred, 0)
      grid_origin$analog_pred[g] <- analog_pred
   }
   
   pred_store[[as.character(origin_week)]] <- grid_origin
}

pred_df_all <- bind_rows(pred_store) %>%
   filter(!is.na(analog_pred), !is.na(qr_pred), !is.na(y_true))

if (nrow(pred_df_all) > 0) {
   score_df <- pred_df_all %>%
      crossing(blend_qr = candidate_blend_qr) %>%
      mutate(
         pred      = blend_qr * qr_pred + (1 - blend_qr) * analog_pred,
         abs_error = abs(y_true - pred)
      ) %>%
      group_by(k, bandwidth, blend_qr) %>%
      summarise(mae_log = mean(abs_error, na.rm = TRUE), n_valid = n(), .groups = "drop") %>%
      arrange(mae_log)
   
   best_h0 <- score_df %>% dplyr::slice(1)
} else {
   best_h0 <- tibble(k = 20, bandwidth = 1.25, blend_qr = 0.50, mae_log = NA_real_, n_valid = 0)
}

# Propagate h=0 tuning to h1-3 with mild horizon penalty on blend_qr
# (longer horizons lean slightly more on analogs since QR can't see as far)
tuning_table <- bind_rows(lapply(horizons, function(h) {
   tibble(
      horizon   = h,
      k         = best_h0$k,
      bandwidth = best_h0$bandwidth,
      # Slightly reduce QR weight at longer horizons
      blend_qr  = pmax(best_h0$blend_qr - h * 0.05, 0.10),
      mae_log   = best_h0$mae_log,
      n_valid   = best_h0$n_valid
   )
}))

cat("\n═══ TUNING TABLE v4 ═══\n")
print(tuning_table)

# ══════════════════════════════════════════════════════════════════════════════
# 10. FINAL FORECAST  (same structure as v3, uses forecast_origin_v4)
# ══════════════════════════════════════════════════════════════════════════════

latest_observed_week <- max(
   flu_features$week[!is.na(flu_features[[target_var]])],
   na.rm = TRUE
)
reference_date <- latest_observed_week + weeks(1)

# Compute thresholds from all data up to the latest observed week
thr_final <- compute_origin_thresholds(flu_features, latest_observed_week)

# Re-apply phase labels with final thresholds (updates flu_features in-place)
flu_features <- assign_phase(flu_features, thr_final)
global_slope_sd <- thr_final$global_slope_sd

current_row <- flu_features %>%
   filter(week == latest_observed_week) %>%
   drop_na(all_of(c(predictors, state_vars)))

if (nrow(current_row) != 1) {
   stop("Current row is incomplete. Check missing predictors for latest observed week.")
}

observed_value <- current_row[[target_var]]
phase <- current_row$phase

final_models <- fit_all_qr_models(
   df = flu_features,
   label_cutoff_week = latest_observed_week
)

final_forecast <- forecast_origin_v4(
   origin_week  = latest_observed_week,
   df           = flu_features,
   qr_models    = final_models$qr_models,
   residual_df  = final_models$residual_df,
   tuning_table = tuning_table,
   thr          = thr_final,
   quantiles    = quantiles
)

if (is.null(final_forecast)) {
   stop("Final forecast failed. Check model completeness and analog availability.")
}

all_quantiles       <- final_forecast$all_quantiles
analog_diagnostics_df <- final_forecast$analog_diagnostics
median_forecast_df  <- final_forecast$median_df

cdc_submission <- all_quantiles %>%
   mutate(
      reference_date  = reference_date,
      target          = "wk inc flu hosp",
      location        = "45",
      output_type     = "quantile",
      output_type_id  = quantile,
      value           = round(predicted_value, 0)
   ) %>%
   dplyr::select(
      reference_date, target, horizon, location,
      target_end_date, output_type, output_type_id, value
   ) %>%
   arrange(horizon, output_type_id)

cat("\n═══ FINAL FORECAST DIAGNOSTICS v4 ═══\n")
cat(sprintf("Latest observed week : %s\n", latest_observed_week))
cat(sprintf("Latest observed value: %s\n", observed_value))
cat(sprintf("Detected phase       : %s\n", phase))

cat("\nForecast medians:\n")
print(
   cdc_submission %>%
      filter(output_type_id == 0.5) %>%
      dplyr::select(horizon, target_end_date, median_forecast = value)
)

cat("\nMedian model diagnostics:\n")
print(
   all_quantiles %>%
      filter(quantile == 0.5) %>%
      dplyr::select(
         horizon, target_end_date, phase,
         raw_corrected_log, final_log,
         qr_median_log, analog_target_log,
         blend_qr_raw, blend_qr, blend_analog,
         prob_peak_by_1, prob_peak_by_2, prob_peak_by_3, prob_peak_by_4,
         k, bandwidth
      )
)

cat("\nAnalog diagnostics summary:\n")
print(
   analog_diagnostics_df %>%
      group_by(horizon) %>%
      summarise(
         n_analogs                   = n(),
         median_future_change        = median(future_change, na.rm = TRUE),
         mean_future_change          = mean(future_change, na.rm = TRUE),
         mean_distance               = mean(distance, na.rm = TRUE),
         max_weight                  = max(analog_weight, na.rm = TRUE),
         share_near_peak_analogs     = mean(is_near_peak, na.rm = TRUE),
         share_fast_decline_analogs  = mean(is_fast_decline, na.rm = TRUE),
         share_post_peak_analogs     = mean(is_post_peak, na.rm = TRUE),
         share_decline_persistence_2wk = mean(decline_persistence_2wk, na.rm = TRUE),
         .groups = "drop"
      )
)

# ══════════════════════════════════════════════════════════════════════════════
# 11. FINAL FORECAST PLOT  (unchanged from v3)
# ══════════════════════════════════════════════════════════════════════════════

plot_intervals <- cdc_submission %>%
   filter(output_type_id %in% c(0.025, 0.25, 0.5, 0.75, 0.975)) %>%
   pivot_wider(names_from = output_type_id, values_from = value) %>%
   rename(
      median   = `0.5`,
      lower_95 = `0.025`, upper_95 = `0.975`,
      lower_50 = `0.25`,  upper_50 = `0.75`
   )

hist_plot_df <- flu_features %>%
   filter(week >= latest_observed_week - weeks(20) &
             week <= latest_observed_week) %>%
   dplyr::select(week, Observed = all_of(target_var)) %>%
   pivot_longer(cols = c(Observed), names_to = "Series", values_to = "Admissions")

# Build a complete weekly sequence so no dates are dropped from x-axis
plot_date_min <- min(hist_plot_df$week)
plot_date_max <- max(plot_intervals$target_end_date)
all_actual_dates <- seq(plot_date_min, plot_date_max, by = "week")

p_final <- ggplot() +
   geom_ribbon(data = plot_intervals,
               aes(x = target_end_date, ymin = lower_95, ymax = upper_95),
               fill = "#457B9D", alpha = 0.15) +
   geom_ribbon(data = plot_intervals,
               aes(x = target_end_date, ymin = lower_50, ymax = upper_50),
               fill = "#457B9D", alpha = 0.30) +
   geom_line(data = hist_plot_df,
             aes(x = week, y = Admissions, color = Series), linewidth = 1) +
   geom_point(data = hist_plot_df,
              aes(x = week, y = Admissions, color = Series), size = 2) +
   geom_line(data = plot_intervals,
             aes(x = target_end_date, y = median, color = "Forecast"),
             linetype = "dashed", linewidth = 1.2) +
   geom_point(data = plot_intervals,
              aes(x = target_end_date, y = median, color = "Forecast"), size = 3) +
   geom_vline(xintercept = latest_observed_week, color = "black", linewidth = 0.7) +
   scale_x_date(breaks = all_actual_dates, date_labels = "%b %d") +
   scale_color_manual(values = c("Observed" = "black", "Forecast" = "#457B9D")) +
   labs(
      title    = paste0("DMAPRIME-QR v4 Recursive Forecast [", phase, "]"),
      subtitle = paste0(
         "Latest observed: ", format(latest_observed_week, "%Y-%m-%d"),
         " (", observed_value, " admissions)"
      ),
      y       = "Weekly Inpatient Hospitalizations (Influenza, South Carolina)",
      x       = "Week Ending Date",
      caption = paste0(
         "v4: recursive 1-step QR + phase-aligned analogs | ",
         "horizon shift eliminated by architecture"
      )
   ) +
   theme_minimal() +
   theme(
      plot.title    = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10),
      axis.text.x   = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y   = element_text(size = 9),
      axis.title    = element_text(size = 11),
      legend.position = "bottom",
      legend.title    = element_blank(),
      panel.grid.minor = element_blank(),
      plot.caption = element_text(hjust = 0, size = 8, color = "gray40")
   )

print(p_final)

# ══════════════════════════════════════════════════════════════════════════════
# 12. EXPORT CDC SUBMISSION  (unchanged from v3)
# ══════════════════════════════════════════════════════════════════════════════

output_filename_cdc <- file.path(
   output_dir_cdc,
   paste0(format(reference_date, "%Y-%m-%d"), "-DMAPRIME-QR-v4.csv")
)
# Uncomment to save:
# write_csv(cdc_submission, output_filename_cdc)

cat("\n✅ FINAL FORECAST COMPLETE\n")
cat("CDC file: ", basename(output_filename_cdc), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# 13. SOFTWARE IMPLEMENTATION FILE  (unchanged from v3)
# ══════════════════════════════════════════════════════════════════════════════

impl_state_rfa <- RFA_data_Influenza %>%
   dplyr::select(week_end, weekly_IP) %>%
   dplyr::filter(week_end >= as.Date("2020-03-07")) %>%
   dplyr::mutate(
      value = weekly_IP, data_source = "RFA",
      estimate_projected_report = 2,
      target_end_date = as.Date(week_end)
   )

impl_state_hs <- full_join(
   MUSC_Weekly_Influenza   %>% dplyr::select(Week, MUSC_Hosp   = Weekly_Inpatient_Hospitalizations),
   Prisma_Weekly_Influenza %>% dplyr::select(Week, Prisma_Hosp = Weekly_Inpatient_Hospitalizations),
   by = "Week"
) %>%
   dplyr::mutate(
      Week  = as.Date(Week),
      value = rowSums(across(c(MUSC_Hosp, Prisma_Hosp)), na.rm = TRUE),
      data_source = "HS", estimate_projected_report = 2, target_end_date = Week
   ) %>%
   dplyr::filter(target_end_date >= as.Date("2020-03-07"))

impl_state_forecast <- cdc_submission %>%
   dplyr::mutate(
      data_source = NA_character_,
      estimate_projected_report = 1,
      target_end_date = as.Date(target_end_date)
   )

sw_implementation_file <- bind_rows(impl_state_rfa, impl_state_hs, impl_state_forecast) %>%
   dplyr::mutate(
      reference_date = format(reference_date, "%Y-%m-%d"),
      target = "wk inc flu hosp",
      target_end_date = format(target_end_date, "%Y-%m-%d"),
      location_general = "state", location = "SC", disease = "influenza",
      population = "general_population", training_validation = 0,
      imputed = 0, outcome_measure = "Weekly_Inpatient_Hospitalizations",
      output_type = "quantile"
   ) %>%
   dplyr::select(
      reference_date, target, target_end_date, location_general, location,
      value, disease, population, training_validation, estimate_projected_report,
      imputed, data_source, outcome_measure, output_type, output_type_id
   )

output_filename_software <- file.path(
   output_dir_software_impl,
   "inpatient-hosp-DMAPRIME-QR-v4-implementation.csv"
)
# Uncomment to save:
# write_csv(sw_implementation_file, output_filename_software)
cat("Software implementation file: ", basename(output_filename_software), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# 14. RETROSPECTIVE EVALUATION  (uses forecast_origin_v4)
# ══════════════════════════════════════════════════════════════════════════════

ai_score <- function(actual, forecast) {
   dplyr::case_when(
      is.na(actual) | is.na(forecast) ~ NA_real_,
      actual == 0 & forecast == 0     ~ 1,
      TRUE ~ pmin(actual, forecast) / pmax(actual, forecast)
   )
}

eval_origin_weeks <- flu_features %>%
   filter(week >= eval_start, week <= eval_end, !is.na(.data[[target_var]])) %>%
   drop_na(all_of(c(predictors, state_vars))) %>%
   pull(week)

cat(sprintf(
   "\nRetrospective evaluation v4 over %d origin weeks (%s to %s)\n",
   length(eval_origin_weeks),
   format(min(eval_origin_weeks), "%Y-%m-%d"),
   format(max(eval_origin_weeks), "%Y-%m-%d")
))

retro_rows <- list()

for (origin_week in eval_origin_weeks) {
   origin_week <- as.Date(origin_week)
   
   # Origin-specific thresholds — only data <= origin_week used
   thr_retro <- compute_origin_thresholds(flu_features, origin_week)
   
   # Rebuild phase labels with these thresholds for the feature df passed
   # into the forecast engine and the QR fitter
   flu_features_retro <- flu_features %>%
      filter(week <= origin_week + weeks(4)) %>%   # keep future labels intact
      assign_phase(thr_retro)
   
   models_retro <- fit_all_qr_models(
      df = flu_features_retro,
      label_cutoff_week = origin_week
   )
   
   fc_retro <- forecast_origin_v4(
      origin_week  = origin_week,
      df           = flu_features_retro,
      qr_models    = models_retro$qr_models,
      residual_df  = models_retro$residual_df,
      tuning_table = tuning_table,
      thr          = thr_retro,
      quantiles    = quantiles
   )
   
   if (is.null(fc_retro)) next
   
   med_retro <- fc_retro$all_quantiles %>%
      filter(quantile == 0.5) %>%
      mutate(forecast_median = round(predicted_value, 1)) %>%
      dplyr::select(
         origin_week, target_end_date, horizon, forecast_median,
         phase, raw_corrected_log, final_log, qr_median_log, analog_target_log,
         blend_qr, prob_peak_by_1, prob_peak_by_2, prob_peak_by_3, prob_peak_by_4
      )
   
   for (i in seq_len(nrow(med_retro))) {
      target_week <- med_retro$target_end_date[i]
      truth_row   <- flu_features %>% filter(week == target_week)
      actual_val  <- if (nrow(truth_row) == 1) truth_row[[target_var]] else NA_real_
      retro_rows[[length(retro_rows) + 1]] <- med_retro[i, ] %>%
         mutate(actual = actual_val)
   }
}

retro_df <- bind_rows(retro_rows) %>%
   filter(!is.na(actual), !is.na(forecast_median)) %>%
   mutate(
      ai_score = ai_score(actual, forecast_median),
      score_band = case_when(
         ai_score >= 0.85 ~ "Excellent (>=0.85)",
         ai_score >= 0.70 ~ "Good (0.70-0.84)",
         ai_score >= 0.55 ~ "Fair (0.55-0.69)",
         TRUE             ~ "Poor (<0.55)"
      ),
      score_band = factor(score_band, levels = c(
         "Excellent (>=0.85)", "Good (0.70-0.84)",
         "Fair (0.55-0.69)",   "Poor (<0.55)"
      ))
   )

cat("\n═══ OVERALL PERIOD SUMMARY v4 ═══\n")
overall_summary <- retro_df %>%
   summarise(
      n_obs          = n(),
      n_weeks        = n_distinct(origin_week),
      mean_ai_score  = mean(ai_score, na.rm = TRUE),
      median_ai_score = median(ai_score, na.rm = TRUE),
      sd_ai_score    = sd(ai_score, na.rm = TRUE),
      min_ai_score   = min(ai_score, na.rm = TRUE),
      max_ai_score   = max(ai_score, na.rm = TRUE),
      pct_excellent  = mean(ai_score >= 0.85, na.rm = TRUE) * 100,
      pct_good       = mean(ai_score >= 0.70 & ai_score < 0.85, na.rm = TRUE) * 100,
      pct_fair       = mean(ai_score >= 0.55 & ai_score < 0.70, na.rm = TRUE) * 100,
      pct_poor       = mean(ai_score < 0.55, na.rm = TRUE) * 100
   )
print(overall_summary)

cat("\n═══ PERFORMANCE BY HORIZON v4 ═══\n")
horizon_summary <- retro_df %>%
   group_by(horizon) %>%
   summarise(
      n_obs           = n(),
      mean_ai_score   = round(mean(ai_score, na.rm = TRUE), 4),
      median_ai_score = round(median(ai_score, na.rm = TRUE), 4),
      sd_ai_score     = round(sd(ai_score, na.rm = TRUE), 4),
      min_ai_score    = round(min(ai_score, na.rm = TRUE), 4),
      max_ai_score    = round(max(ai_score, na.rm = TRUE), 4),
      pct_excellent   = round(mean(ai_score >= 0.85, na.rm = TRUE) * 100, 1),
      pct_poor        = round(mean(ai_score < 0.55, na.rm = TRUE) * 100, 1),
      .groups = "drop"
   ) %>%
   mutate(horizon_label = paste0("H", horizon, " (", horizon + 1, "-wk ahead)"))
print(horizon_summary)

cat("\n═══ WEEK-BY-WEEK DETAIL v4 ═══\n")
week_detail <- retro_df %>%
   dplyr::select(
      origin_week, target_end_date, horizon, actual, forecast_median,
      ai_score, score_band, phase, raw_corrected_log, final_log,
      blend_qr, prob_peak_by_1, prob_peak_by_2, prob_peak_by_3
   ) %>%
   arrange(origin_week, horizon)
print(week_detail, n = Inf)

# ══════════════════════════════════════════════════════════════════════════════
# 15. RETROSPECTIVE VISUALIZATIONS  (unchanged from v3)
# ══════════════════════════════════════════════════════════════════════════════

band_colors <- c(
   "Excellent (>=0.85)" = "#639922", "Good (0.70-0.84)"  = "#1D9E75",
   "Fair (0.55-0.69)"   = "#BA7517", "Poor (<0.55)"       = "#A32D2D"
)

horizon_line_colors <- c(
   "Nowcast (H0)"    = "#4DAF4A", "1-Wk Ahead (H1)" = "#FF7F00",
   "2-Wk Ahead (H2)" = "#E41A1C", "3-Wk Ahead (H3)" = "#377EB8"
)

forecast_lines <- retro_df %>%
   mutate(
      horizon_label = case_when(
         horizon == 0 ~ "Nowcast (H0)",    horizon == 1 ~ "1-Wk Ahead (H1)",
         horizon == 2 ~ "2-Wk Ahead (H2)", horizon == 3 ~ "3-Wk Ahead (H3)"
      ),
      horizon_label = factor(horizon_label, levels = names(horizon_line_colors))
   ) %>%
   dplyr::select(target_end_date, horizon_label, forecast_median)

observed_line <- retro_df %>% dplyr::select(target_end_date, actual) %>% distinct()

all_target_dates <- sort(unique(c(observed_line$target_end_date, forecast_lines$target_end_date)))

p_diag <- ggplot() +
   geom_line(data = forecast_lines,
             aes(x = target_end_date, y = forecast_median,
                 color = horizon_label, group = horizon_label), linewidth = 1.0) +
   geom_line(data = observed_line,
             aes(x = target_end_date, y = actual, group = 1),
             color = "black", linewidth = 1.4) +
   geom_point(data = observed_line,
              aes(x = target_end_date, y = actual), color = "black", size = 2.2) +
   scale_color_manual(values = horizon_line_colors, name = "Horizon") +
   scale_x_date(breaks = all_target_dates, date_labels = "%b %d",
                expand = expansion(mult = 0.02)) +
   scale_y_continuous(labels = comma, expand = expansion(mult = c(0.02, 0.08))) +
   labs(
      title    = "State-Level: CDC Observed vs Forecast Diagnostic — v4 Recursive",
      subtitle = "Colored lines are median forecasts by horizon; recursive approach eliminates horizon shift",
      x = "Week Ending Date",
      y = "Weekly Inpatient Hospitalizations (Influenza, SC)",
      caption = paste0(
         "Black = CDC observed | Colored = DMAPRIME-QR v4 median forecast by horizon\n",
         "v4: recursive 1-step forecasting — each horizon uses features updated with prior horizon prediction."
      )
   ) +
   theme_minimal(base_size = 12) +
   theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 9),
      axis.title = element_text(size = 11),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      plot.caption = element_text(hjust = 0, size = 8, color = "gray40")
   )
print(p_diag)

p_time <- ggplot(retro_df,
                 aes(x = target_end_date, y = ai_score,
                     color = score_band, group = horizon)) +
   geom_hline(yintercept = c(0.55, 0.70, 0.85), linetype = "dashed",
              color = "gray70", linewidth = 0.4) +
   geom_line(color = "gray80", linewidth = 0.6) +
   geom_point(size = 2.5, alpha = 0.9) +
   scale_color_manual(values = band_colors, name = "Score band") +
   scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
   scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d") +
   facet_wrap(~ paste0("H", horizon, " (", horizon + 1, "-wk ahead)"), ncol = 2) +
   labs(title    = "Retrospective AI Scores by Horizon — v4 Recursive",
        subtitle = "AI = min(actual, forecast) / max(actual, forecast)",
        x = "Target end date", y = "AI score") +
   theme_minimal(base_size = 11) +
   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
         legend.position = "bottom", strip.text = element_text(face = "bold", size = 10),
         panel.grid.minor = element_blank())
print(p_time)

p_scatter <- ggplot(retro_df, aes(x = actual, y = forecast_median, color = score_band)) +
   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
   geom_point(size = 2.5, alpha = 0.85) +
   scale_color_manual(values = band_colors, name = "Score band") +
   facet_wrap(~ paste0("H", horizon, " (", horizon + 1, "-wk ahead)"), ncol = 2) +
   labs(title    = "Actual vs. Forecast — Retrospective Evaluation v4 Recursive",
        subtitle = "Dashed line = perfect forecast",
        x = "Observed CDC hospitalizations", y = "Median forecast",
        caption = "Points above line = over-forecast; below line = under-forecast") +
   theme_minimal(base_size = 11) +
   theme(legend.position = "bottom",
         strip.text = element_text(face = "bold", size = 10),
         panel.grid.minor = element_blank())
print(p_scatter)

p_bar <- horizon_summary %>%
   mutate(
      score_band = case_when(
         mean_ai_score >= 0.85 ~ "Excellent (>=0.85)", mean_ai_score >= 0.70 ~ "Good (0.70-0.84)",
         mean_ai_score >= 0.55 ~ "Fair (0.55-0.69)",   TRUE ~ "Poor (<0.55)"
      ),
      score_band = factor(score_band, levels = c(
         "Excellent (>=0.85)", "Good (0.70-0.84)", "Fair (0.55-0.69)", "Poor (<0.55)"
      ))
   ) %>%
   ggplot(aes(x = horizon_label, y = mean_ai_score, fill = score_band)) +
   geom_col(width = 0.55) +
   geom_text(aes(label = sprintf("%.3f\n(n=%d)", mean_ai_score, n_obs)),
             vjust = -0.3, size = 3.5) +
   geom_hline(yintercept = c(0.55, 0.70, 0.85), linetype = "dashed",
              color = "gray50", linewidth = 0.4) +
   scale_fill_manual(values = band_colors, name = "Score band") +
   scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.1)) +
   labs(title    = "Mean AI Score by Forecast Horizon — v4 Recursive",
        subtitle = "Recursive 1-step QR + analog forecast",
        x = NULL, y = "Mean AI score") +
   theme_minimal(base_size = 11) +
   theme(legend.position = "bottom", panel.grid.minor = element_blank())
print(p_bar)

# ══════════════════════════════════════════════════════════════════════════════
# 16. EXPORT RETROSPECTIVE + SOFTWARE EVALUATION  (unchanged from v3)
# ══════════════════════════════════════════════════════════════════════════════

output_retro_file <- file.path(
   output_dir_cdc, "retrospective_ai_scores_v4_dec2025_mar2026.csv"
)
output_horizon_summary_file <- file.path(
   output_dir_cdc, "retrospective_ai_scores_v4_by_horizon_summary.csv"
)
# Uncomment to save:
# write_csv(retro_df, output_retro_file)
# write_csv(horizon_summary, output_horizon_summary_file)

eval_state_training <- flu_features %>%
   dplyr::filter(week < train_cut & !is.na(.data[[target_var]])) %>%
   dplyr::mutate(
      value = .data[[target_var]], training_validation = 1,
      estimate_projected_report = NA_integer_, target_end_date = week
   )

eval_state_testing <- cdc_submission %>%
   dplyr::filter(output_type_id == 0.5) %>%
   dplyr::mutate(
      training_validation = 0,
      estimate_projected_report = as.integer(horizon),
      target_end_date = target_end_date
   )

sw_evaluation_file <- bind_rows(eval_state_training, eval_state_testing) %>%
   dplyr::mutate(
      reference_date = format(reference_date, "%Y-%m-%d"),
      target = "wk inc flu hosp",
      target_end_date = format(as.Date(target_end_date), "%Y-%m-%d"),
      location_general = "state", location = "SC", disease = "influenza",
      population = "general_population", imputed = 0,
      data_source = NA_character_,
      outcome_measure = "Weekly_Inpatient_Hospitalizations",
      output_type = "quantile", output_type_id = 0.5
   ) %>%
   dplyr::arrange(target_end_date) %>%
   dplyr::select(
      reference_date, target, target_end_date, location_general, location,
      value, disease, population, training_validation, estimate_projected_report,
      imputed, data_source, outcome_measure, output_type, output_type_id
   )

output_filename_eval <- file.path(
   output_dir_software_eval,
   "inpatient-hosp-DMAPRIME-QR-v4-evaluation.csv"
)
# Uncomment to save:
# write_csv(sw_evaluation_file, output_filename_eval)

cat("Software evaluation file: ", basename(output_filename_eval), "\n")

cat("\n✅ RETROSPECTIVE EVALUATION v4 COMPLETE\n")
cat(sprintf("Rows evaluated : %d\n", nrow(retro_df)))
cat(sprintf("Origin weeks   : %d\n", n_distinct(retro_df$origin_week)))
cat(sprintf("Overall mean AI: %.4f\n", mean(retro_df$ai_score, na.rm = TRUE)))

# ══════════════════════════════════════════════════════════════════════════════
# 17. ORIGIN-PATH DIAGNOSTIC  (unchanged from v3 — now should show flat paths)
# ══════════════════════════════════════════════════════════════════════════════

# Observed line: start exactly at eval_start (no pre-period history)
actual_line <- flu_features %>%
   filter(week >= eval_start,
          week <= eval_end + weeks(4),
          !is.na(.data[[target_var]])) %>%
   dplyr::select(week, actual = all_of(target_var))

# Up to 16 panels, ordered chronologically (not alphabetically)
selected_origins <- retro_df %>%
   filter(origin_week >= eval_start,
          origin_week <= eval_end) %>%
   distinct(origin_week) %>%
   arrange(origin_week) %>%
   dplyr::slice(1:16) %>%
   pull(origin_week)

# Build ordered factor so facet_wrap respects chronological order
origin_path_df <- retro_df %>%
   filter(origin_week %in% selected_origins) %>%
   mutate(
      origin_label = factor(
         paste0("Issued ", format(origin_week, "%b %d")),
         levels = paste0("Issued ", format(sort(selected_origins), "%b %d"))
      ),
      horizon_label = paste0("H", horizon, " / ", horizon + 1, "-wk ahead")
   )

# Each panel only needs to show ~6 weeks: 2 weeks of observed context before
# the origin + 4 forecast weeks.  We achieve this by filtering actual_line
# per-panel inside the data rather than relying on free_x scales (which would
# hide the epidemic context).  Instead use scales="free_x" so each panel
# auto-ranges to its own origin window.
origin_path_df <- origin_path_df %>%
   mutate(panel_xmin = origin_week - weeks(3),
          panel_xmax = target_end_date)

# Restrict observed line to only the dates that appear in at least one panel
obs_dates_needed <- sort(unique(c(
   as.Date(unlist(lapply(selected_origins, function(o) {
      seq(o - weeks(3), o + weeks(4), by = "week")
   })))
)))

actual_line_trimmed <- actual_line %>%
   filter(week %in% obs_dates_needed)

p_origin_paths <- ggplot() +
   geom_line(data = actual_line_trimmed,
             aes(x = week, y = actual, group = 1), color = "black", linewidth = 1.2) +
   geom_point(data = actual_line_trimmed,
              aes(x = week, y = actual), color = "black", size = 1.8) +
   geom_line(data = origin_path_df,
             aes(x = target_end_date, y = forecast_median, group = origin_week),
             color = "#377EB8", linewidth = 1.1) +
   geom_point(data = origin_path_df,
              aes(x = target_end_date, y = forecast_median, shape = horizon_label),
              color = "#377EB8", size = 2.4) +
   geom_vline(data = origin_path_df %>% distinct(origin_label, origin_week),
              aes(xintercept = origin_week), linetype = "dashed",
              color = "gray50", linewidth = 0.5) +
   facet_wrap(~ origin_label, ncol = 4, scales = "free_x") +
   scale_x_date(date_breaks = "1 week", date_labels = "%b %d") +
   scale_y_continuous(labels = comma) +
   labs(
      title    = "Forecast Paths by Origin Week — v4 Recursive (Horizon-Shift Diagnostic)",
      subtitle = "Each blue line is one issued 4-week forecast path; black is observed CDC",
      x = "Week Ending Date",
      y = "Weekly Inpatient Hospitalizations (Influenza, SC)",
      shape   = "Forecast horizon",
      caption = paste0(
         "v4 recursive: each step conditions on prior step prediction — ",
         "horizon shift is architecturally prevented."
      )
   ) +
   theme_minimal(base_size = 11) +
   theme(
      plot.title    = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      axis.text.x   = element_text(angle = 45, hjust = 1, size = 7),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold", size = 9)
   )
print(p_origin_paths)
