library(dplyr)
library(lubridate)
library(ggplot2)
library(purrr)
library(readr)
library(quantreg)
library(tidyr)
library(slider)

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

# ---------------- 1. LOAD DATA ----------------
RFA_data_Influenza       <- read_csv('/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Data/SC Health Records/RFA/Years 2017-2025/Weekly data/RFA_weekly_influenza_region_incident.csv')


CDC <- read_csv("/Users/tanvirahammed/Downloads/target-hospital-admissions-19.csv") %>%
   filter(location == "45") %>%
   dplyr::select(
      week = date,
      cdc_Total.Influenza.Admissions = value
   ) %>%
   mutate(week = as.Date(week))

MUSC_Weekly_Influenza <- read_csv(
   "/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Data/MUSC/Infectious Disease EHR/Weekly Data/Latest Weekly Data/MUSC_Weekly_Influenza_State_dx_cond_lab_Incident.csv"
)

Prisma_Weekly_Influenza <- read_csv(
   "/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Data/Prisma Health/Infectious Disease EHR/Weekly Data/Latest Weekly Data/Prisma_Health_Weekly_Influenza_State_dx_cond_lab_Incident.csv"
)

prep_source <- function(df, date_col, prefix) {
   df %>%
      rename(week = !!rlang::sym(date_col)) %>%
      mutate(week = as.Date(week)) %>%
      rename_with(~ paste0(prefix, .x), .cols = -week)
}

prisma <- prep_source(Prisma_Weekly_Influenza, "Week", "prisma_")
musc <- prep_source(MUSC_Weekly_Influenza, "Week", "musc_")

flu <- CDC %>%
   full_join(prisma, by = "week") %>%
   full_join(musc, by = "week") %>%
   arrange(week)

target_var <- "cdc_Total.Influenza.Admissions"

safe_log <- function(x) log(pmax(x, 0) + 1)

weeks_since_peak_fun <- function(x) {
   if (all(is.na(x))) return(NA_integer_)
   length(x) - which.max(replace_na(x, -Inf))
}

weighted_quantile_safe <- function(x, w, probs) {
   ok <- !is.na(x) & !is.na(w)
   x <- x[ok]
   w <- w[ok]
   
   if (length(x) == 0 || sum(w) <= 0) {
      return(rep(NA_real_, length(probs)))
   }
   
   ord <- order(x)
   x <- x[ord]
   w <- w[ord] / sum(w[ord])
   cw <- cumsum(w)
   
   sapply(probs, function(p) x[which(cw >= p)[1]])
}

weighted_median_safe <- function(x, w) {
   weighted_quantile_safe(x, w, probs = 0.5)[1]
}

standardize_using <- function(df, ref_df, vars) {
   out <- df
   
   for (v in vars) {
      m <- mean(ref_df[[v]], na.rm = TRUE)
      s <- sd(ref_df[[v]], na.rm = TRUE)
      if (is.na(s) || s == 0) s <- 1
      out[[v]] <- (out[[v]] - m) / s
   }
   
   out
}

# ---------------- 2. FEATURE ENGINEERING ----------------

flu_features <- flu %>%
   arrange(week) %>%
   mutate(
      log_target = safe_log(.data[[target_var]]),
      
      log_tests = safe_log(prisma_Weekly_Tests),
      log_positive_tests = safe_log(prisma_Weekly_Positive_Tests),
      log_hosp = safe_log(prisma_Weekly_Inpatient_Hospitalizations),
      
      lag1 = lag(log_target, 1),
      lag2 = lag(log_target, 2),
      lag3 = lag(log_target, 3),
      lag4 = lag(log_target, 4),
      lag5 = lag(log_target, 5),
      lag8 = lag(log_target, 8),
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
      
      mean_4wk = slide_dbl(log_target, mean, .before = 3, .complete = FALSE, na.rm = TRUE),
      mean_8wk = slide_dbl(log_target, mean, .before = 7, .complete = FALSE, na.rm = TRUE),
      mean_12wk = slide_dbl(log_target, mean, .before = 11, .complete = FALSE, na.rm = TRUE),
      
      sd_4wk = slide_dbl(log_target, sd, .before = 3, .complete = FALSE, na.rm = TRUE),
      sd_8wk = slide_dbl(log_target, sd, .before = 7, .complete = FALSE, na.rm = TRUE),
      sd_12wk = slide_dbl(log_target, sd, .before = 11, .complete = FALSE, na.rm = TRUE),
      
      max_12wk = slide_dbl(log_target, max, .before = 11, .complete = FALSE, na.rm = TRUE),
      min_12wk = slide_dbl(log_target, min, .before = 11, .complete = FALSE, na.rm = TRUE),
      
      position_in_12wk_range = (log_target - min_12wk) / pmax(max_12wk - min_12wk, 1e-6),
      
      weeks_since_12wk_peak = slide_int(
         log_target,
         weeks_since_peak_fun,
         .before = 11,
         .complete = FALSE
      ),
      
      peak_12wk_log = max_12wk,
      drop_from_12wk_peak = log_target - max_12wk,
      drop_from_12wk_peak_pct = (exp(log_target) - exp(max_12wk)) / pmax(exp(max_12wk), 1e-6),
      post_peak_decline_speed = drop_from_12wk_peak / pmax(weeks_since_12wk_peak, 1),
      
      hosp_slope_1wk = log_hosp - lag(log_hosp, 1),
      hosp_slope_2wk = log_hosp - lag(log_hosp, 2),
      pos_slope_1wk = log_positive_tests - lag(log_positive_tests, 1),
      pos_slope_2wk = log_positive_tests - lag(log_positive_tests, 2),
      
      ehr_decline_confirmed = as.integer(
         hosp_slope_1wk < 0 |
            pos_slope_1wk < 0 |
            hosp_slope_2wk < 0
      ),
      
      week_of_year = isoweek(week),
      sin52 = sin(2 * pi * week_of_year / 52),
      cos52 = cos(2 * pi * week_of_year / 52),
      
      month = month(week),
      is_christmas = if_else(month == 12 & day(week) >= 20, 1, 0),
      is_newyear = if_else(month == 1 & day(week) <= 7, 1, 0)
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
      )
   ) %>%
   dplyr::select(-slope_1wk_q25)

for (h in horizons) {
   flu_features[[paste0("future_log_h", h)]] <- lead(flu_features$log_target, h + 1)
   flu_features[[paste0("future_value_h", h)]] <- lead(flu_features[[target_var]], h + 1)
   flu_features[[paste0("future_change_h", h)]] <-
      flu_features[[paste0("future_log_h", h)]] - flu_features$log_target
}

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
   "position_in_12wk_range",
   "weeks_since_12wk_peak",
   "peak_12wk_log",
   "drop_from_12wk_peak",
   "drop_from_12wk_peak_pct",
   "post_peak_decline_speed",
   "is_post_peak",
   "is_fast_decline",
   "sin52", "cos52",
   "is_christmas", "is_newyear"
)

state_vars <- c(
   "log_target", "lag1",
   "slope_1wk", "slope_2wk", "slope_3wk", "slope_4wk",
   "prev_slope_1wk", "prev_slope_2wk",
   "accel_2wk", "accel_3wk",
   "decline_persistence_2wk", "decline_persistence_3wk",
   "mean_4wk", "mean_8wk", "mean_12wk",
   "sd_4wk", "sd_8wk", "sd_12wk",
   "position_in_12wk_range",
   "weeks_since_12wk_peak",
   "peak_12wk_log",
   "drop_from_12wk_peak",
   "drop_from_12wk_peak_pct",
   "post_peak_decline_speed",
   "is_post_peak",
   "is_fast_decline",
   "log_hosp",
   "hosp_slope_1wk",
   "hosp_slope_2wk",
   "log_positive_tests",
   "pos_slope_1wk",
   "pos_slope_2wk",
   "ehr_decline_confirmed",
   "sin52", "cos52"
)

# ---------------- 3. PHASE FUNCTION ----------------

phase_label <- function(row) {
   
   slope <- row$slope_1wk
   accel <- row$accel_2wk
   vol <- row$sd_4wk
   
   slope_cut <- sd(flu_features$slope_1wk, na.rm = TRUE) * 0.25
   accel_cut <- sd(flu_features$accel_2wk, na.rm = TRUE) * 0.25
   vol_cut <- quantile(flu_features$sd_4wk, 0.70, na.rm = TRUE)
   
   case_when(
      row$is_fast_decline == 1 ~ "FAST_POST_PEAK_DECLINE",
      row$is_post_peak == 1 & row$decline_persistence_2wk == 1 ~ "FAST_POST_PEAK_DECLINE",
      row$is_post_peak == 1 & slope < 0 ~ "POST_PEAK_DECLINE",
      slope < -slope_cut & accel < -accel_cut ~ "ACCELERATING_DECLINE",
      slope < -slope_cut ~ "SLOWING_DECLINE",
      slope > slope_cut & accel > accel_cut ~ "ACCELERATING_INCREASE",
      slope > slope_cut ~ "SLOWING_INCREASE",
      abs(slope) <= slope_cut & vol > vol_cut ~ "VOLATILE_PLATEAU",
      abs(slope) <= slope_cut ~ "STABLE_PLATEAU",
      TRUE ~ "UNCERTAIN"
   )
}

flu_features <- flu_features %>%
   rowwise() %>%
   mutate(phase = phase_label(cur_data())) %>%
   ungroup()

decline_qr_floor <- function(row, phase) {
   case_when(
      phase == "FAST_POST_PEAK_DECLINE" ~ 0.65,
      phase == "POST_PEAK_DECLINE" ~ 0.55,
      phase == "ACCELERATING_DECLINE" ~ 0.55,
      phase == "SLOWING_DECLINE" ~ 0.45,
      TRUE ~ 0.00
   )
}

interval_scale_by_phase <- function(phase) {
   case_when(
      phase == "FAST_POST_PEAK_DECLINE" ~ 0.60,
      phase == "POST_PEAK_DECLINE" ~ 0.70,
      phase == "ACCELERATING_DECLINE" ~ 0.75,
      phase == "VOLATILE_PLATEAU" ~ 0.85,
      TRUE ~ 1.00
   )
}

# ---------------- 4. MODEL HELPERS ----------------

prepare_qr_design <- function(model_df, target_h, candidate_predictors) {
   
   model_df <- model_df %>%
      drop_na(all_of(c(target_h, candidate_predictors)))
   
   if (nrow(model_df) < 30) stop("Too few complete rows after removing missing values.")
   
   usable_predictors <- candidate_predictors[
      sapply(candidate_predictors, function(v) {
         x <- model_df[[v]]
         sd(x, na.rm = TRUE) > 1e-8 &&
            length(unique(x[!is.na(x)])) > 1
      })
   ]
   
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

fit_qr_model <- function(train_data, h, tau_value = 0.5) {
   
   target_h <- paste0("future_log_h", h)
   
   prepared <- prepare_qr_design(
      model_df = train_data,
      target_h = target_h,
      candidate_predictors = predictors
   )
   
   model <- rq(
      as.formula(paste(target_h, "~", paste(prepared$predictors, collapse = " + "))),
      tau = tau_value,
      data = prepared$data,
      method = "fn"
   )
   
   attr(model, "used_predictors") <- prepared$predictors
   model
}

predict_qr_safe <- function(model, newdata) {
   as.numeric(predict(model, newdata = newdata))
}

get_analogs <- function(origin_row, analog_pool, vars, k, bandwidth, h = NULL) {
   
   pool_complete <- analog_pool %>% drop_na(all_of(vars))
   if (nrow(pool_complete) < 3) return(pool_complete[0, ])
   
   scaled_pool <- standardize_using(pool_complete, pool_complete, vars)
   scaled_origin <- standardize_using(origin_row, pool_complete, vars)
   
   var_weights <- rep(1, length(vars))
   names(var_weights) <- vars
   
   high_weight_vars <- c(
      "slope_1wk", "slope_2wk", "slope_3wk", "slope_4wk",
      "prev_slope_1wk", "prev_slope_2wk",
      "accel_2wk", "accel_3wk",
      "decline_persistence_2wk", "decline_persistence_3wk",
      "hosp_slope_1wk", "hosp_slope_2wk",
      "pos_slope_1wk", "pos_slope_2wk",
      "ehr_decline_confirmed",
      "drop_from_12wk_peak",
      "drop_from_12wk_peak_pct",
      "post_peak_decline_speed",
      "is_post_peak",
      "is_fast_decline"
   )
   
   medium_weight_vars <- c(
      "position_in_12wk_range",
      "weeks_since_12wk_peak",
      "peak_12wk_log"
   )
   
   var_weights[names(var_weights) %in% high_weight_vars] <- 3.0
   var_weights[names(var_weights) %in% medium_weight_vars] <- 1.5
   
   pool_matrix <- as.matrix(scaled_pool[, vars])
   origin_vector <- as.numeric(scaled_origin[1, vars])
   
   diff_mat <- sweep(pool_matrix, 2, origin_vector, "-")
   weighted_diff <- sweep(diff_mat^2, 2, var_weights, "*")
   
   pool_complete$distance <- sqrt(rowSums(weighted_diff))
   
   selected <- pool_complete %>%
      arrange(distance) %>%
      slice(1:min(k, n())) %>%
      mutate(
         analog_weight = exp(-distance / bandwidth),
         analog_weight = analog_weight / sum(analog_weight, na.rm = TRUE)
      )
   
   if (!is.null(h)) {
      change_col <- paste0("future_change_h", h)
      
      decline_now <- origin_row$is_fast_decline == 1 ||
         origin_row$is_post_peak == 1 ||
         origin_row$decline_persistence_2wk == 1
      
      if (decline_now && change_col %in% names(selected)) {
         selected <- selected %>%
            mutate(
               rebound_penalty = if_else(.data[[change_col]] > 0, 0.20, 1.00),
               weak_decline_penalty = if_else(.data[[change_col]] > -0.05, 0.60, 1.00),
               analog_weight = analog_weight * rebound_penalty * weak_decline_penalty,
               analog_weight = analog_weight / sum(analog_weight, na.rm = TRUE)
            )
      }
   }
   
   selected
}

# ---------------- 5. ROLLING VALIDATION ----------------

validation_origins <- flu_features %>%
   filter(
      week >= train_cut,
      week < valid_end,
      !is.na(log_target)
   ) %>%
   drop_na(all_of(c(predictors, state_vars))) %>%
   pull(week)

tune_results <- list()

for (h in horizons) {
   
   cat("\nTuning horizon", h, "...\n")
   
   pred_store <- list()
   
   for (origin_week in validation_origins) {
      
      train_origin <- flu_features %>%
         filter(week < origin_week) %>%
         drop_na(all_of(c(paste0("future_log_h", h), predictors)))
      
      origin_row <- flu_features %>%
         filter(week == origin_week) %>%
         drop_na(all_of(c(predictors, state_vars)))
      
      y_true <- flu_features %>%
         filter(week == origin_week) %>%
         pull(paste0("future_log_h", h))
      
      if (nrow(train_origin) < 80 || nrow(origin_row) != 1 || is.na(y_true)) next
      
      qr_model_h <- tryCatch(
         fit_qr_model(train_origin, h, tau_value = 0.5),
         error = function(e) NULL
      )
      
      if (is.null(qr_model_h)) next
      
      qr_pred <- tryCatch(
         predict_qr_safe(qr_model_h, origin_row),
         error = function(e) NA_real_
      )
      
      if (is.na(qr_pred)) next
      
      analog_pool <- flu_features %>%
         filter(week < origin_week) %>%
         drop_na(all_of(c(state_vars, paste0("future_change_h", h))))
      
      grid_origin <- expand_grid(
         k = candidate_k,
         bandwidth = candidate_bandwidth
      ) %>%
         mutate(
            origin_week = origin_week,
            y_true = y_true,
            qr_pred = qr_pred,
            analog_pred = NA_real_
         )
      
      for (g in seq_len(nrow(grid_origin))) {
         
         analogs_g <- get_analogs(
            origin_row = origin_row,
            analog_pool = analog_pool,
            vars = state_vars,
            k = grid_origin$k[g],
            bandwidth = grid_origin$bandwidth[g],
            h = h
         )
         
         if (nrow(analogs_g) < 5) next
         
         change_col <- paste0("future_change_h", h)
         
         decline_now <- origin_row$is_fast_decline == 1 ||
            origin_row$is_post_peak == 1 ||
            origin_row$decline_persistence_2wk == 1
         
         analog_change <- if (decline_now) {
            weighted_quantile_safe(
               analogs_g[[change_col]],
               analogs_g$analog_weight,
               probs = 0.35
            )[1]
         } else {
            weighted_median_safe(
               analogs_g[[change_col]],
               analogs_g$analog_weight
            )
         }
         
         grid_origin$analog_pred[g] <- origin_row$log_target + analog_change
      }
      
      pred_store[[as.character(origin_week)]] <- grid_origin
   }
   
   if (length(pred_store) == 0) {
      tune_results[[paste0("h", h)]] <- tibble(
         horizon = h,
         k = 20,
         bandwidth = 1.25,
         blend_qr = 0.50,
         mae_log = NA_real_,
         n_valid = 0
      )
      next
   }
   
   pred_df <- bind_rows(pred_store) %>%
      filter(!is.na(analog_pred), !is.na(qr_pred), !is.na(y_true))
   
   if (nrow(pred_df) == 0) {
      tune_results[[paste0("h", h)]] <- tibble(
         horizon = h,
         k = 20,
         bandwidth = 1.25,
         blend_qr = 0.50,
         mae_log = NA_real_,
         n_valid = 0
      )
      next
   }
   
   score_df <- pred_df %>%
      crossing(blend_qr = candidate_blend_qr) %>%
      mutate(
         pred = blend_qr * qr_pred + (1 - blend_qr) * analog_pred,
         abs_error = abs(y_true - pred)
      ) %>%
      group_by(k, bandwidth, blend_qr) %>%
      summarise(
         mae_log = mean(abs_error, na.rm = TRUE),
         n_valid = n(),
         .groups = "drop"
      ) %>%
      arrange(mae_log)
   
   best <- score_df %>%
      slice(1) %>%
      mutate(horizon = h)
   
   tune_results[[paste0("h", h)]] <- best
}

tuning_table <- bind_rows(tune_results) %>%
   dplyr::select(horizon, k, bandwidth, blend_qr, mae_log, n_valid)

print(tuning_table)

# ---------------- 6. FINAL QR MODELS ----------------

final_train_df <- flu_features %>%
   filter(week < valid_end) %>%
   drop_na(all_of(c(predictors, state_vars)))

qr_models <- list()
residual_store <- list()

for (h in horizons) {
   
   target_h <- paste0("future_log_h", h)
   
   model_df <- final_train_df %>%
      drop_na(all_of(c(target_h, predictors)))
   
   qr_models[[paste0("h", h)]] <- fit_qr_model(model_df, h, tau_value = 0.5)
   
   model_df$pred <- predict(qr_models[[paste0("h", h)]], newdata = model_df)
   model_df$residual <- model_df[[target_h]] - model_df$pred
   model_df$horizon <- h
   
   residual_store[[paste0("h", h)]] <- model_df %>%
      dplyr::select(week, horizon, residual, phase)
}

residual_df <- bind_rows(residual_store)

# ---------------- 7. FINAL FORECAST ----------------

latest_observed_week <- max(flu_features$week[!is.na(flu_features[[target_var]])], na.rm = TRUE)
reference_date <- latest_observed_week + weeks(1)

current_row <- flu_features %>%
   filter(week == latest_observed_week) %>%
   drop_na(all_of(c(predictors, state_vars)))

if (nrow(current_row) != 1) {
   stop("Current row is incomplete. Check missing predictors for latest observed week.")
}

observed_value <- current_row[[target_var]]
phase <- current_row$phase

forecast_rows <- list()
analog_diagnostics <- list()

for (h in horizons) {
   
   tune_h <- tuning_table %>%
      filter(horizon == h) %>%
      slice(1)
   
   k_h <- tune_h$k
   bandwidth_h <- tune_h$bandwidth
   blend_qr_raw <- tune_h$blend_qr
   
   qr_floor <- decline_qr_floor(current_row, phase)
   blend_qr_h <- max(blend_qr_raw, qr_floor)
   blend_analog_h <- 1 - blend_qr_h
   
   qr_median_log <- predict_qr_safe(qr_models[[paste0("h", h)]], current_row)
   
   analog_pool <- flu_features %>%
      filter(week < latest_observed_week) %>%
      drop_na(all_of(c(state_vars, paste0("future_change_h", h))))
   
   selected_analogs <- get_analogs(
      origin_row = current_row,
      analog_pool = analog_pool,
      vars = state_vars,
      k = k_h,
      bandwidth = bandwidth_h,
      h = h
   )
   
   change_col <- paste0("future_change_h", h)
   
   decline_now <- current_row$is_fast_decline == 1 ||
      current_row$is_post_peak == 1 ||
      current_row$decline_persistence_2wk == 1
   
   analog_change <- if (decline_now) {
      weighted_quantile_safe(
         selected_analogs[[change_col]],
         selected_analogs$analog_weight,
         probs = 0.35
      )[1]
   } else {
      weighted_median_safe(
         selected_analogs[[change_col]],
         selected_analogs$analog_weight
      )
   }
   
   analog_target_log <- current_row$log_target + analog_change
   
   corrected_median_log <- blend_qr_h * qr_median_log +
      blend_analog_h * analog_target_log
   
   corrected_median_log <- pmax(corrected_median_log, 0)
   
   analog_error <- selected_analogs[[change_col]] - analog_change
   
   analog_error_q <- weighted_quantile_safe(
      analog_error,
      selected_analogs$analog_weight,
      quantiles
   )
   
   residual_phase_h <- residual_df %>%
      filter(horizon == h, phase == !!phase) %>%
      pull(residual)
   
   if (length(residual_phase_h) < 25) {
      residual_phase_h <- residual_df %>%
         filter(horizon == h) %>%
         pull(residual)
   }
   
   residual_error_q <- quantile(
      residual_phase_h,
      probs = quantiles,
      na.rm = TRUE,
      names = FALSE
   )
   
   combined_error_q <- blend_analog_h * analog_error_q +
      blend_qr_h * residual_error_q
   
   combined_error_q <- combined_error_q * interval_scale_by_phase(phase)
   
   lower_cap <- quantile(combined_error_q, 0.025, na.rm = TRUE)
   upper_cap <- quantile(combined_error_q, 0.975, na.rm = TRUE)
   
   combined_error_q <- pmin(pmax(combined_error_q, lower_cap), upper_cap)
   
   log_quantiles <- corrected_median_log + combined_error_q
   values <- pmax(exp(log_quantiles) - 1, 0)
   
   forecast_rows[[paste0("h", h)]] <- data.frame(
      week = reference_date + weeks(h),
      horizon = h,
      quantile = quantiles,
      predicted_value = values,
      qr_median_log = qr_median_log,
      analog_target_log = analog_target_log,
      corrected_median_log = corrected_median_log,
      blend_qr_raw = blend_qr_raw,
      blend_qr = blend_qr_h,
      blend_analog = blend_analog_h,
      phase = phase,
      k = k_h,
      bandwidth = bandwidth_h
   )
   
   analog_diagnostics[[paste0("h", h)]] <- selected_analogs %>%
      mutate(
         horizon = h,
         future_change = .data[[change_col]]
      ) %>%
      dplyr::select(
         horizon,
         week,
         distance,
         analog_weight,
         future_change,
         all_of(target_var),
         is_post_peak,
         is_fast_decline,
         decline_persistence_2wk,
         decline_persistence_3wk,
         slope_1wk,
         slope_2wk,
         slope_3wk,
         drop_from_12wk_peak,
         post_peak_decline_speed
      )
}

all_quantiles <- bind_rows(forecast_rows) %>%
   group_by(week) %>%
   arrange(quantile, .by_group = TRUE) %>%
   mutate(predicted_value = isoreg(predicted_value)$yf) %>%
   ungroup()

analog_diagnostics_df <- bind_rows(analog_diagnostics)

cdc_submission <- all_quantiles %>%
   mutate(
      reference_date = reference_date,
      target = "wk inc flu hosp",
      location = "45",
      target_end_date = week,
      output_type = "quantile",
      output_type_id = quantile,
      value = round(predicted_value, 0)
   ) %>%
   dplyr::select(
      reference_date,
      target,
      horizon,
      location,
      target_end_date,
      output_type,
      output_type_id,
      value
   ) %>%
   arrange(horizon, output_type_id)

cat("\n═══ FINAL FORECAST DIAGNOSTICS ═══\n")
cat(sprintf("Latest observed week: %s\n", latest_observed_week))
cat(sprintf("Latest observed value: %s\n", observed_value))
cat(sprintf("Detected phase: %s\n", phase))

cat("\nLatest row phase variables:\n")
print(
   current_row %>%
      dplyr::select(
         week,
         all_of(target_var),
         log_target,
         lag1,
         slope_1wk,
         slope_2wk,
         prev_slope_1wk,
         accel_2wk,
         weeks_since_12wk_peak,
         drop_from_12wk_peak,
         is_post_peak,
         is_fast_decline,
         decline_persistence_2wk,
         phase
      )
)

cat("\nForecast medians:\n")
print(
   cdc_submission %>%
      filter(output_type_id == 0.5) %>%
      dplyr::select(horizon, target_end_date, median_forecast = value)
)

cat("\nBlend used:\n")
print(
   all_quantiles %>%
      filter(quantile == 0.5) %>%
      dplyr::select(
         horizon,
         week,
         phase,
         blend_qr_raw,
         blend_qr,
         blend_analog,
         k,
         bandwidth
      )
)

cat("\nAnalog diagnostics:\n")
print(
   analog_diagnostics_df %>%
      group_by(horizon) %>%
      summarise(
         n_analogs = n(),
         median_future_change = median(future_change, na.rm = TRUE),
         mean_future_change = mean(future_change, na.rm = TRUE),
         mean_distance = mean(distance, na.rm = TRUE),
         max_weight = max(analog_weight, na.rm = TRUE),
         share_fast_decline_analogs = mean(is_fast_decline, na.rm = TRUE),
         share_post_peak_analogs = mean(is_post_peak, na.rm = TRUE),
         share_decline_persistence_2wk = mean(decline_persistence_2wk, na.rm = TRUE),
         .groups = "drop"
      )
)

# ---------------- 8. VISUALIZATION ----------------
# Important: removed misleading historical nowcast line.
# h0 model predicts next week, not current week.

plot_intervals <- cdc_submission %>%
   filter(output_type_id %in% c(0.025, 0.25, 0.5, 0.75, 0.975)) %>%
   pivot_wider(names_from = output_type_id, values_from = value) %>%
   rename(
      median = `0.5`,
      lower_95 = `0.025`,
      upper_95 = `0.975`,
      lower_50 = `0.25`,
      upper_50 = `0.75`
   )

hist_plot_df <- flu_features %>%
   filter(week >= as.Date("2025-10-01") & week <= latest_observed_week) %>%
   dplyr::select(
      week,
      Observed = all_of(target_var)
   ) %>%
   pivot_longer(
      cols = c(Observed),
      names_to = "Series",
      values_to = "Admissions"
   )

all_actual_dates <- sort(unique(c(hist_plot_df$week, plot_intervals$target_end_date)))

p1 <- ggplot() +
   geom_ribbon(
      data = plot_intervals,
      aes(x = target_end_date, ymin = lower_95, ymax = upper_95),
      fill = "#457B9D",
      alpha = 0.15
   ) +
   geom_ribbon(
      data = plot_intervals,
      aes(x = target_end_date, ymin = lower_50, ymax = upper_50),
      fill = "#457B9D",
      alpha = 0.30
   ) +
   geom_line(
      data = hist_plot_df,
      aes(x = week, y = Admissions, color = Series),
      linewidth = 1
   ) +
   geom_point(
      data = hist_plot_df,
      aes(x = week, y = Admissions, color = Series),
      size = 2
   ) +
   geom_line(
      data = plot_intervals,
      aes(x = target_end_date, y = median, color = "Forecast"),
      linetype = "dashed",
      linewidth = 1.2
   ) +
   geom_point(
      data = plot_intervals,
      aes(x = target_end_date, y = median, color = "Forecast"),
      size = 3
   ) +
   geom_vline(
      xintercept = latest_observed_week,
      color = "black",
      linewidth = 0.7
   ) +
   scale_x_date(
      breaks = all_actual_dates,
      date_labels = "%b %d"
   ) +
   scale_color_manual(
      values = c(
         "Observed" = "black",
         "Forecast" = "#457B9D"
      )
   ) +
   labs(
      title = paste0("Corrected Decline-Sensitive Direct QR + Analog Forecast [", phase, "]"),
      subtitle = paste0(
         "Latest: ", format(latest_observed_week, "%Y-%m-%d"),
         " (", observed_value, " admissions)"
      ),
      y = "Weekly Inpatient Hospitalizations (Influenza, South Carolina)",
      x = "Week Ending Date",
      caption = "Model: direct horizon-specific QR + decline-aware analog correction + phase-conditioned intervals"
   ) +
   theme_minimal() +
   theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 9),
      axis.title = element_text(size = 11),
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.grid.minor = element_blank(),
      plot.caption = element_text(hjust = 0, size = 8, color = "gray40")
   )

print(p1)

# ---------------- 9. EXPORT CDC FILE ----------------

output_filename_cdc <- paste0(
   "/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Tanvir/FluSight/Tanvir_QR/State/CDC_submission/",
   format(reference_date, "%Y-%m-%d"),
   "-DMAPRIME-QR.csv"
)

write_csv(cdc_submission, output_filename_cdc)

cat("\n✅ PIPELINE COMPLETE\n")
cat("CDC file: ", basename(output_filename_cdc), "\n")














# ══════════════════════════════════════════════════════════════════════════════
# 9. GENERATE SW IMPLEMENTATION FILE
# ══════════════════════════════════════════════════════════════════════════════

# A. Prepare Historical RFA Rows (State Level)
impl_state_rfa <- RFA_data_Influenza %>%
   dplyr::select(week_end, weekly_IP) %>%
   dplyr::filter(week_end >= as.Date("2020-03-07")) %>%
   dplyr::mutate(
      value = weekly_IP,
      data_source = "RFA",
      estimate_projected_report = 2, # Set to 2 for historical 
      target_end_date = week_end
   )

# B. Prepare Historical Healthsystem (HS) Rows (State Level)
# Aggregating MUSC and Prisma to get a single State total
impl_state_hs <- full_join(
   MUSC_Weekly_Influenza %>% dplyr::select(Week, MUSC_Hosp = Weekly_Inpatient_Hospitalizations),
   Prisma_Weekly_Influenza %>% dplyr::select(Week, Prisma_Hosp = Weekly_Inpatient_Hospitalizations),
   by = "Week"
) %>%
   dplyr::mutate(
      value = rowSums(across(c(MUSC_Hosp, Prisma_Hosp)), na.rm = TRUE),
      data_source = "HS",
      estimate_projected_report = 2, # Set to 2 for historical 
      target_end_date = Week
   ) %>%
   dplyr::filter(target_end_date >= as.Date("2020-03-07"))

# C. Prepare Future Forecast Rows (Quantile 0.5 only for implementation)
impl_state_forecast <- cdc_submission %>%
   #dplyr::filter(output_type_id == 0.5) %>%
   dplyr::mutate(
      data_source = NA_character_,
      estimate_projected_report = 1, # Set to 1 for forecasts 
      target_end_date = target_end_date
   )

# D. Combine and Apply Fixed SW Formatting
sw_implementation_file <- bind_rows(impl_state_rfa, impl_state_hs, impl_state_forecast) %>%
   dplyr::mutate(
      # reference_date: Saturday of current week [cite: 6]
      reference_date   = format(latest_observed_week + weeks(1), "%Y-%m-%d"), 
      
      # target: Fixed character string [cite: 11]
      target           = "wk inc flu hosp",           
      
      # target_end_date: ISO format [cite: 16]
      target_end_date  = format(target_end_date, "%Y-%m-%d"),
      
      # location_general: Fixed for state level [cite: 19]
      location_general = "state",                    
      
      # location: "SC" for state model
      location         = "SC",      
      
      # disease: influenza [cite: 35]
      disease          = "influenza",                 
      
      # population: general_population [cite: 38]
      population       = "general_population",        
      
      # training_validation: 0 for implementation 
      training_validation = 0,                        
      
      # imputed: Numeric 0 [cite: 58]
      imputed          = 0,                           
      
      # outcome_measure: Fixed label [cite: 70]
      outcome_measure  = "Weekly_Inpatient_Hospitalizations", 
      
      # output_type: quantile [cite: 74]
      output_type      = "quantile",                  
      
      # output_type_id: 0.5 for point estimates [cite: 82]
      #output_type_id   = 0.5                          
   ) %>%
   dplyr::select(
      reference_date, target, target_end_date, location_general, location, 
      value, disease, population, training_validation, 
      estimate_projected_report, imputed, data_source, 
      outcome_measure, output_type, output_type_id
   )
# Software Submission
# output_filename_software <- paste0(
#    "/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Tanvir/FluSight/Tanvir_QR/State/software.submission/",
#    format(reference_date, "%Y-%m-%d"),
#    "-Ahammed-Tanvir-17-state-flu-implementation.csv"
# )




output_filename_software <- paste0(
   "/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Forecasting Resources/Forecast-Drop-Off/Software/Implementation/inpatient-hosp-DMAPRIME-QR-implementation.csv"
)





# Optional: Export files (uncomment to save)
write_csv(sw_implementation_file, output_filename_software)


cat("\n✅ PIPELINE COMPLETE | CDC File: ", basename(output_filename_cdc), "\n")







# ══════════════════════════════════════════════════════════════════════════════
# 10. GENERATE STATE-LEVEL EVALUATION FILE
# ══════════════════════════════════════════════════════════════════════════════

# A. Prepare Training Data (Historical baseline from flu_features)
eval_state_training <- flu_features %>%
   # Filter to data before your training cut-off [cite: 43, 47]
   dplyr::filter(week < train_cut & !is.na(!!rlang::sym(target_var))) %>%
   dplyr::mutate(
      value = !!rlang::sym(target_var),
      training_validation = 1,           # Training: 1 
      estimate_projected_report = NA,     # NA for training baseline
      target_end_date = week
   )

# B. Prepare Testing Data (Recent Nowcasts + Future Forecasts)
eval_state_testing <- cdc_submission %>%
   # Use the median (0.5) for the evaluation file [cite: 80, 82]
   dplyr::filter(output_type_id == 0.5) %>%
   dplyr::mutate(
      training_validation = 0,           # Testing: 0 [cite: 46]
      estimate_projected_report = as.integer(horizon), # 0-3 horizon [cite: 50, 52]
      target_end_date = target_end_date
   )

# C. Combine and Apply Standard SW Formatting
sw_evaluation_file <- bind_rows(eval_state_training, eval_state_testing) %>%
   dplyr::mutate(
      # reference_date: Saturday of submission week [cite: 6, 8]
      reference_date   = format(latest_observed_week + weeks(1), "%Y-%m-%d"), 
      
      # target: Fixed character string [cite: 11]
      target           = "wk inc flu hosp",           
      
      # target_end_date: ISO format [cite: 16]
      target_end_date  = format(target_end_date, "%Y-%m-%d"),
      
      # location_general: Fixed for state level [cite: 19]
      location_general = "state",                    
      
      # location: "SC" for state model
      location         = "SC",                        
      
      # disease: influenza [cite: 35]
      disease          = "influenza",                 
      
      # population: general_population [cite: 38]
      population       = "general_population",        
      
      # imputed: Numeric 0 [cite: 58]
      imputed          = 0,                           
      
      # data_source: NA for evaluation model 
      data_source      = NA_character_,               
      
      # outcome_measure: Fixed label [cite: 70]
      outcome_measure  = "Weekly_Inpatient_Hospitalizations", 
      
      # output_type: Fixed as quantile [cite: 74]
      output_type      = "quantile",                  
      
      # output_type_id: Fixed at 0.5 for evaluation [cite: 80, 82]
      output_type_id   = 0.5                          
   ) %>%
   dplyr::arrange(target_end_date) %>%
   dplyr::select(
      reference_date, target, target_end_date, location_general, location, 
      value, disease, population, training_validation, 
      estimate_projected_report, imputed, data_source, 
      outcome_measure, output_type, output_type_id
   )

# D. Export Evaluation File
output_filename_eval <- paste0(
   "/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Forecasting Resources/Forecast-Drop-Off/Software/Evaluation/inpatient-hosp-DMAPRIME-QR-evaluation.csv"
)

write_csv(sw_evaluation_file, output_filename_eval)

cat(sprintf("✅ Evaluation File generated: %s\n", basename(output_filename_eval)))
