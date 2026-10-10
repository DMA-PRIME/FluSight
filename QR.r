

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
# 0. CONFIGURATION
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
eval_end   <- as.Date("2026-03-07")

# train_cut <- as.Date("2022-09-01")
# valid_end <- as.Date("2024-10-01")
# eval_start <- as.Date("2024-11-06")
# eval_end   <- as.Date("2025-03-07")

RUN_LEAKAGE_TEST       <- TRUE   # section 14a; stops the run on failure
LEAKAGE_TEST_N_ORIGINS <- 3

path_rfa <- "/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Data/SC Health Records/RFA/Years 2017-2025/Weekly data/RFA_weekly_influenza_region_incident.csv"
path_cdc <- '/Users/tanvirahammed/Downloads/target-hospital-admissions-22.csv'
path_musc <- "/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Data/MUSC/Infectious Disease EHR/Weekly Data/Latest Weekly Data/MUSC_Weekly_Influenza_State_dx_cond_lab_Incident.csv"
path_prisma <- "/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Data/Prisma Health/Infectious Disease EHR/Weekly Data/Latest Weekly Data/Prisma_Health_Weekly_Influenza_State_dx_cond_lab_Incident.csv"

output_dir_cdc <- "/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Tanvir/FluSight/Tanvir_QR/State/CDC_submission"
output_dir_software_impl <- "/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Forecasting Resources/Forecast-Drop-Off/Software/Implementation"
output_dir_software_eval <- "/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Forecasting Resources/Forecast-Drop-Off/Software/Evaluation"

target_var <- "cdc_Total.Influenza.Admissions"

TARGET_LOCATION <- "45"     # SC FIPS
TARGET_STATE    <- "SC"

#  "none"  = v4-LP behaviour. RUN THIS FIRST.
#  "pct"   = adds pct_states_rising.
EXTERNAL_SET <- "pct"

# ── v5.1 anti-lag switches (all off = v5-LP exactly) ─────────────────────────
DIRECT_BLEND   <- FALSE      # TRUE: blend recursive median with direct H1-H3 QR
EHR_MODE       <- "model"   # "freeze" (v5-LP) or "model" (forecast EHR forward)
PEAK_PULL      <- FALSE      # TRUE: analog peak timing pulls H1-H3 down
PATH_TUNE_STEP <- 2          # use every k-th validation origin for path tuning

# ── v5.2 switches (both off = v5.1 exactly) ──────────────────────────────────
SEASON_GATE      <- TRUE    # TRUE: peak/decline labels need an in-season level
SEASON_GATE_FRAC <- 0.25     # ... recent max >= this x typical season peak
POOL_STATES      <- TRUE    # TRUE: extra QR models trained on all states

if (!EHR_MODE %in% c("freeze", "model")) stop("EHR_MODE must be 'freeze' or 'model'")

.cfg_tag <- paste0(EXTERNAL_SET,
                   if (DIRECT_BLEND) "_direct" else "",
                   if (EHR_MODE == "model") "_ehrmodel" else "",
                   if (PEAK_PULL) "_peakpull" else "",
                   if (SEASON_GATE) "_seasongate" else "",
                   if (POOL_STATES) "_pool" else "")
MODEL_LABEL <- paste0("DMAPRIME-QR v5.2-LP [", gsub("_", " + ", .cfg_tag), "]")
MODEL_TAG   <- paste0("v52LP_", .cfg_tag, "_dec2025_mar2026")

# needed by the external-signal block; section 2 redefines it identically
safe_log <- function(x) log(pmax(x, 0) + 1)

# ══════════════════════════════════════════════════════════════════════════════
# 1. LOAD DATA - multi-state CDC panel + cross-state signal
# ══════════════════════════════════════════════════════════════════════════════

RFA_data_Influenza <- read_csv(path_rfa, show_col_types = FALSE)

MUSC_Weekly_Influenza   <- read_csv(path_musc,   show_col_types = FALSE)
Prisma_Weekly_Influenza <- read_csv(path_prisma, show_col_types = FALSE)

prep_source <- function(df, date_col, prefix) {
   df %>%
      rename(week = !!rlang::sym(date_col)) %>%
      mutate(week = as.Date(week)) %>%
      rename_with(~ paste0(prefix, .x), .cols = -week)
}

prisma <- prep_source(Prisma_Weekly_Influenza, "Week", "prisma_")
musc   <- prep_source(MUSC_Weekly_Influenza,   "Week", "musc_")

for (.src in c("prisma", "musc")) {
   if (anyDuplicated(get(.src)$week) > 0) stop(.src, " has duplicate weeks.")
}

# ── read the CDC file as a panel ─────────────────────────────────────────────

cdc_raw <- read_csv(path_cdc, show_col_types = FALSE)

.pick_col <- function(df, candidates) {
   hit <- candidates[candidates %in% names(df)]
   if (length(hit)) hit[1] else NA_character_
}

.col_week  <- .pick_col(cdc_raw, c("Week.Ending.Date", "Week Ending Date",
                                   "week_ending_date", "target_end_date",
                                   "date", "week"))
.col_flu   <- .pick_col(cdc_raw, c("Total.Influenza.Admissions",
                                   "Total Influenza Admissions",
                                   "total_influenza_admissions", "value"))
.col_loc   <- .pick_col(cdc_raw, c("location", "Location", "fips", "FIPS"))
.col_rate  <- .pick_col(cdc_raw, c("weekly_rate", "Weekly.Rate", "rate"))   # v5.2 pooling
.col_state <- .pick_col(cdc_raw, c("state", "State", "location_name",
                                   "Geographic.aggregation",
                                   "Geographic aggregation", "jurisdiction"))

if (is.na(.col_week)) stop("Cannot find a week column in: ",
                           paste(head(names(cdc_raw), 30), collapse = ", "))
if (is.na(.col_flu))  stop("Cannot find a flu-admissions column in: ",
                           paste(head(names(cdc_raw), 30), collapse = ", "))

.has_panel <- !is.na(.col_loc) || !is.na(.col_state)

cat(sprintf("[cdc] week col: '%s' | flu col: '%s' | geo col: '%s'\n",
            .col_week, .col_flu,
            if (.has_panel) paste(na.omit(c(.col_loc, .col_state)), collapse = "/")
            else "NONE - file appears pre-filtered"))

if (!all(c("Week.Ending.Date", "Total.Influenza.Admissions") %in% names(cdc_raw))) {
   cat("[cdc] NOTE: this is not the file v4 read; v4 scores are not directly comparable.\n")
   if (!is.na(.col_flu) && .col_flu == "value" && "weekly_rate" %in% names(cdc_raw)) {
      cat("[cdc]       Using 'value' (count); 'weekly_rate' is per 100k.\n")
   }
}

cdc_panel <- cdc_raw %>%
   transmute(
      week = as.Date(.data[[.col_week]]),
      value = as.numeric(.data[[.col_flu]]),
      rate  = if (!is.na(.col_rate)) as.numeric(.data[[.col_rate]]) else NA_real_,
      geo = if (!is.na(.col_loc)) as.character(.data[[.col_loc]])
      else if (!is.na(.col_state)) as.character(.data[[.col_state]])
      else TARGET_STATE
   ) %>%
   filter(!is.na(week), !is.na(geo)) %>%
   arrange(geo, week)

if (.has_panel) {
   cat(sprintf("[cdc] %d jurisdictions | %d rows | %s to %s\n",
               dplyr::n_distinct(cdc_panel$geo), nrow(cdc_panel),
               format(min(cdc_panel$week)), format(max(cdc_panel$week))))
}

.target_key <- if (TARGET_LOCATION %in% cdc_panel$geo) TARGET_LOCATION else
   if (TARGET_STATE %in% cdc_panel$geo) TARGET_STATE else
      if ("South Carolina" %in% cdc_panel$geo) "South Carolina" else
         NA_character_

if (is.na(.target_key)) {
   stop("Neither '", TARGET_LOCATION, "' nor '", TARGET_STATE,
        "' found in the geo column. Values present: ",
        paste(head(sort(unique(cdc_panel$geo)), 60), collapse = ", "))
}

# "US" matches ^[A-Z]{2}$ but is the national row - exclude aggregates
.aggregate_keys <- c("US", "USA", "National", "United States", "0", "00")

.is_state_key <- function(x) {
   (grepl("^[0-9]{1,2}$", x) | grepl("^[A-Z]{2}$", x)) & !(x %in% .aggregate_keys)
}

cat(sprintf("[cdc] target key: '%s'\n", .target_key))

if (.has_panel) {
   .keys <- sort(unique(cdc_panel$geo))
   .kept <- .keys[.is_state_key(.keys) & .keys != .target_key]
   .drop <- setdiff(.keys, c(.kept, .target_key))
   cat(sprintf("[cdc] %d keys counted in pct_states_rising\n", length(.kept)))
   if (length(.drop) > 0) {
      cat(sprintf("[cdc] %d excluded as aggregates/non-states: %s\n",
                  length(.drop), paste(head(.drop, 15), collapse = ", ")))
   }
   if (length(.kept) < 30) {
      cat("[cdc] WARNING: fewer than 30 states counted - pct_states_rising is noisy.\n")
   }
}

# helper: sum that keeps an all-missing week missing (the original turned it into 0)
.sum_or_na <- function(x) if (all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)

# ── target: SC, under the v4 column name ─────────────────────────────────────

CDC <- cdc_panel %>%
   filter(geo == .target_key) %>%
   group_by(week) %>%
   summarise(cdc_Total.Influenza.Admissions = .sum_or_na(value), .groups = "drop") %>%
   arrange(week)

cat(sprintf("[cdc] %s target: %d weeks | %s to %s\n",
            .target_key, nrow(CDC), format(min(CDC$week)), format(max(CDC$week))))

# ── external signal: share of OTHER states rising ────────────────────────────
# Uses each state's value at week t and week t-7 only (date-matched, not
# row-lagged, so a state's missing week cannot turn a 2-week change into a
# "1-week" one). Nothing after week t enters week t's value.

external_signals <- cdc_panel %>%
   filter(geo != .target_key, .is_state_key(geo)) %>%
   group_by(geo, week) %>%
   summarise(value = .sum_or_na(value), .groups = "drop") %>%
   group_by(geo) %>%
   mutate(
      prev_value = value[match(week - 7, week)],
      st_slope   = safe_log(value) - safe_log(prev_value)
   ) %>%
   ungroup() %>%
   group_by(week) %>%
   summarise(
      pct_states_rising  = if (all(is.na(st_slope))) NA_real_ else mean(st_slope > 0, na.rm = TRUE),
      n_states_reporting = sum(!is.na(st_slope)),
      .groups = "drop"
   ) %>%
   arrange(week)

ext_vars <- switch(
   EXTERNAL_SET,
   none = character(0),
   pct  = c("pct_states_rising"),
   stop("EXTERNAL_SET must be 'none' or 'pct'")
)

if (length(ext_vars) > 0 && !.has_panel) {
   stop("EXTERNAL_SET='pct' needs a multi-state CDC file, but no state or ",
        "location column was found in:\n  ", path_cdc,
        "\nSet EXTERNAL_SET <- \"none\" or point path_cdc at a panel file.")
}

if (length(ext_vars) > 0) {
   cat(sprintf("[cdc] pct_states_rising: %d weeks | median %.0f states reporting\n",
               nrow(external_signals), median(external_signals$n_states_reporting)))
}
cat(sprintf("[cdc] EXTERNAL_SET='%s' -> %d external predictors\n",
            EXTERNAL_SET, length(ext_vars)))

# ── gap-free weekly grid anchored on the SC target weeks (fix L5) ────────────

week_grid <- tibble(week = seq(min(CDC$week), max(CDC$week), by = "week"))
if (!all(CDC$week %in% week_grid$week)) {
   stop("CDC weeks do not fall on a regular 7-day grid (mixed weekdays?).")
}

.prisma_cov <- mean(week_grid$week %in% prisma$week)
cat(sprintf("[ehr] Prisma covers %.0f%% of CDC weeks\n", 100 * .prisma_cov))
if (.prisma_cov == 0) {
   stop("No Prisma week matches a CDC week - check that both use the same week-ending weekday.")
}

flu <- week_grid %>%
   left_join(CDC, by = "week") %>%
   left_join(dplyr::select(external_signals, week, all_of(ext_vars)), by = "week") %>%
   left_join(prisma, by = "week") %>%
   left_join(musc,   by = "week") %>%
   arrange(week)

optional_cols <- c(
   "prisma_Weekly_Tests",
   "prisma_Weekly_Positive_Tests",
   "prisma_Weekly_Inpatient_Hospitalizations"
)
for (nm in optional_cols) {
   if (!nm %in% names(flu)) flu[[nm]] <- NA_real_
}

if (anyDuplicated(flu$week) > 0) stop("Duplicate weeks in `flu` - a join multiplied rows.")

if (length(ext_vars) > 0) {
   .lost <- flu %>%
      filter(!is.na(.data[[target_var]]), is.na(.data[["pct_states_rising"]])) %>%
      nrow()
   if (.lost > 0)
      cat(sprintf("[cdc] NOTE: %d weeks have a target but no external signal - these drop as origins.\n",
                  .lost))
}

# ── v5.2 POOL_STATES: every other state as extra training data ──────────────

pool_flu_list <- list()
if (POOL_STATES) {
   if (!.has_panel) stop("POOL_STATES needs the multi-state CDC file.")
   if (is.na(.col_rate)) {
      stop("POOL_STATES needs a per-100k rate column (e.g. 'weekly_rate') to ",
           "rescale states to SC's population.")
   }
   .pop_src <- cdc_panel %>%
      filter(week < eval_start, !is.na(value), !is.na(rate), value > 0, rate > 0)
   pop_by_geo <- .pop_src %>%
      group_by(geo) %>%
      summarise(pop = median(value / rate * 1e5), .groups = "drop")
   .pop_sc <- pop_by_geo$pop[pop_by_geo$geo == .target_key]
   if (length(.pop_sc) != 1 || !is.finite(.pop_sc)) {
      stop("Cannot derive SC's population from the rate column.")
   }
   .pool_geos <- pop_by_geo %>%
      filter(geo != .target_key, .is_state_key(geo), is.finite(pop), pop > 0) %>%
      pull(geo)
   .skipped <- character(0)
   for (g in .pool_geos) {
      s_g <- cdc_panel %>%
         filter(geo == g) %>%
         group_by(week) %>%
         summarise(v = .sum_or_na(value), .groups = "drop") %>%
         arrange(week)
      if (nrow(s_g) < 60) { .skipped <- c(.skipped, g); next }
      grid_g <- tibble(week = seq(min(s_g$week), max(s_g$week), by = "week"))
      if (!all(s_g$week %in% grid_g$week)) { .skipped <- c(.skipped, g); next }
      fac <- .pop_sc / pop_by_geo$pop[pop_by_geo$geo == g]
      d_g <- grid_g %>%
         left_join(s_g, by = "week") %>%
         mutate(v = v * fac)
      names(d_g)[names(d_g) == "v"] <- target_var
      d_g <- d_g %>%
         left_join(dplyr::select(external_signals, week, all_of(ext_vars)), by = "week")
      for (nm in optional_cols) d_g[[nm]] <- NA_real_
      d_g$pool_geo <- g
      pool_flu_list[[g]] <- d_g
   }
   cat(sprintf("[pool] %d states pooled (rescaled to SC population %.0f)%s\n",
               length(pool_flu_list), .pop_sc,
               if (length(.skipped)) paste0(" | skipped: ", paste(.skipped, collapse = ", ")) else ""))
   if (length(pool_flu_list) < 10) {
      cat("[pool] WARNING: fewer than 10 states pooled - pooling adds little.\n")
   }
}

cat(sprintf("[data] flu: %d weeks (%s to %s) | %d cols\n",
            nrow(flu), format(min(flu$week)), format(max(flu$week)), ncol(flu)))

# ══════════════════════════════════════════════════════════════════════════════
# 2. GENERAL HELPERS  (unchanged, plus assert_weekly_grid)
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

# v5.2: flu seasons run August-July; season_id is the year the season starts.
season_id <- function(d) as.integer(year(d) - as.integer(month(d) < 8))

# Median peak (in counts) of seasons that ENDED before the origin's season,
# from data <= cutoff only. Seasons need >= 20 observed weeks to count.
typical_season_peak <- function(df, cutoff_week) {
   d <- df %>% filter(week <= cutoff_week, !is.na(.data[[target_var]]))
   if (nrow(d) == 0) return(NA_real_)
   pk <- d %>%
      mutate(.s = season_id(week)) %>%
      filter(.s < season_id(cutoff_week)) %>%
      group_by(.s) %>%
      summarise(n = n(), peak = max(.data[[target_var]]), .groups = "drop") %>%
      filter(n >= 20)
   if (nrow(pk) == 0) NA_real_ else median(pk$peak)
}

# LEAKAGE FIX L5: every lag()/lead() below is ROW-based. They only mean
# "previous week" / "next week" if rows are a gap-free 7-day sequence.
assert_weekly_grid <- function(df, where = "") {
   if (anyDuplicated(df$week) > 0)
      stop(where, ": duplicate weeks - a join multiplied rows.")
   if (nrow(df) > 1 && any(as.numeric(diff(sort(df$week))) != 7))
      stop(where, ": weeks are not a gap-free 7-day grid, so row-based lag()/",
           "lead() would not mean previous/next week.")
   invisible(TRUE)
}

# ══════════════════════════════════════════════════════════════════════════════
# 3. FEATURE ENGINEERING
# ══════════════════════════════════════════════════════════════════════════════

regime_flag_vars <- c("is_fast_decline", "any_decline", "any_peak_or_decline")

build_flu_features <- function(flu_df) {
   
   flu_df <- flu_df %>% arrange(week)
   assert_weekly_grid(flu_df, "build_flu_features")
   
   flu_df %>%
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
         
         # Row-local, backward-looking flags: safe to compute globally
         is_post_peak = if_else(
            weeks_since_12wk_peak >= 1 &
               log_target < max_12wk &
               slope_1wk < 0,
            1, 0
         ),
         
         is_near_peak = as.integer(
            recent_peak_ratio >= 0.80 &
               (accel_2wk <= 0 | slope_1wk <= 0.05)
         ),
         
         # Quantile-dependent flags: deliberately NA here (see L3)
         is_fast_decline     = NA_real_,
         any_decline         = NA_integer_,
         any_peak_or_decline = NA_integer_
      )
}

label_cols_for_h <- function(h) {
   paste0(c("future_log_h", "future_value_h", "future_change_h"), h)
}
peak_label_cols <- c("future_peak_offset_4wk",
                     "peak_by_1wk", "peak_by_2wk", "peak_by_3wk", "peak_by_4wk")

# Future labels. Built on the full series, but they never reach a model
# directly: make_origin_view() masks every label whose date is after the origin.
add_future_labels <- function(ff) {
   ff <- ff %>% arrange(week)
   assert_weekly_grid(ff, "add_future_labels")
   for (h in horizons) {
      ff[[paste0("future_log_h", h)]]    <- lead(ff$log_target, h + 1)
      ff[[paste0("future_value_h", h)]]  <- lead(ff[[target_var]], h + 1)
      ff[[paste0("future_change_h", h)]] <- ff[[paste0("future_log_h", h)]] - ff$log_target
   }
   ff %>%
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
}

build_all_features <- function(flu_df) add_future_labels(build_flu_features(flu_df))

flu_features <- build_all_features(flu)

# v5.2 POOL_STATES: features and labels for every pooled state, each built on
# its own weekly grid. Never used directly: make_pool_view() truncates and
# masks them per origin exactly like make_origin_view() does for SC.
build_pool_features <- function(lst) {
   if (length(lst) == 0) return(NULL)
   bind_rows(lapply(lst, build_all_features))
}
pool_features <- if (POOL_STATES) build_pool_features(pool_flu_list) else NULL

# LEAKAGE FIX L5: set every label that would only be known AFTER the origin to
# NA. Date-based (week + 7*(h+1)), so it holds even if a filter elsewhere is
# wrong. future_peak_offset_4wk looks at weeks t..t+4, so it needs week+28.
mask_future_labels <- function(df, origin_week) {
   origin_week <- as.Date(origin_week)
   for (h in horizons) {
      unknown <- df$week + 7 * (h + 1) > origin_week
      for (col in label_cols_for_h(h)) {
         if (col %in% names(df)) df[[col]][unknown] <- NA
      }
   }
   unknown_peak <- df$week + 28 > origin_week
   for (col in peak_label_cols) {
      if (col %in% names(df)) df[[col]][unknown_peak] <- NA
   }
   df
}

# ══════════════════════════════════════════════════════════════════════════════
# 4. PREDICTOR SETS AND PHASE LABELS
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

# External predictors (empty in v4; pct_states_rising in v5). Deliberately NOT
# added to state_vars, which drive analog matching.
predictors <- c(predictors, ext_vars)

cat(sprintf("[model] predictors: %d (%d external)\n",
            length(predictors), length(ext_vars)))

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

# Columns that must be non-NA for a week to be usable as an origin, excluding
# the per-origin flags (which are NA in the global flu_features by design).
origin_ready_vars <- setdiff(union(predictors, state_vars), regime_flag_vars)

# ══════════════════════════════════════════════════════════════════════════════
# THRESHOLDS, REGIME FLAGS, PHASE, AND THE ORIGIN VIEW
#
#  make_origin_view(features_df, origin_week) is the ONLY way data reaches a
#  model, a threshold, an analog pool or a forecast. It:
#    1. drops every row after origin_week
#    2. masks every label column dated after origin_week
#    3. computes thresholds from what is left
#    4. fills the quantile-dependent regime flags with those thresholds
#    5. assigns phase with those thresholds
#  flu_features itself carries NO phase and NO regime flags, so nothing can
#  silently use a full-sample version of them.
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
      vol_cut           = quantile(past$sd_4wk,    0.70, na.rm = TRUE, names = FALSE),
      slope_q25         = quantile(past$slope_1wk, 0.25, na.rm = TRUE, names = FALSE),
      global_slope_sd   = if (is.na(slope_sd) || slope_sd <= 0) 1 else slope_sd,
      # v5.1 PEAK_PULL inputs, origin-specific like everything else here
      post_peak_decline = {
         pp <- past$slope_1wk[!is.na(past$is_post_peak) & past$is_post_peak == 1 &
                                 past$slope_1wk < 0]
         if (length(pp) >= 5) median(-pp) else 0.15
      },
      high_level_cut    = quantile(past$log_target, 0.60, na.rm = TRUE, names = FALSE),
      # v5.2 SEASON_GATE input
      typical_peak      = typical_season_peak(df, cutoff_week)
   )
}

# LEAKAGE FIX L3: origin-specific version of the flags the original computed
# with a full-sample quantile.
assign_regime_flags <- function(df, thr, typical_peak = NULL) {
   q25 <- thr$slope_q25
   if (!SEASON_GATE) {
      return(df %>%
                mutate(
                   season_gate = 1L,
                   is_fast_decline = if_else(
                      slope_1wk < q25 & slope_2wk < 0 & drop_from_12wk_peak < 0, 1, 0),
                   any_decline = as.integer(
                      is_fast_decline == 1 | is_post_peak == 1 | decline_persistence_2wk == 1),
                   any_peak_or_decline = as.integer(is_near_peak == 1 | any_decline == 1)
                ))
   }
   # v5.2 SEASON_GATE: peak/decline labels only where the 8-week maximum is at
   # least SEASON_GATE_FRAC x the typical season peak (origin-specific; a
   # vector when pooled rows carry their own state's typical peak). Unknown
   # typical peak -> gate passes (= v5.1 behaviour).
   tp <- if (is.null(typical_peak)) rep(thr$typical_peak, nrow(df)) else typical_peak
   recent_max <- exp(df$max_8wk) - 1
   gate <- is.na(tp) | (!is.na(recent_max) & recent_max >= SEASON_GATE_FRAC * tp)
   df %>%
      mutate(
         season_gate  = as.integer(gate),
         is_near_peak = as.integer(ifelse(gate, is_near_peak, 0L)),
         is_post_peak = as.numeric(ifelse(gate, is_post_peak, 0)),
         is_fast_decline = if_else(
            gate & slope_1wk < q25 & slope_2wk < 0 & drop_from_12wk_peak < 0, 1, 0),
         any_decline = as.integer(
            is_fast_decline == 1 | is_post_peak == 1 |
               (gate & decline_persistence_2wk == 1)),
         any_peak_or_decline = as.integer(is_near_peak == 1 | any_decline == 1)
      )
}

# v5.2 POOL_STATES: per-origin view of the pooled states. Same rules as
# make_origin_view(): rows <= origin, labels dated after the origin masked,
# flags from origin thresholds (SC's), typical peak per state.
make_pool_view <- function(pool_df, origin_week, thr) {
   origin_week <- as.Date(origin_week)
   p <- pool_df %>% filter(week <= origin_week)
   if (nrow(p) == 0) return(NULL)
   p <- mask_future_labels(p, origin_week)
   tp <- NULL
   if (SEASON_GATE) {
      tpt <- p %>%
         filter(!is.na(.data[[target_var]])) %>%
         mutate(.s = season_id(week)) %>%
         filter(.s < season_id(origin_week)) %>%
         group_by(pool_geo, .s) %>%
         summarise(n = n(), peak = max(.data[[target_var]]), .groups = "drop") %>%
         filter(n >= 20) %>%
         group_by(pool_geo) %>%
         summarise(tp = median(peak), .groups = "drop")
      tp <- tpt$tp[match(p$pool_geo, tpt$pool_geo)]
   }
   assign_regime_flags(p, thr, typical_peak = tp)
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

make_origin_view <- function(features_df, origin_week) {
   origin_week <- as.Date(origin_week)
   past <- features_df %>%
      filter(week <= origin_week) %>%
      arrange(week)
   if (nrow(past) == 0 || max(past$week) != origin_week) return(NULL)
   past <- mask_future_labels(past, origin_week)
   thr  <- compute_origin_thresholds(past, origin_week)
   past <- past %>%
      assign_regime_flags(thr) %>%
      assign_phase(thr)
   list(df = past, thr = thr)
}

# NOTE: there is intentionally no global `global_slope_sd` and no global phase
# assignment on flu_features (original L2 / L4). Any stray reference errors.

# ══════════════════════════════════════════════════════════════════════════════
# 5. QR MODEL HELPERS  (unchanged)
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

# H0 is the only model used for the median; H1/H2/H3 apply it to the
# recursively updated row. Per-horizon fits are kept only for residuals.
fit_qr_model <- function(train_data, h, tau_value = 0.5, candidate = predictors) {
   target_h <- paste0("future_log_h", h)
   prepared <- prepare_qr_design(
      model_df = train_data,
      target_h = target_h,
      candidate_predictors = candidate
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

# ── v5.1 EHR_MODE = "model": 1-step models for the Prisma series ────────────
# Each series s gets a median regression of s(t+1) on its own level and
# 1-week change plus the CDC state (log_target, slopes) and season. Fitted on
# the origin view only: training rows need week + 7 <= origin, and lead()
# inside the truncated view is NA at the origin row anyway.
EHR_SERIES      <- c("log_hosp", "log_positive_tests", "log_tests")
EHR_FORMULA_RHS <- "lvl + d1 + log_target + slope_1wk + slope_2wk + sin52 + cos52"

fit_ehr_models <- function(view_df, origin_week) {
   out <- list()
   for (s in EHR_SERIES) {
      d <- view_df %>%
         arrange(week) %>%
         transmute(week,
                   lvl = .data[[s]],
                   d1  = .data[[s]] - lag(.data[[s]]),
                   y   = lead(.data[[s]]),
                   log_target, slope_1wk, slope_2wk, sin52, cos52) %>%
         filter(week + 7 <= origin_week) %>%
         drop_na()
      if (nrow(d) < 40) next
      fit <- tryCatch(
         suppressWarnings(rq(as.formula(paste("y ~", EHR_FORMULA_RHS)),
                             tau = 0.5, data = d, method = "br")),
         error = function(e) NULL)
      if (!is.null(fit)) {
         attr(fit, "y_range") <- range(d$y)
         out[[s]] <- fit
      }
   }
   out
}

# Next-week value of each modelled series, from the current (possibly
# synthetic) row. Clamped to the training range +/- 0.5 log units so a
# recursion can never run away. Series without a model are simply absent,
# which append_forecast_row() treats as frozen.
predict_ehr_next <- function(ehr_models, working_df, current_row) {
   if (is.null(ehr_models) || length(ehr_models) == 0) return(NULL)
   n   <- nrow(working_df)
   res <- list()
   for (s in names(ehr_models)) {
      lvl  <- current_row[[s]]
      prev <- if (n >= 2) working_df[[s]][n - 1] else NA_real_
      nd <- tibble(lvl = lvl, d1 = lvl - prev,
                   log_target = current_row$log_target,
                   slope_1wk  = current_row$slope_1wk,
                   slope_2wk  = current_row$slope_2wk,
                   sin52 = current_row$sin52, cos52 = current_row$cos52)
      p <- tryCatch(as.numeric(predict(ehr_models[[s]], newdata = nd)),
                    error = function(e) NA_real_)
      if (length(p) != 1 || is.na(p)) next
      r <- attr(ehr_models[[s]], "y_range")
      res[[s]] <- min(max(p, r[1] - 0.5), r[2] + 0.5)
   }
   res
}

# v5.2 POOL_STATES: direct H0-H3 models on SC + all pooled states. Prisma
# predictors are dropped (other states have none).
EHR_PREDICTORS  <- c("log_tests", "log_positive_tests", "log_hosp",
                     "hosp_slope_1wk", "hosp_slope_2wk",
                     "pos_slope_1wk", "pos_slope_2wk", "ehr_decline_confirmed")
predictors_pool <- setdiff(predictors, EHR_PREDICTORS)

fit_pool_models <- function(view, origin_week, pool_df) {
   if (is.null(pool_df)) return(NULL)
   pv <- make_pool_view(pool_df, origin_week, view$thr)
   if (is.null(pv)) return(NULL)
   train_all <- bind_rows(view$df %>% mutate(pool_geo = .target_key), pv)
   out <- list()
   for (h in horizons) {
      target_h <- paste0("future_log_h", h)
      md <- get_label_known_training(train_all, origin_week, h) %>%
         drop_na(all_of(c(target_h, predictors_pool)))
      if (nrow(md) < 80) next
      fit_h <- tryCatch(fit_qr_model(md, h, tau_value = 0.5, candidate = predictors_pool),
                        error = function(e) NULL)
      if (!is.null(fit_h)) out[[paste0("h", h)]] <- fit_h
   }
   out
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
# 6. ANALOG HELPERS  (unchanged except decline_analog_quantile_v4)
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
   # CEILING on QR weight during peak/decline phases so analogs get more say.
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

# LEAKAGE FIX L2: slope_sd is now an explicit argument. The original read a
# top-level `global_slope_sd`; the "local override" inside forecast_origin_v4
# never reached this function because R scopes lexically, so tuning used the
# all-data value and the backtest used the value through the latest week.
decline_analog_quantile_v4 <- function(row, h, peak_probs_row, slope_sd) {
   if (missing(slope_sd) || length(slope_sd) != 1 || is.na(slope_sd))
      stop("decline_analog_quantile_v4(): pass the origin-specific slope_sd (thr$global_slope_sd).")
   slope <- row$slope_1wk
   slope_severity <- pmin(abs(pmin(slope, 0)) / (slope_sd + 1e-8), 2) / 2
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
# 7. RECURSIVE ROW EXTENSION
# ══════════════════════════════════════════════════════════════════════════════

append_forecast_row <- function(df, new_week, pred_log, origin_row, thr,
                                ehr_next = NULL) {
   
   new_week <- as.Date(new_week)
   
   # LEAKAGE GUARD (L1): df must end exactly one week before new_week.
   if (nrow(df) == 0 || as.numeric(new_week - max(df$week)) != 7) {
      stop(sprintf(
         "append_forecast_row(): new_week %s is not one week after the last row %s. ",
         format(new_week), if (nrow(df)) format(max(df$week)) else "<empty>"),
         "df must not contain rows after the forecast origin.")
   }
   
   pred_count <- pmax(exp(pred_log) - 1, 0)
   
   # Template from the last row (<= current week by the guard above)
   new_row <- df %>% dplyr::slice(n())
   
   new_row$week           <- new_week
   new_row[[target_var]]  <- pred_count
   new_row$log_target     <- pred_log
   
   if (is.null(ehr_next)) {
      # EHR carry-forward (constant = last known) -> EHR slopes 0  (v5-LP)
      new_row$log_tests          <- origin_row$log_tests
      new_row$log_positive_tests <- origin_row$log_positive_tests
      new_row$log_hosp           <- origin_row$log_hosp
      new_row$hosp_slope_1wk     <- 0
      new_row$hosp_slope_2wk     <- 0
      new_row$pos_slope_1wk      <- 0
      new_row$pos_slope_2wk      <- 0
      new_row$ehr_decline_confirmed <- origin_row$ehr_decline_confirmed
   } else {
      # v5.1 EHR_MODE = "model": modelled next-week values; unmodelled series
      # stay frozen. Slopes are recomputed exactly as build_flu_features does.
      nd <- nrow(df)
      prev1 <- function(s) df[[s]][nd]
      prev2 <- function(s) if (nd >= 2) df[[s]][nd - 1] else NA_real_
      for (s in c("log_tests", "log_positive_tests", "log_hosp")) {
         new_row[[s]] <- if (!is.null(ehr_next[[s]])) ehr_next[[s]] else origin_row[[s]]
      }
      new_row$hosp_slope_1wk <- new_row$log_hosp - prev1("log_hosp")
      new_row$hosp_slope_2wk <- new_row$log_hosp - prev2("log_hosp")
      new_row$pos_slope_1wk  <- new_row$log_positive_tests - prev1("log_positive_tests")
      new_row$pos_slope_2wk  <- new_row$log_positive_tests - prev2("log_positive_tests")
      new_row$ehr_decline_confirmed <- as.integer(
         new_row$hosp_slope_1wk < 0 | new_row$pos_slope_1wk < 0 |
            new_row$hosp_slope_2wk < 0 | new_row$pos_slope_2wk < 0)
   }
   
   # External covariates (v5) held at their origin value; unknown in future.
   if (length(ext_vars) > 0) {
      for (.v in ext_vars) new_row[[.v]] <- origin_row[[.v]]
   }
   
   woy <- isoweek(new_week)
   new_row$week_of_year <- woy
   new_row$sin52        <- sin(2 * pi * woy / 52)
   new_row$cos52        <- cos(2 * pi * woy / 52)
   new_row$month        <- month(new_week)
   new_row$is_christmas <- if_else(month(new_week) == 12 & day(new_week) >= 20, 1, 0)
   new_row$is_newyear   <- if_else(month(new_week) == 1  & day(new_week) <= 7,  1, 0)
   
   # Observed history up to the origin + earlier predictions + this prediction
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
   
   new_row$is_post_peak <- if_else(
      new_row$weeks_since_12wk_peak >= 1 &
         pred_log < new_row$max_12wk &
         new_row$slope_1wk < 0, 1, 0)
   new_row$is_near_peak <- as.integer(
      new_row$recent_peak_ratio >= 0.80 &
         (new_row$accel_2wk <= 0 | new_row$slope_1wk <= 0.05))
   
   # Labels unknown for synthetic rows
   for (h_lab in horizons) {
      for (col in label_cols_for_h(h_lab)) new_row[[col]] <- NA_real_
   }
   new_row$future_peak_offset_4wk <- NA_integer_
   for (p in c("peak_by_1wk", "peak_by_2wk", "peak_by_3wk", "peak_by_4wk"))
      new_row[[p]] <- NA_integer_
   
   # Origin-specific regime flags and phase
   new_row <- new_row %>%
      assign_regime_flags(thr) %>%
      assign_phase(thr)
   
   bind_rows(df, new_row)
}

# ══════════════════════════════════════════════════════════════════════════════
# 8. FORECAST ENGINE - recursive loop
# ══════════════════════════════════════════════════════════════════════════════

# v5.1 output-stage options: direct-model weight per horizon (index h+1) and
# the peak-pull strength. make_opts() with defaults = v5-LP behaviour.
make_opts <- function(direct_w = rep(0, length(horizons)), peak_lambda = 0,
                      pool_w = rep(0, length(horizons))) {
   list(direct_w = direct_w, peak_lambda = peak_lambda, pool_w = pool_w)
}

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
                               quantiles,
                               ehr_models = NULL,
                               opts = make_opts(),
                               pool_models = NULL) {
   
   origin_week <- as.Date(origin_week)
   
   # LEAKAGE GUARD (L1/L5): whatever the caller passed, only rows <= origin
   # survive and no label dated after the origin is visible.
   df <- df %>% filter(week <= origin_week) %>% arrange(week)
   df <- mask_future_labels(df, origin_week)
   if (nrow(df) == 0 || max(df$week) != origin_week) return(NULL)
   if (!"phase" %in% names(df) || all(is.na(df$is_fast_decline))) {
      stop("forecast_origin_v4(): df must come from make_origin_view().")
   }
   
   required_cols <- union(predictors, state_vars)
   
   origin_row_0 <- df %>%
      dplyr::slice(n()) %>%
      drop_na(all_of(required_cols))
   if (nrow(origin_row_0) != 1) return(NULL)
   
   model_h0 <- qr_models[["h0"]]
   if (is.null(model_h0)) return(NULL)
   
   # Analog pool: past rows whose one-week-ahead change was known by the origin.
   # Built once from the past-only df; identical for every horizon.
   change_col  <- "future_change_h0"
   analog_pool <- get_label_known_training(df, origin_week, h = 0) %>%
      drop_na(all_of(c(state_vars, change_col)))
   
   working_df <- df
   
   # v5.1 DIRECT_BLEND: direct H1-H3 predictions from the REAL origin row
   dir_vec <- rep(NA_real_, length(horizons))
   for (hh in horizons[horizons >= 1]) {
      m_h <- qr_models[[paste0("h", hh)]]
      if (!is.null(m_h)) {
         dir_vec[hh + 1] <- tryCatch(predict_qr_safe(m_h, origin_row_0),
                                     error = function(e) NA_real_)
      }
   }
   # v5.2 POOL_STATES: pooled direct H0-H3 predictions from the REAL origin row
   pool_vec <- rep(NA_real_, length(horizons))
   if (!is.null(pool_models)) {
      for (hh in horizons) {
         m_p <- pool_models[[paste0("h", hh)]]
         if (!is.null(m_p)) {
            pool_vec[hh + 1] <- tryCatch(predict_qr_safe(m_p, origin_row_0),
                                         error = function(e) NA_real_)
         }
      }
   }
   # v5.1 PEAK_PULL: P(K <= k), k = 0..3, from the origin's analogs (set at h = 0)
   origin_peak_cdf <- rep(0, 4)
   
   median_rows          <- list()
   quantile_error_store <- list()
   analog_diag_store    <- list()
   
   for (h in horizons) {
      
      current_week_for_h <- origin_week + weeks(h)
      target_end_date_h  <- origin_week + weeks(h + 1)
      
      # ── 1. Current feature row: ALWAYS the last row of working_df ─────────
      # (The original filtered by week and took slice(1). With future rows in
      # df that returned the REAL observed row, not the synthetic one.)
      if (max(working_df$week) != current_week_for_h) {
         stop(sprintf("forecast_origin_v4(): expected last row %s at h=%d, found %s.",
                      format(current_week_for_h), h, format(max(working_df$week))))
      }
      current_row_raw <- working_df %>% dplyr::slice(n())
      
      # For synthetic rows: impute any NAs from origin_row_0
      if (h > 0) {
         na_cols <- required_cols[sapply(required_cols, function(v) {
            v %in% names(current_row_raw) && is.na(current_row_raw[[v]][1])
         })]
         for (col in na_cols) {
            if (col %in% names(origin_row_0))
               current_row_raw[[col]] <- origin_row_0[[col]]
         }
      }
      
      still_na <- required_cols[sapply(required_cols, function(v) {
         v %in% names(current_row_raw) && is.na(current_row_raw[[v]][1])
      })]
      
      if (length(still_na) > 0) {
         cat(sprintf("  [v4] h=%d: %d NAs remain after imputation (%s) - carry-forward\n",
                     h, length(still_na), paste(head(still_na, 3), collapse = ",")))
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
      
      qr_limit      <- decline_qr_floor_v4(current_row)
      decline_phase <- is_one(current_row$any_peak_or_decline)
      blend_qr_h    <- if (decline_phase) min(blend_qr_raw, qr_limit) else
         max(blend_qr_raw, qr_limit)
      blend_analog_h <- 1 - blend_qr_h
      
      # ── 3. QR prediction ──────────────────────────────────────────────────
      qr_median_log <- tryCatch(
         predict_qr_safe(model_h0, current_row),
         error = function(e) NA_real_
      )
      
      if (is.na(qr_median_log)) {
         cat(sprintf("  [v4] h=%d: QR predict failed - carry-forward\n", h))
         last_log <- current_row$log_target
         working_df <- append_forecast_row(working_df, target_end_date_h,
                                           last_log, current_row, thr)
         next
      }
      
      # ── 4. Analog selection ───────────────────────────────────────────────
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
         cat(sprintf("  [v4] h=%d: only %d analogs - QR-only fallback\n",
                     h, nrow(selected_analogs)))
      
      if (h == 0 && use_analog) {
         origin_peak_cdf <- vapply(0:3, function(k) {
            v <- weighted_mean_safe(selected_analogs$future_peak_offset_4wk <= k,
                                    selected_analogs$analog_weight)
            if (is.na(v)) 0 else v
         }, numeric(1))
      }
      
      # ── 5. Median forecast ────────────────────────────────────────────────
      peak_probs_h <- if (use_analog) {
         analog_peak_probs(selected_analogs)
      } else {
         tibble(prob_peak_by_1 = NA_real_, prob_peak_by_2 = NA_real_,
                prob_peak_by_3 = NA_real_, prob_peak_by_4 = NA_real_)
      }
      
      if (use_analog) {
         analog_change <- if (is_one(current_row$any_peak_or_decline)) {
            dq <- decline_analog_quantile_v4(current_row, h, peak_probs_h,
                                             slope_sd = thr$global_slope_sd)
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
      
      # ── 8. Extend working_df with the prediction ──────────────────────────
      ehr_next <- if (!is.null(ehr_models)) {
         predict_ehr_next(ehr_models, working_df, current_row)
      } else NULL
      
      working_df <- append_forecast_row(
         df         = working_df,
         new_week   = target_end_date_h,
         pred_log   = raw_combined_log,
         origin_row = current_row,
         thr        = thr,
         ehr_next   = ehr_next
      )
      
      cat(sprintf("  [v4] h=%d done: pred_log=%.3f  phase=%s  use_analog=%s\n",
                  h, raw_combined_log, current_row$phase, use_analog))
   }
   
   median_df <- bind_rows(median_rows)
   if (nrow(median_df) == 0) return(NULL)
   
   # ── v5.1 output stage: direct blend and peak pull ──────────────────────
   # The recursion above is untouched (rec_log). Final median for h >= 1:
   #   (1 - w_h) * rec + w_h * direct  -  lambda * E[weeks past peak] * decline
   # With make_opts() defaults this equals rec_log exactly (v5-LP).
   in_season <- if (SEASON_GATE) {
      isTRUE(origin_row_0$season_gate == 1)          # v5.2: season-relative
   } else {
      isTRUE(origin_row_0$log_target >= thr$high_level_cut)
   }
   pull_unit_vec <- vapply(horizons, function(hh) {
      if (hh >= 1 && in_season) sum(origin_peak_cdf[seq_len(hh + 1)]) * thr$post_peak_decline else 0
   }, numeric(1))
   
   median_df <- median_df %>%
      mutate(
         rec_log   = raw_corrected_log,
         dir_log   = dir_vec[horizon + 1],
         pull_unit = pull_unit_vec[horizon + 1],
         w_direct  = opts$direct_w[horizon + 1],
         pool_log  = pool_vec[horizon + 1],
         w_pool    = opts$pool_w[horizon + 1],
         final_log = pmax((1 - w_direct - w_pool) * rec_log +
                             w_direct * coalesce(dir_log, rec_log) +
                             w_pool * coalesce(pool_log, rec_log) -
                             opts$peak_lambda * pull_unit, 0)
      )
   
   error_df <- bind_rows(quantile_error_store)
   
   all_quantiles <- median_df %>%
      dplyr::select(
         origin_week, target_end_date, horizon, phase,
         final_log, raw_corrected_log, qr_median_log, analog_target_log,
         blend_qr_raw, blend_qr, blend_analog, k, bandwidth,
         rec_log, dir_log, pull_unit, w_direct, pool_log, w_pool,
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
# 9. ROLLING VALIDATION / TUNING  (h=0, reused recursively for all horizons)
#
#  Every origin now goes through make_origin_view(): QR training rows, the
#  origin row, AND the analog pool all come from the past-only view (the
#  original drew the pool from the globally phase-labelled flu_features).
# ══════════════════════════════════════════════════════════════════════════════

validation_origins <- flu_features %>%
   filter(week >= train_cut, week < valid_end, !is.na(log_target)) %>%
   drop_na(all_of(origin_ready_vars)) %>%
   pull(week)

if (length(validation_origins) == 0) stop("No validation origins in [train_cut, valid_end).")

# The tuning labels (origin + 1 week) must all precede the evaluation window,
# otherwise hyperparameters are chosen with evaluation-period outcomes.
if (max(validation_origins) + weeks(1) >= eval_start) {
   stop("Tuning window overlaps evaluation: last tuning label ",
        format(max(validation_origins) + weeks(1)),
        " is not before eval_start ", format(eval_start), ".")
}

cat("\nTuning h=0 (1-step) - used recursively for all horizons ...\n")

pred_store <- list()
change_col <- "future_change_h0"
target_h   <- "future_log_h0"

for (origin_week in validation_origins) {
   origin_week <- as.Date(origin_week)
   
   view_tune <- make_origin_view(flu_features, origin_week)
   if (is.null(view_tune)) next
   ov       <- view_tune$df
   thr_tune <- view_tune$thr
   
   train_origin <- get_label_known_training(ov, origin_week, 0) %>%
      drop_na(all_of(c(target_h, predictors)))
   
   origin_row <- ov %>%
      filter(week == origin_week) %>%
      drop_na(all_of(c(predictors, state_vars)))
   
   # Scoring label only - read from the unmasked series, never fed to a model.
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
   
   analog_pool <- get_label_known_training(ov, origin_week, 0) %>%
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
         dq <- decline_analog_quantile_v4(origin_row, 0, peak_probs_g,
                                          slope_sd = thr_tune$global_slope_sd)
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

tuning_table <- bind_rows(lapply(horizons, function(h) {
   tibble(
      horizon   = h,
      k         = best_h0$k,
      bandwidth = best_h0$bandwidth,
      blend_qr  = pmax(best_h0$blend_qr - h * 0.05, 0.10),
      mae_log   = best_h0$mae_log,
      n_valid   = best_h0$n_valid
   )
}))

cat("\n═══ TUNING TABLE ═══\n")
print(tuning_table)

# ══════════════════════════════════════════════════════════════════════════════
# 9b. v5.1 PIPELINE HELPERS + PATH TUNING (DIRECT_BLEND weight, PEAK_PULL lambda)
# ══════════════════════════════════════════════════════════════════════════════

# Every forecast (tuning paths, final, leakage test, backtest) goes through
# these two functions, so all of them use identical machinery.
fit_origin_models <- function(view, origin_week, pool_df = pool_features) {
   m <- fit_all_qr_models(view$df, origin_week)
   m$ehr_models  <- if (EHR_MODE == "model") fit_ehr_models(view$df, origin_week) else NULL
   m$pool_models <- if (POOL_STATES) fit_pool_models(view, origin_week, pool_df) else NULL
   m
}

run_forecast_at <- function(view, origin_week, opts, pool_df = pool_features) {
   m <- fit_origin_models(view, origin_week, pool_df)
   forecast_origin_v4(
      origin_week  = origin_week,
      df           = view$df,
      qr_models    = m$qr_models,
      residual_df  = m$residual_df,
      tuning_table = tuning_table,
      thr          = view$thr,
      quantiles    = quantiles,
      ehr_models   = m$ehr_models,
      opts         = opts,
      pool_models  = m$pool_models
   )
}

opts_final   <- make_opts()
path_tune_df <- NULL

if (DIRECT_BLEND || PEAK_PULL || POOL_STATES) {
   
   path_origins <- validation_origins[seq(1, length(validation_origins), by = PATH_TUNE_STEP)]
   
   # Every path target (origin + 4 weeks at H3) must precede the eval window
   if (max(path_origins) + weeks(max(horizons) + 1) >= eval_start) {
      stop("Path tuning overlaps evaluation: last H3 target ",
           format(max(path_origins) + weeks(max(horizons) + 1)),
           " is not before eval_start ", format(eval_start), ".")
   }
   
   cat(sprintf("\nPath tuning on %d validation origins (every %d) ...\n",
               length(path_origins), PATH_TUNE_STEP))
   
   path_rows <- list()
   for (o in path_origins) {
      o <- as.Date(o)
      v <- make_origin_view(flu_features, o)
      if (is.null(v)) next
      fc <- run_forecast_at(v, o, make_opts())
      if (is.null(fc)) next
      # Scoring labels only - read from the unmasked series after forecasting
      truth_o <- flu_features %>% filter(week == o)
      path_rows[[length(path_rows) + 1]] <- fc$median_df %>%
         transmute(origin_week, horizon, rec_log, dir_log, pull_unit, pool_log,
                   y_true = vapply(horizon, function(hh)
                      truth_o[[paste0("future_log_h", hh)]][1], numeric(1)))
   }
   path_tune_df <- bind_rows(path_rows) %>% filter(!is.na(y_true))
   
   if (nrow(path_tune_df) == 0) {
      cat("[v5.1] no usable validation paths - switches left at neutral values\n")
   } else {
      w_grid  <- if (DIRECT_BLEND) seq(0, 1, by = 0.25) else 0
      wp_grid <- if (POOL_STATES)  seq(0, 1, by = 0.25) else 0
      l_grid  <- if (PEAK_PULL) c(0, 0.5, 1, 1.5, 2) else 0
      
      # H0 has no SC direct model and no pull, so only the pool weight can
      # change it; the tie-break below keeps H0 = recursive when pooling is off.
      grid_scores <- path_tune_df %>%
         crossing(w = w_grid, wp = wp_grid, lambda = l_grid) %>%
         filter(w + wp <= 1 + 1e-9) %>%
         mutate(
            pred    = pmax((1 - w - wp) * rec_log +
                              w  * coalesce(dir_log, rec_log) +
                              wp * coalesce(pool_log, rec_log) -
                              lambda * pull_unit, 0),
            abs_err = abs(y_true - pred)
         ) %>%
         group_by(horizon, w, wp, lambda) %>%
         summarise(mae_log = mean(abs_err), n = n(), .groups = "drop")
      
      # For each lambda, best (w, wp) per horizon (ties -> less blending);
      # then the lambda with the lowest total MAE (ties -> smaller lambda).
      best_by_lambda <- grid_scores %>%
         group_by(lambda, horizon) %>%
         arrange(mae_log, w + wp, w, .by_group = TRUE) %>%
         dplyr::slice(1) %>%
         ungroup()
      lambda_best <- best_by_lambda %>%
         group_by(lambda) %>%
         summarise(total = sum(mae_log), .groups = "drop") %>%
         arrange(total, lambda) %>%
         dplyr::slice(1) %>%
         pull(lambda)
      chosen <- best_by_lambda %>% filter(lambda == lambda_best)
      
      w_vec  <- rep(0, length(horizons))
      wp_vec <- rep(0, length(horizons))
      for (i in seq_len(nrow(chosen))) {
         w_vec[chosen$horizon[i] + 1]  <- chosen$w[i]
         wp_vec[chosen$horizon[i] + 1] <- chosen$wp[i]
      }
      opts_final <- make_opts(direct_w = w_vec, peak_lambda = lambda_best, pool_w = wp_vec)
      
      cat("\n═══ PATH TUNING RESULT (validation MAE, log scale) ═══\n")
      print(
         grid_scores %>%
            filter(w == 0, wp == 0, lambda == 0) %>%
            dplyr::select(horizon, mae_recursive_only = mae_log, n) %>%
            left_join(chosen %>% dplyr::select(horizon, w_direct = w, w_pool = wp,
                                               mae_chosen = mae_log), by = "horizon") %>%
            mutate(peak_lambda = lambda_best)
      )
   }
}

cat(sprintf("[v5.2] final options: direct_w = (%s) | pool_w = (%s) | peak_lambda = %.2f | EHR_MODE = %s | SEASON_GATE = %s\n",
            paste(sprintf("%.2f", opts_final$direct_w), collapse = ", "),
            paste(sprintf("%.2f", opts_final$pool_w), collapse = ", "),
            opts_final$peak_lambda, EHR_MODE, SEASON_GATE))

# ══════════════════════════════════════════════════════════════════════════════
# 10. FINAL FORECAST
# ══════════════════════════════════════════════════════════════════════════════

latest_observed_week <- max(
   flu_features$week[!is.na(flu_features[[target_var]])],
   na.rm = TRUE
)
reference_date <- latest_observed_week + weeks(1)

final_view <- make_origin_view(flu_features, latest_observed_week)
if (is.null(final_view)) stop("Could not build the origin view for the latest observed week.")
thr_final <- final_view$thr

current_row <- final_view$df %>%
   filter(week == latest_observed_week) %>%
   drop_na(all_of(c(predictors, state_vars)))

if (nrow(current_row) != 1) {
   stop("Current row is incomplete. Check missing predictors for latest observed week.")
}

observed_value <- current_row[[target_var]]
phase <- current_row$phase

final_models <- fit_origin_models(final_view, latest_observed_week)

final_forecast <- forecast_origin_v4(
   origin_week  = latest_observed_week,
   df           = final_view$df,
   qr_models    = final_models$qr_models,
   residual_df  = final_models$residual_df,
   tuning_table = tuning_table,
   thr          = thr_final,
   quantiles    = quantiles,
   ehr_models   = final_models$ehr_models,
   opts         = opts_final,
   pool_models  = final_models$pool_models
)

if (is.null(final_forecast)) {
   stop("Final forecast failed. Check model completeness and analog availability.")
}

all_quantiles         <- final_forecast$all_quantiles
analog_diagnostics_df <- final_forecast$analog_diagnostics
median_forecast_df    <- final_forecast$median_df

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

cat("\n═══ FINAL FORECAST DIAGNOSTICS ═══\n")
cat(sprintf("Latest observed week : %s\n", latest_observed_week))
cat(sprintf("Latest observed value: %s\n", observed_value))
cat(sprintf("Detected phase       : %s\n", phase))
cat(sprintf("Typical season peak  : %s  | season gate at origin: %s\n",
            if (is.na(thr_final$typical_peak)) "NA (no completed season)" else format(round(thr_final$typical_peak)),
            if (SEASON_GATE) (if (isTRUE(current_row$season_gate == 1)) "open (in season)" else "closed (peak/decline labels off)") else "switch off"))

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
# 11. FINAL FORECAST PLOT  (unchanged)
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
      title    = paste0(MODEL_LABEL, " Recursive Forecast [", phase, "]"),
      subtitle = paste0(
         "Latest observed: ", format(latest_observed_week, "%Y-%m-%d"),
         " (", observed_value, " admissions)"
      ),
      y       = "Weekly Inpatient Hospitalizations (Influenza, South Carolina)",
      x       = "Week Ending Date",
      caption = "Recursive 1-step QR + phase-aligned analogs | origin-only data (leakage-proof)"
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
# 12. EXPORT CDC SUBMISSION  (unchanged)
# ══════════════════════════════════════════════════════════════════════════════

output_filename_cdc <- file.path(
   output_dir_cdc,
   paste0(format(reference_date, "%Y-%m-%d"), "-DMAPRIME-QR.csv")
)
# Uncomment to save:
#write_csv(cdc_submission, output_filename_cdc)

cat("\n✅ FINAL FORECAST COMPLETE\n")
cat("CDC file: ", basename(output_filename_cdc), "\n")

# 13. SOFTWARE IMPLEMENTATION FILE  (unchanged) ------------------------------

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

sw_implementation_file <- bind_rows(
   impl_state_rfa %>% dplyr::mutate(reference_date = .env$reference_date),
   impl_state_hs %>% dplyr::mutate(reference_date = .env$reference_date),
   impl_state_forecast %>% dplyr::mutate(reference_date = .env$reference_date)
) %>%
   dplyr::mutate(
      reference_date = reference_date,
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
   "inpatient-hosp-DMAPRIME-QR-implementation.csv"
)
# Uncomment to save:
#write_csv(sw_implementation_file, output_filename_software)
cat("Software implementation file: ", basename(output_filename_software), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# 14. RETROSPECTIVE EVALUATION
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
   drop_na(all_of(origin_ready_vars)) %>%
   pull(week)

if (length(eval_origin_weeks) == 0) stop("No evaluation origins in [eval_start, eval_end].")

# ══════════════════════════════════════════════════════════════════════════════
# 14a. LEAKAGE PROOF - FUTURE-PERTURBATION TEST
#
#  If a forecast issued at origin T uses only data <= T, then overwriting every
#  raw value AFTER T with noise (including filling future NAs) and rebuilding
#  ALL features from scratch must leave that forecast unchanged to machine
#  precision. This catches bad row filtering, unmasked labels, AND full-sample
#  statistics inside feature engineering (the original L1, L2 and L3 all fail
#  it). Tuning is held fixed: tuning_table was chosen on data before
#  eval_start, which section 9 asserts.
# ══════════════════════════════════════════════════════════════════════════════

run_origin_forecast <- function(features_df, origin_week, pool_df = pool_features) {
   v <- make_origin_view(features_df, origin_week)
   if (is.null(v)) return(NULL)
   fc <- run_forecast_at(v, origin_week, opts_final, pool_df = pool_df)   # same pipeline as backtest
   if (is.null(fc)) return(NULL)
   fc$all_quantiles %>% arrange(horizon, quantile)
}

# Overwrite every numeric value dated after the origin with noise (and fill
# future NAs), leaving the past untouched.
perturb_future <- function(d, origin_week) {
   fut <- d$week > origin_week
   if (!any(fut)) return(d)
   num_cols <- names(d)[vapply(d, is.numeric, logical(1))]
   for (cc in num_cols) {
      x <- as.numeric(d[[cc]])
      xf <- x[fut]
      xf <- ifelse(is.na(xf), runif(length(xf), 0, 500),
                   xf * runif(length(xf), 0.1, 10) + runif(length(xf), 0, 500))
      x[fut] <- xf
      d[[cc]] <- x
   }
   d
}

leakage_perturbation_test <- function(origin_week, raw_flu = flu, seed = 2026) {
   origin_week <- as.Date(origin_week)
   set.seed(seed)
   if (!any(raw_flu$week > origin_week)) return(NA)
   flu_bad  <- perturb_future(raw_flu, origin_week)
   # v5.2: the pooled states' future is corrupted too, then rebuilt from scratch
   pool_bad <- if (POOL_STATES) {
      build_pool_features(lapply(pool_flu_list, perturb_future, origin_week = origin_week))
   } else NULL
   a <- run_origin_forecast(flu_features, origin_week)
   b <- run_origin_forecast(build_all_features(flu_bad), origin_week, pool_df = pool_bad)
   if (is.null(a) && is.null(b)) return(NA)
   if (is.null(a) != is.null(b)) return(FALSE)
   identical(a$horizon, b$horizon) &&
      identical(a$phase, b$phase) &&
      isTRUE(all.equal(a$predicted_value, b$predicted_value, tolerance = 1e-10)) &&
      isTRUE(all.equal(a$final_log,       b$final_log,       tolerance = 1e-10))
}

if (RUN_LEAKAGE_TEST) {
   n_test <- min(LEAKAGE_TEST_N_ORIGINS, length(eval_origin_weeks))
   test_origins <- unique(eval_origin_weeks[
      round(seq(1, length(eval_origin_weeks), length.out = n_test))])
   cat(sprintf("\n═══ LEAKAGE PERTURBATION TEST (%d origins) ═══\n", length(test_origins)))
   leak_results <- vapply(test_origins, function(o) {
      res <- leakage_perturbation_test(o)
      cat(sprintf("  origin %s : %s\n", format(as.Date(o)),
                  if (is.na(res)) "skipped (no forecast)" else if (res) "PASS" else "FAIL"))
      res
   }, logical(1))
   if (any(leak_results %in% FALSE)) {
      stop("LEAKAGE DETECTED: corrupting data after the origin changed a forecast.")
   }
   cat("  No future information reaches any tested forecast.\n")
}

cat(sprintf(
   "\nRetrospective evaluation over %d origin weeks (%s to %s)\n",
   length(eval_origin_weeks),
   format(min(eval_origin_weeks), "%Y-%m-%d"),
   format(max(eval_origin_weeks), "%Y-%m-%d")
))

retro_rows <- list()

for (origin_week in eval_origin_weeks) {
   origin_week <- as.Date(origin_week)
   
   # Past-only view: rows <= origin, future labels masked, origin thresholds,
   # origin regime flags, origin phases. (The original kept rows to origin+4.)
   view_retro <- make_origin_view(flu_features, origin_week)
   if (is.null(view_retro)) next
   
   fc_retro <- run_forecast_at(view_retro, origin_week, opts_final)
   
   if (is.null(fc_retro)) next
   
   med_retro <- fc_retro$all_quantiles %>%
      filter(quantile == 0.5) %>%
      mutate(forecast_median = round(predicted_value, 1)) %>%
      dplyr::select(
         origin_week, target_end_date, horizon, forecast_median,
         phase, raw_corrected_log, final_log, qr_median_log, analog_target_log,
         blend_qr, prob_peak_by_1, prob_peak_by_2, prob_peak_by_3, prob_peak_by_4,
         rec_log, dir_log, pull_unit, w_direct, pool_log, w_pool
      )
   
   # Interval quantiles for coverage checks
   q_retro <- fc_retro$all_quantiles %>%
      filter(quantile %in% c(0.025, 0.25, 0.75, 0.975)) %>%
      dplyr::select(horizon, quantile, predicted_value) %>%
      mutate(quantile = paste0("q", quantile)) %>%
      pivot_wider(names_from = quantile, values_from = predicted_value)
   med_retro <- med_retro %>% left_join(q_retro, by = "horizon")
   
   # Truth lookup for SCORING ONLY (after the forecast is fixed)
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

cat("\n═══ OVERALL PERIOD SUMMARY ═══\n")
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

cat("\n═══ PERFORMANCE BY HORIZON ═══\n")
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

cat("\n═══ INTERVAL COVERAGE (nominal 50% and 95%) ═══\n")
coverage_summary <- retro_df %>%
   group_by(horizon) %>%
   summarise(
      cover_50 = round(100 * mean(actual >= q0.25  & actual <= q0.75,  na.rm = TRUE), 1),
      cover_95 = round(100 * mean(actual >= q0.025 & actual <= q0.975, na.rm = TRUE), 1),
      .groups = "drop"
   )
print(coverage_summary)

cat("\n═══ VS. PERSISTENCE (forecast = value observed at origin) ═══\n")
persistence_summary <- retro_df %>%
   left_join(flu_features %>% dplyr::select(week, origin_value = all_of(target_var)),
             by = c("origin_week" = "week")) %>%
   group_by(horizon) %>%
   summarise(
      model_ai       = round(mean(ai_score, na.rm = TRUE), 4),
      persistence_ai = round(mean(ai_score(actual, origin_value), na.rm = TRUE), 4),
      .groups = "drop"
   ) %>%
   mutate(beats_persistence = model_ai > persistence_ai)
print(persistence_summary)

cat("\n═══ TIMING LAG (weeks the forecast trails the observed curve) ═══\n")
lag_summary <- retro_df %>%
   dplyr::select(horizon, target_end_date, forecast_median) %>%
   crossing(shift = 0:5) %>%
   mutate(week = target_end_date - weeks(shift)) %>%
   left_join(flu_features %>% dplyr::select(week, obs = all_of(target_var)), by = "week") %>%
   group_by(horizon, shift) %>%
   summarise(r = cor(log1p(forecast_median), log1p(obs), use = "complete.obs"),
             .groups = "drop") %>%
   group_by(horizon) %>%
   slice_max(r, n = 1, with_ties = FALSE) %>%
   ungroup() %>%
   rename(best_shift_wk = shift)
print(lag_summary)
cat("  0 = on time; close to h+1 = the forecast mostly carries the origin level forward.\n")

# One-row summary per configuration - append across runs to compare switches
config_summary <- tibble(
   run_time      = format(Sys.time(), "%Y-%m-%d %H:%M"),
   config        = .cfg_tag,
   eval_window   = paste(eval_start, eval_end, sep = " to "),
   direct_w      = paste(sprintf("%.2f", opts_final$direct_w), collapse = "/"),
   pool_w        = paste(sprintf("%.2f", opts_final$pool_w), collapse = "/"),
   peak_lambda   = opts_final$peak_lambda,
   mean_ai       = round(mean(retro_df$ai_score, na.rm = TRUE), 4)
) %>%
   bind_cols(
      horizon_summary %>%
         dplyr::select(horizon, mean_ai_score) %>%
         mutate(name = paste0("ai_h", horizon)) %>%
         dplyr::select(name, mean_ai_score) %>%
         pivot_wider(names_from = name, values_from = mean_ai_score),
      coverage_summary %>%
         summarise(cover_50 = mean(cover_50), cover_95 = mean(cover_95)),
      lag_summary %>%
         mutate(name = paste0("lag_h", horizon)) %>%
         dplyr::select(name, best_shift_wk) %>%
         pivot_wider(names_from = name, values_from = best_shift_wk)
   )
cat("\n═══ CONFIG SUMMARY ═══\n")
print(as.data.frame(config_summary))

config_log_file <- file.path(output_dir_cdc, "v52LP_config_comparison.csv")
# Uncomment to append this run to the comparison log:
# write_csv(config_summary, config_log_file, append = file.exists(config_log_file))

cat("\n═══ WEEK-BY-WEEK DETAIL ═══\n")
week_detail <- retro_df %>%
   dplyr::select(
      origin_week, target_end_date, horizon, actual, forecast_median,
      ai_score, score_band, phase, raw_corrected_log, final_log,
      blend_qr, prob_peak_by_1, prob_peak_by_2, prob_peak_by_3
   ) %>%
   arrange(origin_week, horizon)
print(week_detail, n = Inf)

# ══════════════════════════════════════════════════════════════════════════════
# 15. RETROSPECTIVE VISUALIZATIONS  (unchanged apart from titles)
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
      title    = paste0("State-Level: CDC Observed vs Forecast Diagnostic - ", MODEL_LABEL),
      subtitle = "Colored lines are median forecasts by horizon (origin-only data)",
      x = "Week Ending Date",
      y = "Weekly Inpatient Hospitalizations (Influenza, SC)",
      caption = paste0(
         "Black = CDC observed | Colored = median forecast by horizon\n",
         "Each horizon uses features updated with the prior horizon's prediction, never observed future data."
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
   labs(title    = paste0("Retrospective AI Scores by Horizon - ", MODEL_LABEL),
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
   labs(title    = paste0("Actual vs. Forecast - Retrospective Evaluation ", MODEL_LABEL),
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
   labs(title    = paste0("Mean AI Score by Forecast Horizon - ", MODEL_LABEL),
        subtitle = "Recursive 1-step QR + analog forecast",
        x = NULL, y = "Mean AI score") +
   theme_minimal(base_size = 11) +
   theme(legend.position = "bottom", panel.grid.minor = element_blank())
print(p_bar)

# ══════════════════════════════════════════════════════════════════════════════
# 16. EXPORT RETROSPECTIVE + SOFTWARE EVALUATION  (unchanged)
# ══════════════════════════════════════════════════════════════════════════════

output_retro_file <- file.path(
   output_dir_cdc, paste0("retrospective_ai_scores_", MODEL_TAG, ".csv")
)
output_horizon_summary_file <- file.path(
   output_dir_cdc, paste0("retrospective_ai_scores_", MODEL_TAG, "_by_horizon_summary.csv")
)
# Uncomment to save:
# write_csv(retro_df, output_retro_file)
# write_csv(horizon_summary, output_horizon_summary_file)

eval_state_training <- flu_features %>%
   dplyr::filter(week <= train_cut & !is.na(.data[[target_var]])) %>%
   dplyr::mutate(
      value = .data[[target_var]], training_validation = 1,
      target_end_date = week
   )

eval_state_validation <- flu_features %>%
   dplyr::filter(
      week >= train_cut &
         week <= valid_end &
         !is.na(.data[[target_var]])
   ) %>%
   dplyr::mutate(
      value = .data[[target_var]], training_validation = 2,
      target_end_date = week
   )

eval_state_testing <- retro_df %>%
   dplyr::filter(
      origin_week >= eval_start &
         origin_week <= eval_end &
         !is.na(.data[["forecast_median"]])
   ) %>%
   dplyr::mutate(
      value = .data[["forecast_median"]],
      training_validation = 0,
      target_end_date = target_end_date
   )

sw_evaluation_file <- bind_rows(
   eval_state_training   %>% dplyr::mutate(reference_date = .env$reference_date),
   eval_state_validation %>% dplyr::mutate(reference_date = .env$reference_date),
   eval_state_testing    %>% dplyr::mutate(reference_date = .env$reference_date)
) %>%
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
      value, disease, population, training_validation, horizon,
      imputed, data_source, outcome_measure, output_type, output_type_id
   )

output_filename_eval <- file.path(
   output_dir_software_eval,
   "inpatient-hosp-DMAPRIME-QR-evaluation.csv"
)
# Uncomment to save:
#write_csv(sw_evaluation_file, output_filename_eval)

cat("Software evaluation file: ", basename(output_filename_eval), "\n")

cat("\n✅ RETROSPECTIVE EVALUATION COMPLETE\n")
cat(sprintf("Rows evaluated : %d\n", nrow(retro_df)))
cat(sprintf("Origin weeks   : %d\n", n_distinct(retro_df$origin_week)))
cat(sprintf("Overall mean AI: %.4f\n", mean(retro_df$ai_score, na.rm = TRUE)))

# ══════════════════════════════════════════════════════════════════════════════
# 17. ORIGIN-PATH DIAGNOSTIC  (unchanged)
# ══════════════════════════════════════════════════════════════════════════════

actual_line <- flu_features %>%
   filter(week >= eval_start,
          week <= eval_end + weeks(4),
          !is.na(.data[[target_var]])) %>%
   dplyr::select(week, actual = all_of(target_var))

selected_origins <- retro_df %>%
   filter(origin_week >= eval_start,
          origin_week <= eval_end) %>%
   distinct(origin_week) %>%
   arrange(origin_week) %>%
   dplyr::slice(1:16) %>%
   pull(origin_week)

origin_path_df <- retro_df %>%
   filter(origin_week %in% selected_origins) %>%
   mutate(
      origin_label = factor(
         paste0("Issued ", format(origin_week, "%b %d")),
         levels = paste0("Issued ", format(sort(selected_origins), "%b %d"))
      ),
      horizon_label = paste0("H", horizon, " / ", horizon + 1, "-wk ahead")
   ) %>%
   mutate(panel_xmin = origin_week - weeks(3),
          panel_xmax = target_end_date)

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
      title    = paste0("Forecast Paths by Origin Week - ", MODEL_LABEL),
      subtitle = "Each blue line is one issued 4-week forecast path; black is observed CDC",
      x = "Week Ending Date",
      y = "Weekly Inpatient Hospitalizations (Influenza, SC)",
      shape   = "Forecast horizon",
      caption = "Each step conditions on the prior step's prediction only."
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
