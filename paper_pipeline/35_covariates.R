# ==============================================================================
# 35_covariates.R
# Tract-year covariates, from the same block files the segregation index uses,
# plus ACS per-capita income.
#
# Produces:
#   clean/p2_tract_lodes_cov_panel_v3.rds
#     tract_id, year,
#     workers_rac, workers_wac,
#     n_white_rac, n_black_rac, pct_white_rac, pct_black_rac,
#     pct_lowincome_rac, pct_highincome_rac,
#     pct_white_wac, pct_black_wac, pct_lowincome_wac, pct_manuf_wac
#   clean/p2_tract_income_panel.rds
#     tract_id, year, income_percapita
#
# THE v3 SUFFIX IS LOAD-BEARING. 60_co_setup.R refuses to run against the
# un-versioned file because that one built "White" as the non-Black
# complement. Here White is CR01 (White alone) and Black is CR02 (Black
# alone), and the raw counts are carried alongside the shares so the
# correction is auditable downstream. 63 asserts pct_white_rac,
# pct_white_wac, n_white_rac and pct_black_rac are all present.
#
# LOW-WAGE = CE01, jobs paying <= $1,250/month. That is the definition the
# manuscript states; do not swap it for CE01+CE02 without changing the text.
#
# INCOME: ACS 5-year B19301_001 (per-capita income in survey-year dollars),
# one release per panel year. Needs a Census API key --
#   tidycensus::census_api_key("YOUR KEY", install = TRUE)
# See the tract-vintage note at the income section: releases on 2010-vintage
# tracts are bridged to the 2020 tracts LODES 8 uses.
# ==============================================================================

source("30_p2_setup.R")

pct_of <- function(num, den) ifelse(den > 0, 100 * num / den, NA_real_)

## =============================================================================
## 1. LODES covariates
## =============================================================================
cov_file <- file.path(DIR_CLEAN, "p2_tract_lodes_cov_panel_v3.rds")

if (file.exists(cov_file)) {
  message("  p2_tract_lodes_cov_panel_v3.rds exists")
} else {
  tract_side <- function(yr, side) {
    f <- file.path(DIR_CLEAN, sprintf("blocks_%s_%s.rds", side, yr))
    if (!file.exists(f)) stop("Missing ", f, " -- run 20_lodes_blocks.R.")
    b <- readRDS(f)
    keep <- intersect(c("C000", "CR01", "CR02", "CE01", "CE03", "CNS05"),
                      names(b))
    b |> group_by(tract_id) |>
      summarise(across(all_of(keep), ~ sum(.x, na.rm = TRUE)), .groups = "drop")
  }

  lodes_cov <- map(P2_YEARS, function(yr) {
    message("Covariates ", yr, " ...")
    r <- tract_side(yr, "rac")
    w <- tract_side(yr, "wac")
    rr <- r |> transmute(
      tract_id,
      workers_rac       = C000,
      n_white_rac       = CR01,                    # White ALONE (the v3 fix)
      n_black_rac       = CR02,
      pct_white_rac     = pct_of(CR01, C000),
      pct_black_rac     = pct_of(CR02, C000),
      pct_lowincome_rac = pct_of(CE01, C000),      # <= $1,250/month
      pct_highincome_rac = pct_of(CE03, C000))     # >  $3,333/month
    ww <- w |> transmute(
      tract_id,
      workers_wac       = C000,
      pct_white_wac     = pct_of(CR01, C000),
      pct_black_wac     = pct_of(CR02, C000),
      pct_lowincome_wac = pct_of(CE01, C000),
      pct_manuf_wac     = if ("CNS05" %in% names(w)) pct_of(CNS05, C000)
                          else NA_real_)
    full_join(rr, ww, by = "tract_id") |>
      mutate(year = as.integer(yr)) |>
      relocate(tract_id, year)
  }) |> bind_rows() |> arrange(tract_id, year)

  # the guard 63 mirrors
  stopifnot(all(c("pct_white_rac", "pct_white_wac", "n_white_rac",
                  "pct_black_rac") %in% names(lodes_cov)))
  saveRDS(lodes_cov, cov_file)
  message("  wrote p2_tract_lodes_cov_panel_v3.rds (", nrow(lodes_cov),
          " tract-years)")
  write_diag(
    lodes_cov |> group_by(year) |>
      summarise(n = n(),
                mean_pct_white_rac = round(mean(pct_white_rac, na.rm = TRUE), 2),
                mean_pct_black_rac = round(mean(pct_black_rac, na.rm = TRUE), 2),
                mean_pct_low_rac   = round(mean(pct_lowincome_rac, na.rm = TRUE), 2)),
    "35_lodes_covariates_by_year")
}

## =============================================================================
## 2. ACS per-capita income
## =============================================================================
# TRACT VINTAGE. LODES 8 is on 2020 census tracts for every year. ACS 5-year
# releases before the 2020 release are on 2010 tracts. Those releases are
# mapped onto 2020 tracts with the Census 2010->2020 tract relationship file
# by DOMINANT PARENT: each 2020 tract takes the value of the single 2010 tract
# that contributes the most of its land area. Per-capita income is an intensive
# quantity, so a tract whose boundary did not change -- the large majority --
# passes through exactly, and a split tract inherits its parent's value. The
# vintage of each release is decided from the data (share of ids that are 2020
# tract ids), not assumed from the release year.
inc_file <- file.path(DIR_CLEAN, "p2_tract_income_panel.rds")

if (file.exists(inc_file)) {
  message("  p2_tract_income_panel.rds exists")
} else {
  suppressPackageStartupMessages(library(tidycensus))
  if (Sys.getenv("CENSUS_API_KEY") == "")
    stop("No Census API key. Run tidycensus::census_api_key(\"<key>\", ",
         "install = TRUE) once, restart R, then re-run.\n",
         "Free key: https://api.census.gov/data/key_signup.html")

  ## 2010 -> 2020 tract crosswalk (land area of the overlap)
  xw_cache <- file.path(DIR_RAW, "tract_rel_2010_2020_co.rds")
  xw <- if (file.exists(xw_cache)) readRDS(xw_cache) else {
    message("Downloading the 2010->2020 tract relationship file ...")
    url <- paste0("https://www2.census.gov/geo/docs/maps-data/data/rel2020/",
                  "tract/tab20_tract20_tract10_st08.txt")
    d <- readr::read_delim(url, delim = "|", show_col_types = FALSE)
    d <- d |>
      transmute(tract20 = as.character(GEOID_TRACT_20),
                tract10 = as.character(GEOID_TRACT_10),
                w = as.numeric(AREALAND_PART)) |>
      filter(!is.na(tract20), !is.na(tract10), w > 0)
    saveRDS(d, xw_cache); d
  }

  # the 2020-vintage tract ids LODES8 uses, for the vintage test below
  ids_2020 <- readRDS(file.path(DIR_CLEAN, "tract_centroids_km.rds"))$tract_id

  # dominant parent: the 2010 tract contributing the most land to each 2020
  # tract; the 2020 tract stays missing when that parent is missing
  dominant_parent <- xw |>
    group_by(tract20) |>
    slice_max(w, n = 1, with_ties = FALSE) |>
    ungroup()

  to_2020 <- function(df) {              # df: tract_id (2010), income_percapita
    dominant_parent |>
      left_join(df, by = c("tract10" = "tract_id")) |>
      transmute(tract_id = tract20, income_percapita) |>
      filter(!is.na(income_percapita))
  }

  # RELEASE. Panel year Y uses the ACS 5-year release ending in Y + 1, capped at
  # the latest release that exists (ACS_LATEST_RELEASE, default 2024). The cap
  # is explicit so that the vintage never depends on when the script is run or
  # on a transient download failure: any other error stops the script.
  ACS_OFFSET <- 1L
  ACS_LATEST <- as.integer(Sys.getenv("ACS_LATEST_RELEASE", "2024"))

  # one ACS release, cached by the RELEASE year
  get_release <- function(acs_yr) {
    ck <- file.path(DIR_RAW, sprintf("acs_b19301_co_%s.rds", acs_yr))
    if (file.exists(ck)) return(readRDS(ck))
    message("ACS B19301, ", acs_yr, " release ...")
    x <- suppressMessages(
      tidycensus::get_acs(geography = "tract", variables = "B19301_001",
                          state = "08", year = acs_yr, survey = "acs5"))
    x <- x |> transmute(tract_id = as.character(GEOID),
                        income_percapita = estimate)
    saveRDS(x, ck)
    x
  }

  income <- map(P2_YEARS, function(yr) {
    acs_yr <- min(yr + ACS_OFFSET, ACS_LATEST)
    d <- get_release(acs_yr)

    # Decide the tract vintage from the data, not from the year: the 2019
    # release matches 85% of 2020 tract ids, the 2020 release 100%.
    share20 <- mean(d$tract_id %in% ids_2020)
    if (share20 < 0.95) {
      message("  panel ", yr, " <- ", acs_yr, " release: ",
              round(100 * share20), "% 2020 ids -- crosswalking from 2010")
      d <- to_2020(d)
    } else {
      message("  panel ", yr, " <- ", acs_yr, " release")
    }
    d |> mutate(year = as.integer(yr)) |> relocate(tract_id, year)
  }) |> bind_rows() |> arrange(tract_id, year)

  saveRDS(income, inc_file)
  message("  wrote p2_tract_income_panel.rds (", nrow(income), " tract-years)")
  write_diag(
    income |> group_by(year) |>
      summarise(n = n(),
                pctna = round(100 * mean(is.na(income_percapita)), 2),
                median_income = round(median(income_percapita, na.rm = TRUE))),
    "35_income_by_year")
}

message("35_covariates.R complete.")
