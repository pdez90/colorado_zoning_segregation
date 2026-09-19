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
# See the tract-vintage note at the income section: 2011-2020 ACS releases are
# on 2010-vintage tracts and are crosswalked to the 2020 tracts LODES8 uses.
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
# TRACT VINTAGE. LODES8 is on 2020 census tracts for every year. ACS 5-year
# releases through 2020 are published on 2010 tracts; 2021+ are on 2020
# tracts. Pre-2021 years are therefore mapped onto 2020 tracts with the
# Census 2010->2020 tract relationship file, taking the land-area-weighted
# mean of the 2010 tracts overlapping each 2020 tract. For a tract whose
# boundary did not change -- the large majority -- the mapping is one-to-one
# and the value passes through unchanged. Splits and merges are approximated.
# This bridge is an approximation for the minority of tracts that were split
# or merged in 2020; 90_validate_rebuild.R reports its effect on the estimates.
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

  # DOMINANT-PARENT assignment: each 2020 tract takes the value of the single
  # 2010 tract contributing the most of its land, and stays missing when that
  # parent is missing. Chosen on evidence -- 39_income_crosswalk_test.R
  # compares missingness against the published panel across the eight bridged
  # years (2011-2018; panel 2019 reads the 2020 release and needs no bridge):
  #
  #   dominant parent              total gap 1.3   (1.0 / 1.0 / 0.8 ...)
  #   area-weighted mean           total gap 4.9   (0.4 every year)
  #   area-weighted, all parents   total gap 5.7   (1.8 / 1.8 / 1.7 ...)
  #
  # against a reference of 1.1 / 1.0 / 1.0 ... It is also the same rule the
  # SLD turned out to use (see 53_sld.R), which is what a pipeline with one
  # crosswalk helper would do.
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

  # RELEASE OFFSET. Panel year Y reads the ACS 5-year release ending Y+1, not
  # Y. This is not a guess: 36_income_variant_test.R rebuilt the published SI
  # tercile means (56.286 / 57.023 / 66.362) from eight candidate
  # constructions, and B19301 from the *2024* release reproduced all three to
  # 0.000 for panel year 2023, while every 2023-release construction came in
  # $2.0k-$2.9k low and the 2022 release lower still. The original pipeline
  # was built in August 2026, by which point the 2024 release (December 2025)
  # was the newest available. If a release is not published yet, fall back to
  # the panel year itself.
  ACS_OFFSET <- 1L

  # one ACS release, cached by the RELEASE year; NULL if it is not published
  get_release <- function(acs_yr) {
    ck <- file.path(DIR_RAW, sprintf("acs_b19301_co_%s.rds", acs_yr))
    if (file.exists(ck)) return(readRDS(ck))
    message("ACS B19301, ", acs_yr, " release ...")
    x <- tryCatch(
      suppressMessages(
        tidycensus::get_acs(geography = "tract", variables = "B19301_001",
                            state = "08", year = acs_yr, survey = "acs5")),
      error = function(e) NULL)
    if (is.null(x)) return(NULL)
    x <- x |> transmute(tract_id = as.character(GEOID),
                        income_percapita = estimate)
    saveRDS(x, ck)
    x
  }

  income <- map(P2_YEARS, function(yr) {
    acs_yr <- yr + ACS_OFFSET
    d <- get_release(acs_yr)
    if (is.null(d)) {                       # release not out yet
      message("  ", acs_yr, " release unavailable -- using the ", yr, " one")
      acs_yr <- yr
      d <- get_release(acs_yr)
    }
    if (is.null(d)) stop("No ACS release available for panel year ", yr)

    # Decide the tract vintage from the data, not from the year. The switch to
    # 2020 tracts happened with the 2020 ACS 5-year release, NOT 2021 as the
    # documentation order suggests -- checked against tract_centroids_km.rds:
    # the 2019 release matches 85% of 2020 tract ids, the 2020 release 100%.
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
