# ==============================================================================
# 66_co_group_flows.R      [COLORADO ZONING case study, step 6]
# WHO travels how far, into WHAT kind of workhood -- group-specific commute
# and mixing measures from the LODES groups that have REAL flows.
#
# DATA FACTS driving this design (LODES 8.4, verified):
#   * The OD file's ONLY group flows are: earnings SE01 (<=$1250/mo), SE02,
#     SE03 (>$3333/mo); industry supergroup SI01 (goods), SI02 (trade/
#     transport/utilities), SI03 (all other services); age SA01 (<=29),
#     SA02 (30-54), SA03 (55+). NO occupation (not in LODES at all -- CTPP
#     is the flow source for occupation), NO education, NO race in OD.
#   * Education (CD01-04, workers 30+), detailed industry (CNS01-20), and
#     full earnings (CE01-03) DO exist as RAC/WAC MARGINS. They enter here
#     as attributes of DESTINATION tracts (what kind of workforce mix a
#     commuter lands in), which needs no flow assumption, and as HOME-tract
#     population weights for descriptives (shared-destination assumption,
#     same as Paper 3's race handling -- labeled wherever used).
#
# Measures, per HOME tract x year x group (group in S000/SE/SI/SA):
#   * n_commuters, mean/share distance measures, eff_n_dest  (real flows)
#   * wexp (flow-weighted destination workplace D, whiteblack half)
#   * destination-mix exposure: flow-weighted WORKFORCE DIVERSITY of the
#     destination tracts -- normalized entropy of the WAC margins:
#       - earnings mix   (CE01-03)
#       - education mix  (CD01-04)
#       - industry mix   (CNS01-20)
#     "does this tract's workforce commute into class-diverse workhoods or
#      class-monoculture ones?" -- the mixing question, assumption-free.
#
# Needs: 62's cached OD parts (run 62 first), Paper 2 seg panel + centroids.
# New downloads: Colorado WAC + RAC, tract-agg, per year (~5 MB/yr, cached).
#
# NATIONAL PATHWAY (noted for the broader paper): 51's national OD download
# already carries SE/SI/SA columns, so the aggregation below lifts into 52
# unchanged; the WAC/RAC margin pulls generalize by looping states (or by
# widening KEEP_BLOCK in 31 and re-aggregating blocks).
#
# Output: clean/co_wac_diversity_panel.rds   (destination-mix attributes)
#         clean/co_group_flows_panel.rds     (tract x year x group, long)
# ==============================================================================

source("60_co_setup.R")
library(lehdr)

seg_panel <- readRDS(file.path(DIR_P2_OUT, "p2_tract_segregation_panel.rds"))
tract_centroids <- readRDS(file.path(DIR_CLEAN, "tract_centroids_km.rds"))
cent <- tract_centroids |>
  filter(CBSA_Code %in% CO_CBSA_KEEP) |>
  select(GEOID, CBSA_Code, X_km, Y_km)

## ---- group definitions -------------------------------------------------------
CO_FLOW_GROUPS <- c(all = "S000",
                    earn_low = "SE01", earn_mid = "SE02", earn_high = "SE03",
                    ind_goods = "SI01", ind_tradetrans = "SI02",
                    ind_services = "SI03",
                    age_u30 = "SA01", age_30_54 = "SA02", age_55p = "SA03")

CE_COLS  <- sprintf("CE%02d", 1:3)
CD_COLS  <- sprintf("CD%02d", 1:4)
CNS_COLS <- sprintf("CNS%02d", 1:20)

# normalized Shannon entropy of a count vector (0 = monoculture, 1 = even)
norm_entropy <- function(m) {           # m: matrix, rows = tracts
  tot <- rowSums(m)
  p <- m / pmax(tot, .Machine$double.eps)
  h <- -rowSums(ifelse(p > 0, p * log(p), 0))
  ifelse(tot > 0, h / log(ncol(m)), NA_real_)
}

## =============================================================================
## A. destination-tract workforce-mix attributes (WAC margins) + home weights
## =============================================================================
grab_margin <- function(yr, type) {     # type: "wac" or "rac"
  cache <- file.path(DIR_CO_CLEAN, sprintf("co_%s_tract_%s.rds", type, yr))
  if (file.exists(cache)) return(readRDS(cache))
  message(sprintf("Downloading %s tract co %s", toupper(type), yr))
  df <- grab_lodes(state = "co", year = yr, version = "LODES8",
                   lodes_type = type, job_type = "JT01", agg_geo = "tract")
  id_col <- if (type == "wac") "w_tract" else "h_tract"
  df <- df |>
    mutate(tract_id = as.character(.data[[id_col]])) |>
    select(tract_id, any_of(c("C000", CE_COLS, CD_COLS, CNS_COLS)))
  saveRDS(df, cache)
  df
}

wac_div_file <- file.path(DIR_CO_CLEAN, "co_wac_diversity_panel.rds")
if (file.exists(wac_div_file)) {
  message(basename(wac_div_file), " exists")
  wac_div <- readRDS(wac_div_file)
} else {
  wac_div <- map(CO_YEARS, function(yr) {
    w <- grab_margin(yr, "wac")
    need <- c(CE_COLS, CD_COLS, CNS_COLS)
    miss <- setdiff(need, names(w))
    if (length(miss) > 0) stop("WAC ", yr, " missing columns: ",
                               paste(miss, collapse = ", "))
    tibble(
      tract_id = w$tract_id, year = yr,
      workers_wac_margin = w$C000,
      earn_entropy_wac = norm_entropy(as.matrix(w[CE_COLS])),
      edu_entropy_wac  = norm_entropy(as.matrix(w[CD_COLS])),
      ind_entropy_wac  = norm_entropy(as.matrix(w[CNS_COLS])))
  }) |> bind_rows()
  saveRDS(wac_div, wac_div_file)
  message("Saved WAC diversity panel: ", nrow(wac_div), " tract-years")
}

# home-tract education/earnings composition (RAC margins) -- used ONLY as
# population weights in 67's descriptive rows (shared-destination
# assumption, education has no OD flows; labeled there)
rac_wt_file <- file.path(DIR_CO_CLEAN, "co_rac_weights_panel.rds")
if (file.exists(rac_wt_file)) {
  message(basename(rac_wt_file), " exists")
} else {
  rac_wt <- map(CO_YEARS, function(yr) {
    r <- grab_margin(yr, "rac")
    r |> select(tract_id, any_of(c("C000", CE_COLS, CD_COLS))) |>
      mutate(year = yr)
  }) |> bind_rows()
  saveRDS(rac_wt, rac_wt_file)
  message("Saved RAC weights panel: ", nrow(rac_wt), " tract-years")
}

## =============================================================================
## B. per-group OD aggregation (real flows)
## =============================================================================
out_file <- file.path(DIR_CO_CLEAN, "co_group_flows_panel.rds")
if (file.exists(out_file)) {
  message(basename(out_file), " exists -- delete to recompute.")
} else {

read_co_od <- function(yr) {            # 62's caches (either location)
  per_part <- map(P3_OD_PARTS, function(part) {
    for (f in c(p3_od_cache(part, yr, "co"), co_od_cache(part, yr)))
      if (file.exists(f)) return(readRDS(f))
    NULL
  }) |> compact()
  if (length(per_part) == 0)
    stop("No cached OD parts for ", yr, " -- run 62 first.")
  bind_rows(per_part) |>
    group_by(w_tract, h_tract) |>
    summarise(across(any_of(P3_OD_COLS), ~ sum(.x, na.rm = TRUE)),
              .groups = "drop")
}

wmean_na <- function(d, w) {
  ok <- !is.na(d) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(d[ok] * w[ok]) / sum(w[ok])
}

group_panel <- map(CO_YEARS, function(yr) {
  od <- read_co_od(yr) |>
    inner_join(cent, by = c("h_tract" = "GEOID")) |>
    rename(cbsa_h = CBSA_Code, xh = X_km, yh = Y_km) |>
    inner_join(cent, by = c("w_tract" = "GEOID")) |>
    rename(cbsa_w = CBSA_Code, xw = X_km, yw = Y_km) |>
    filter(cbsa_h == cbsa_w) |>                    # same-MSA rule, as 52/62
    mutate(dist_km = sqrt((xh - xw)^2 + (yh - yw)^2))

  seg_yr <- seg_panel |> filter(year == yr) |>
    select(tract_id, d_whiteblack_wac_half)
  div_yr <- wac_div |> filter(year == yr) |>
    select(tract_id, earn_entropy_wac, edu_entropy_wac, ind_entropy_wac)
  od <- od |>
    left_join(seg_yr, by = c("w_tract" = "tract_id")) |>
    left_join(div_yr, by = c("w_tract" = "tract_id"))

  per_group <- imap(CO_FLOW_GROUPS, function(col, gname) {
    if (!col %in% names(od)) return(NULL)
    w <- od[[col]]
    od |>
      mutate(.w = w) |>
      filter(.w > 0) |>
      group_by(tract_id = h_tract) |>
      summarise(
        n_commuters = sum(.w),
        mean_dist_km = weighted.mean(dist_km, .w),
        pct_lt5km  = 100 * sum(.w[dist_km < 5])  / sum(.w),
        pct_gt24km = 100 * sum(.w[dist_km > 24]) / sum(.w),
        eff_n_dest = 1 / sum((.w / sum(.w))^2),
        wexp_whiteblack = wmean_na(d_whiteblack_wac_half, .w),
        dest_earn_entropy = wmean_na(earn_entropy_wac, .w),
        dest_edu_entropy  = wmean_na(edu_entropy_wac, .w),
        dest_ind_entropy  = wmean_na(ind_entropy_wac, .w),
        .groups = "drop") |>
      mutate(group = gname, year = yr)
  }) |> compact() |> bind_rows()

  message(yr, ": ", n_distinct(per_group$group), " groups aggregated")
  rm(od); gc()
  per_group
}) |> bind_rows()

saveRDS(group_panel, out_file)
message("Saved group flows panel: ", nrow(group_panel), " tract-year-groups")

## ---- DIAGNOSTICS -------------------------------------------------------------
# (a) group flows must (nearly) partition totals: SE and SI and SA each sum
# to S000 in the raw file; deviations = suppression noise, should be tiny
chk <- group_panel |>
  filter(year == CO_ANCHOR_YEAR) |>
  group_by(group) |>
  summarise(workers = sum(n_commuters), .groups = "drop")
tot <- chk$workers[chk$group == "all"]
write_codiag(chk |>
  mutate(family = case_when(grepl("^earn", group) ~ "earnings",
                            grepl("^ind",  group) ~ "industry",
                            grepl("^age",  group) ~ "age",
                            TRUE ~ "all"),
         pct_of_total = round(100 * workers / tot, 2)),
  "66_group_partition_check")
# Eyeball: earnings/industry/age families should each sum to ~100%.

# (b) the headline gradient: distance and destination-mix by group by year
write_codiag(
  group_panel |>
    group_by(year, group) |>
    summarise(
      workers = sum(n_commuters),
      mean_dist  = weighted.mean(mean_dist_km, n_commuters, na.rm = TRUE),
      pct_gt24km = weighted.mean(pct_gt24km, n_commuters, na.rm = TRUE),
      eff_n_dest = weighted.mean(eff_n_dest, n_commuters, na.rm = TRUE),
      wexp       = weighted.mean(wexp_whiteblack, n_commuters, na.rm = TRUE),
      dest_earn_entropy = weighted.mean(dest_earn_entropy, n_commuters,
                                        na.rm = TRUE),
      .groups = "drop"),
  "66_group_gradients_by_year")
}
message("66 complete.")
