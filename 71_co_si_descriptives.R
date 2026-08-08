# ==============================================================================
# 71_co_si_descriptives.R      [PAPER 4, step 11]
# SI Tables S1/S2: WHERE workers in each LODES flow group LIVE and WORK in
# the Denver MSA, 2023.
#
# Groups: earnings brackets SE01 (<=$1,250/mo), SE02, SE03 (>$3,333/mo) and
# industry SUPERGROUPS SI01 (goods), SI02 (trade/transport/utilities),
# SI03 (services). NOTE FOR THE TEXT: LODES carries NO occupation data --
# industry supergroups are the closest available construct (CTPP is the flow
# source for occupation). Residence and workplace distributions are taken
# from the OD file's own group columns (summed by home tract and by work
# tract respectively), which keeps them exactly consistent with the flow
# analyses and avoids hand-mapping CNS sectors onto supergroups.
#
# Frame: Denver MSA (19740) tract pairs, same-MSA restriction, 2023 -- the
# manuscript's frame. Per group x side:
#   * share of workers by county (Adams/Arapahoe/Broomfield/Denver/Douglas/
#     Jefferson; fringe counties pooled as Other)
#   * share living in each exclusionary-zoning tercile (residence side; the
#     terciles of 69, defined on the in-scope model frame; workers in tracts
#     outside the zoning sample are reported as such, never silently dropped)
#   * share working in a major employment center (workplace side; top 2% of
#     regional tracts by 2023 jobs, 68's definition)
#   * mean distance of residence / workplace from the CBD; mean commute
#
# Needs: 62's OD caches, 63's panel (terciles), 66's WAC margin cache,
#        68's accessibility file.
# Output: output/models/p4_si_tabS1_residence.csv
#         output/models/p4_si_tabS2_workplace.csv
# ==============================================================================

source("60_co_setup.R")

## ---- pieces ------------------------------------------------------------------
cent <- readRDS(file.path(DIR_CLEAN, "tract_centroids_km.rds")) |>
  transmute(tract_id = as.character(GEOID), CBSA_Code)
acc  <- readRDS(file.path(DIR_CO_CLEAN, "co_accessibility_2023.rds")) |>
  select(tract_id, dist_cbd_km)
dat  <- readRDS(file.path(DIR_CO_OUT, "co_analysis_home_panel.rds"))

# exclusionary terciles on the model frame (exactly 69's construction)
terc <- dat |>
  filter(year == CO_ANCHOR_YEAR, in_scope, denver_msa,
         n_commuters >= P3_MIN_COMMUTERS) |>
  transmute(tract_id, tercile = ntile(pct_reslow_of_res, 3))

# employment centers: top 2% of REGIONAL tracts by 2023 WAC jobs (as in 68)
wac <- readRDS(file.path(DIR_CO_CLEAN,
                         sprintf("co_wac_tract_%d.rds", CO_ANCHOR_YEAR))) |>
  transmute(tract_id = as.character(tract_id), jobs = C000) |>
  inner_join(cent |> filter(CBSA_Code %in% CO_CBSA_KEEP), by = "tract_id")
centers <- wac$tract_id[wac$jobs >= quantile(wac$jobs, 0.98)]
message(length(centers), " employment-center tracts")

# 2023 OD, Denver MSA, same-MSA restriction (62's caches, 62's rule)
od <- map(P3_OD_PARTS, function(part) {
  for (f in c(p3_od_cache(part, CO_ANCHOR_YEAR, "co"),
              co_od_cache(part, CO_ANCHOR_YEAR)))
    if (file.exists(f)) return(readRDS(f))
  NULL
}) |> compact() |> bind_rows() |>
  group_by(w_tract, h_tract) |>
  summarise(across(any_of(P3_OD_COLS), ~ sum(.x, na.rm = TRUE)),
            .groups = "drop") |>
  inner_join(cent, by = c("h_tract" = "tract_id")) |>
  rename(cbsa_h = CBSA_Code) |>
  inner_join(cent, by = c("w_tract" = "tract_id")) |>
  rename(cbsa_w = CBSA_Code) |>
  filter(cbsa_h == CO_CBSA_MAIN, cbsa_w == CO_CBSA_MAIN)

# commute distance (same construction as 62)
xy <- readRDS(file.path(DIR_CLEAN, "tract_centroids_km.rds")) |>
  transmute(tract_id = as.character(GEOID), X_km, Y_km)
od <- od |>
  left_join(xy, by = c("h_tract" = "tract_id")) |>
  rename(xh = X_km, yh = Y_km) |>
  left_join(xy, by = c("w_tract" = "tract_id")) |>
  rename(xw = X_km, yw = Y_km) |>
  mutate(dist_km = sqrt((xh - xw)^2 + (yh - yw)^2))

COUNTY_NAMES <- c(`08001` = "Adams", `08005` = "Arapahoe",
                  `08014` = "Broomfield", `08031` = "Denver",
                  `08035` = "Douglas", `08059` = "Jefferson")
county_of <- function(tract) {
  cf <- substr(tract, 1, 5)
  ifelse(cf %in% names(COUNTY_NAMES), COUNTY_NAMES[cf], "Other")
}

GROUPS <- c(`All workers` = "S000",
  `Earnings <= $1,250/mo` = "SE01", `Earnings $1,251-3,333/mo` = "SE02",
  `Earnings > $3,333/mo` = "SE03",
  `Goods-producing` = "SI01", `Trade/transport/utilities` = "SI02",
  `All other services` = "SI03")

side_table <- function(side_tract) {
  imap(GROUPS, function(col, gname) {
    w <- od[[col]]
    df <- tibble(tract = od[[side_tract]], w = w, dist_km = od$dist_km) |>
      filter(w > 0)
    tot <- sum(df$w)
    shares <- df |>
      mutate(county = county_of(tract)) |>
      group_by(county) |>
      summarise(pct = 100 * sum(w) / tot, .groups = "drop") |>
      pivot_wider(names_from = county, values_from = pct)
    extra <- if (side_tract == "h_tract") {
      tl <- df |> left_join(terc, by = c("tract" = "tract_id"))
      tibble(
        pct_terc_least = 100 * sum(tl$w[tl$tercile == 1], na.rm = TRUE) / tot,
        pct_terc_mid   = 100 * sum(tl$w[tl$tercile == 2], na.rm = TRUE) / tot,
        pct_terc_most  = 100 * sum(tl$w[tl$tercile == 3], na.rm = TRUE) / tot,
        pct_outside_zoning_sample =
          100 * sum(tl$w[is.na(tl$tercile)]) / tot)
    } else {
      tibble(pct_in_employment_center =
               100 * sum(df$w[df$tract %in% centers]) / tot)
    }
    dcbd <- df |> left_join(acc, by = c("tract" = "tract_id"))
    bind_cols(tibble(group = gname, workers = tot), shares, extra,
              tibble(mean_dist_cbd_km =
                       weighted.mean(dcbd$dist_cbd_km, dcbd$w, na.rm = TRUE),
                     mean_commute_km = weighted.mean(df$dist_km, df$w)))
  }) |> bind_rows() |>
    mutate(across(where(is.numeric), ~ round(.x, 1)))
}

s1 <- side_table("h_tract")
s2 <- side_table("w_tract")
write.csv(s1, file.path(DIR_CO_MOD, "p4_si_tabS1_residence.csv"),
          row.names = FALSE)
write.csv(s2, file.path(DIR_CO_MOD, "p4_si_tabS2_workplace.csv"),
          row.names = FALSE)
print(as.data.frame(s1)); print(as.data.frame(s2))

## ---- DIAGNOSTICS -------------------------------------------------------------
# partition check: SE and SI families must each sum to ~ the S000 total
fam <- s1 |> mutate(fam = case_when(grepl("Earnings", group) ~ "SE",
                                    group == "All workers" ~ "all",
                                    TRUE ~ "SI")) |>
  group_by(fam) |> summarise(workers = sum(workers), .groups = "drop")
write_codiag(fam, "71_partition_check")
message("71 complete. SI tables in ", DIR_CO_MOD)
