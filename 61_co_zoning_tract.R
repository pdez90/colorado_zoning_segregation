# ==============================================================================
# 61_co_zoning_tract.R      [COLORADO ZONING case study, step 1]
# Zoning polygons -> 2020-tract-level zoning measures.
#
# For every Colorado tract touched by the zoning file:
#   * cover_zoned      : UNIONED zoned area / tract POLYGON area (land+water).
#                        Denominator is the polygon, not ALAND: jurisdictions
#                        zone their reservoirs (Standley Lake sits in a tract
#                        with 5 km2 of water), so ALAND denominators push
#                        coverage far above 1. Uncovered share = land whose
#                        jurisdiction is NOT in the file = UNKNOWN zoning,
#                        never zero.
#   * pct_<group>      : share of the tract's zoned area in each analysis
#                        group. Areas are DISSOLVED (st_union) per tract x
#                        group before measuring: the layer contains stacked/
#                        duplicate district polygons, and raw piece sums
#                        double-count them (worst tract: 2.5x its area).
#                        Cross-group double claims (county and town zoning
#                        the same ground differently, e.g. Foxfield/Arapahoe)
#                        survive the dissolve; shares are normalized by the
#                        group-area total, i.e. split proportionally, and
#                        tracts with >10% cross-group excess are logged.
#   * pct_reslow_of_res: Residential_Low share of RESIDENTIAL land -- the
#                        exclusionary-zoning headline (parallels "% of
#                        residential land zoned single-family-only" in the
#                        zoning-atlas literature)
#   * pct_job_zone     : Commercial + Industrial + mixed_res share
#   * zoning_entropy   : normalized Shannon entropy of the 7 non-uncertain
#                        groups (0 = monoculture, 1 = even mix)
#   * pct_adu_res      : share of residential land where ADUs are allowed
#   * pct_grouphome_res: share of residential land where group homes allowed
#   * n_jurisd, jurisd_main : jurisdictions with >=5% of the tract's zoned
#                        area; the largest = the model clustering unit
#                        (zoning is WRITTEN at the jurisdiction level)
#
# Geometry: TIGER 2024 CO tracts (= 2020-vintage boundaries/GEOIDs, matching
# LODES 2020 blocks); everything projected to CRS_METERS (EPSG:5070) before
# areas are taken, as in the rest of the pipeline.
#
# VERIFIED against the raw layer (Aug 2026): 100,837 polygons, 9 invalid
# (fixed by st_make_valid); 818 CO tracts intersect; ~730 pass the 80%
# coverage screen (Denver MSA ~660, Boulder ~70, Weld 3); fringe counties
# (Elbert/Park/Clear Creek/Gilpin/El Paso/Larimer edges) drop out at ~0%.
#
# Output: clean/co_tract_zoning.rds       (measures, one row per tract)
#         clean/co_tract_geom.rds         (sf: tract geometry for figures)
#         clean/co_zoning_jurisd.rds      (jurisdiction-level summary)
# ==============================================================================

source("60_co_setup.R")

out_file <- file.path(DIR_CO_CLEAN, "co_tract_zoning.rds")
if (file.exists(out_file)) {
  message(basename(out_file), " exists -- delete it to recompute.")
} else {

## ---- 1. read + clean the zoning layer ---------------------------------------
if (!file.exists(CO_ZONING_SHP)) stop("Zoning shapefile not found: ",
                                      CO_ZONING_SHP)
message("Reading zoning shapefile (100k polygons, ~1 min) ...")
zon <- st_read(CO_ZONING_SHP, quiet = TRUE) |>
  st_zm(drop = TRUE) |>
  st_transform(CRS_METERS)
n_polygons_raw <- nrow(zon)
zon <- st_make_valid(zon)          # 9 invalid polygons in the Oct 2023 layer
zon <- zon[!st_is_empty(zon), ]
write_codiag(tibble(polygons_in_layer = n_polygons_raw,
                    polygons_with_valid_geometry = nrow(zon),
                    jurisdictions = dplyr::n_distinct(zon$jurisd)),
             "61_zoning_layer_counts")

# normalize the messy free-text Yes/No flags ("Yes'", "YEs", "N", "" ...)
norm_yn <- function(x) {
  x <- gsub("[^a-z]", "", tolower(trimws(x)))
  dplyr::case_when(x %in% c("yes", "y") ~ "yes",
                   x %in% c("no", "n")  ~ "no",
                   TRUE                 ~ NA_character_)
}

# GenZone2 -> analysis group; STOP on any value not in the verified list
# (spelling drift would silently misclassify land)
zone_lookup <- tibble(
  GenZone2 = unlist(CO_ZONE_GROUPS, use.names = FALSE),
  zone_group = rep(names(CO_ZONE_GROUPS), lengths(CO_ZONE_GROUPS)))

zon <- zon |>
  mutate(GenZone2 = trimws(GenZone2),
         jurisd   = trimws(jurisd),
         adu_yn   = norm_yn(ADU),
         gh_yn    = norm_yn(GroupHome))
unknown <- setdiff(unique(zon$GenZone2), zone_lookup$GenZone2)
if (length(unknown) > 0)
  stop("GenZone2 values not in CO_ZONE_GROUPS (fix 60 before proceeding): ",
       paste(unknown, collapse = ", "))
zon <- zon |> left_join(zone_lookup, by = "GenZone2")

write_codiag(zon |> st_drop_geometry() |> count(zone_group, GenZone2),
             "61_zoning_class_counts")
write_codiag(zon |> st_drop_geometry() |>
               count(jurisd, sort = TRUE), "61_zoning_jurisdiction_counts")

## ---- 2. tract geometry (2020 vintage) ---------------------------------------
if (!dir.exists(CO_TIGER_DIR))
  stop("TIGER 2024 tract shapefile not found at ", CO_TIGER_DIR,
       " -- run paper_pipeline/10_geography.R, which writes it.")
tr_shp <- st_read(CO_TIGER_DIR, quiet = TRUE)
tracts <- tr_shp |>
  st_transform(CRS_METERS) |>
  transmute(tract_id  = as.character(GEOID),
            aland_km2 = as.numeric(ALAND) / 1e6)
tracts$poly_km2 <- as.numeric(st_area(tracts)) / 1e6  # land + water

# keep only tracts that can touch the zoning layer (fast bbox prefilter,
# then true intersects)
tracts <- tracts[lengths(st_intersects(
  tracts, st_as_sfc(st_bbox(zon)))) > 0, ]
tracts <- tracts[lengths(st_intersects(tracts, zon)) > 0, ]
message(nrow(tracts), " tracts intersect the zoning layer.")

## ---- 3. intersect ------------------------------------------------------------
message("Intersecting zoning x tracts (several minutes) ...")
pieces_sf <- st_intersection(
  tracts |> select(tract_id),
  zon |> select(zone_group, jurisd, adu_yn, gh_yn))
pieces_sf <- st_collection_extract(pieces_sf, "POLYGON")
pieces_sf$piece_km2 <- as.numeric(st_area(pieces_sf)) / 1e6
pieces_sf <- pieces_sf |> filter(piece_km2 > 0)

## ---- 4. aggregate (DISSOLVED areas -- overlap-proof) -------------------------
# (a) tract x group: union first, then measure (kills stacked duplicates)
by_group <- pieces_sf |>
  group_by(tract_id, zone_group) |>
  summarise(.groups = "drop") |>                # sf summarise = st_union
  mutate(km2 = as.numeric(st_area(geometry)) / 1e6) |>
  st_drop_geometry()

# (b) true zoned coverage per tract: union of ALL pieces
cover <- pieces_sf |>
  group_by(tract_id) |>
  summarise(.groups = "drop") |>
  mutate(union_km2 = as.numeric(st_area(geometry)) / 1e6) |>
  st_drop_geometry()

zoned_tot <- by_group |>
  group_by(tract_id) |>
  summarise(zoned_km2 = sum(km2), .groups = "drop") |>   # sum of group areas
  left_join(cover, by = "tract_id") |>
  # cross-group excess: >1 where different jurisdictions zone the SAME
  # ground into DIFFERENT classes (shares below split it proportionally)
  mutate(excess_crossgroup = zoned_km2 / union_km2)

shares <- by_group |>
  left_join(zoned_tot |> select(tract_id, zoned_km2), by = "tract_id") |>
  mutate(share = 100 * km2 / zoned_km2) |>
  select(tract_id, zone_group, share) |>
  pivot_wider(names_from = zone_group, values_from = share,
              values_fill = 0, names_prefix = "pct_")
# guarantee every group column exists even if absent in the data
for (g in names(CO_ZONE_GROUPS)) {
  cn <- paste0("pct_", g)
  if (!cn %in% names(shares)) shares[[cn]] <- 0
}

# residential-land measures: reslow share (group areas), ADU / group homes
# (dissolved per flag level; unknown rows dropped from num AND denom --
# blank = not coded, not "no")
res_flag_area <- function(flag_col) {
  pieces_sf |>
    filter(zone_group %in% CO_RES_GROUPS, !is.na(.data[[flag_col]])) |>
    select(tract_id, flag = all_of(flag_col)) |>
    group_by(tract_id, flag) |>
    summarise(.groups = "drop") |>              # dissolve: stacked polygons count once
    mutate(km2 = as.numeric(st_area(geometry)) / 1e6) |>
    st_drop_geometry() |>
    group_by(tract_id) |>
    summarise(pct_yes = pct_of(sum(km2[flag == "yes"]), sum(km2)),
              .groups = "drop")
}
res_measures <- by_group |>
  filter(zone_group %in% CO_RES_GROUPS) |>
  group_by(tract_id) |>
  summarise(res_zoned_km2     = sum(km2),
            pct_reslow_of_res = pct_of(sum(km2[zone_group == "res_low"]),
                                       sum(km2)),
            .groups = "drop") |>
  left_join(res_flag_area("adu_yn") |> rename(pct_adu_res = pct_yes),
            by = "tract_id") |>
  left_join(res_flag_area("gh_yn") |> rename(pct_grouphome_res = pct_yes),
            by = "tract_id")

# zoning-mix entropy over the non-uncertain groups (normalized 0..1)
entropy <- by_group |>
  filter(zone_group %in% CO_ENTROPY_GROUPS) |>
  group_by(tract_id) |>
  summarise(zoning_entropy = {
    p <- km2 / sum(km2); p <- p[p > 0]
    if (length(p) <= 1) 0 else -sum(p * log(p)) / log(length(CO_ENTROPY_GROUPS))
  }, .groups = "drop")

# jurisdictions: dissolved areas; count only those with >=5% of zoned area
jur_areas <- pieces_sf |>
  group_by(tract_id, jurisd) |>
  summarise(.groups = "drop") |>
  mutate(km2 = as.numeric(st_area(geometry)) / 1e6) |>
  st_drop_geometry()
jur <- jur_areas |>
  group_by(tract_id) |>
  summarise(n_jurisd = sum(km2 >= 0.05 * sum(km2)),
            jurisd_main = jurisd[which.max(km2)], .groups = "drop")

co_zoning <- tracts |>
  st_drop_geometry() |>
  inner_join(zoned_tot,   by = "tract_id") |>
  left_join(shares,       by = "tract_id") |>
  left_join(res_measures, by = "tract_id") |>
  left_join(entropy,      by = "tract_id") |>
  left_join(jur,          by = "tract_id") |>
  mutate(
    cover_zoned  = pmin(union_km2 / poly_km2, 1),
    pct_job_zone = pct_commercial + pct_industrial + pct_mixed_res)

saveRDS(co_zoning, out_file)
saveRDS(tracts, file.path(DIR_CO_CLEAN, "co_tract_geom.rds"))

# jurisdiction-level summary (for the scatter figure; dissolved by class)
jur_sum <- zon |>
  group_by(jurisd, zone_group) |>
  summarise(.groups = "drop") |>
  mutate(km2 = as.numeric(st_area(geometry)) / 1e6) |>
  st_drop_geometry() |>
  group_by(jurisd) |>
  summarise(zoned_km2 = sum(km2),
            pct_res_low = pct_of(sum(km2[zone_group == "res_low"]),
                                 sum(km2)), .groups = "drop") |>
  left_join(
    zon |>
      filter(zone_group %in% CO_RES_GROUPS, !is.na(adu_yn)) |>
      group_by(jurisd, adu_yn) |>
      summarise(.groups = "drop") |>            # dissolve within jurisdiction x flag
      mutate(km2 = as.numeric(st_area(geometry)) / 1e6) |>
      st_drop_geometry() |>
      group_by(jurisd) |>
      summarise(pct_adu_res = pct_of(sum(km2[adu_yn == "yes"]), sum(km2)),
                .groups = "drop"),
    by = "jurisd")
saveRDS(jur_sum, file.path(DIR_CO_CLEAN, "co_zoning_jurisd.rds"))

message("Saved ", nrow(co_zoning), " tracts with zoning measures.")

## ---- DIAGNOSTICS -------------------------------------------------------------
write_codiag(
  co_zoning |>
    mutate(county_fips = substr(tract_id, 1, 5)) |>
    group_by(county_fips) |>
    summarise(n_tracts = n(),
              n_covered80 = sum(cover_zoned >= CO_MIN_ZONED_COVER),
              mean_cover = round(mean(cover_zoned), 3),
              mean_pct_res_low = round(mean(pct_res_low, na.rm = TRUE), 1)),
  "61_coverage_by_county")
# Verified expectation: high coverage in 08001 Adams, 08005 Arapahoe, 08031
# Denver, 08035 Douglas, 08059 Jefferson, 08013 Boulder, 08014 Broomfield;
# 08123 Weld partial (SW towns only); fringe (08039 Elbert, 08093 Park,
# 08019 Clear Creek, 08047 Gilpin) ~zero -> dropped by the coverage screen.

xg <- co_zoning |> filter(excess_crossgroup > 1.10)
if (nrow(xg) > 0)
  write_codiag(xg |> select(tract_id, poly_km2, union_km2, zoned_km2,
                            excess_crossgroup) |>
                 arrange(desc(excess_crossgroup)), "61_crossgroup_conflicts")
# Cross-group excess = overlapping jurisdictions assigning DIFFERENT classes
# to the same ground (e.g. Foxfield vs Arapahoe County). Shares split these
# proportionally. Verified Aug 2026: p99 ~ 1.04, a handful of tracts >1.10.
if (any(co_zoning$excess_crossgroup > 1.5))
  warning("Some tracts have >50% cross-group overlap -- inspect ",
          "diagnostics/61_crossgroup_conflicts.csv before trusting shares.")

write_codiag(
  co_zoning |>
    filter(cover_zoned >= CO_MIN_ZONED_COVER) |>
    summarise(across(c(pct_res_low, pct_reslow_of_res, pct_adu_res,
                       pct_grouphome_res, zoning_entropy, pct_job_zone),
                     list(q25 = ~ quantile(.x, .25, na.rm = TRUE),
                          med = ~ quantile(.x, .50, na.rm = TRUE),
                          q75 = ~ quantile(.x, .75, na.rm = TRUE)))) |>
    pivot_longer(everything(), names_to = c("measure", "stat"),
                 names_pattern = "(.*)_(q25|med|q75)") |>
    pivot_wider(names_from = stat, values_from = value),
  "61_measure_distributions")
# Verified rough landmarks (covered tracts): pct_res_low p10/50/90 ~
# 0/57/92; pct_adu_res p10/50/90 ~ 0/89/100 (Denver's 2023 ADU expansion).

co_diag_map(readRDS(file.path(DIR_CO_CLEAN, "co_tract_geom.rds")),
            co_zoning |> select(tract_id, pct_res_low),
            "tract_id", "pct_res_low", "61_map_pct_res_low",
            "Share of zoned area: Residential_Low (exclusionary proxy)")
}
message("61 complete.")
