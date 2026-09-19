# ==============================================================================
# 53_sld.R
# EPA Smart Location Database v3 (2010 block groups) -> 2020 census tracts.
#
# Produces
#   clean/p3_tract_sld.rds            the PRIMARY construction, read by 63/64/75
#     tract_id, sld_D5AR, sld_D5BR, sld_NatWalkInd, sld_pop
#   clean/p3_tract_sld_variants.rds   every construction, long, read by 75b
#     construction, bridge, nodata, tract_id, sld_D5AR, sld_D5BR
#   diagnostics/53_sld_summary.csv    one row per construction
#
#   D5AR  jobs within a 45-minute AUTO commute, time-decay weighted
#   D5BR  jobs within a 45-minute TRANSIT commute, time-decay weighted
#   NatWalkInd  EPA National Walkability Index
#
# ------------------------------------------------------------------------------
# METHOD (primary construction)
# ------------------------------------------------------------------------------
# The SLD is published on 2010 block groups; LODES 8 and everything downstream
# is on 2020 tracts. The bridge is population-weighted areal interpolation at
# the 2020 census block, the finest unit available:
#
#   1. every 2020 block is assigned to the 2010 block group that contains its
#      internal point (INTPTLAT20 / INTPTLON20, which TIGER guarantees to lie
#      inside the block);
#   2. a 2020 tract's value is the mean of its blocks' block-group values,
#      weighted by 2020 block population (POP20).
#
# This is the standard block-based crosswalk design (as in the NHGIS
# crosswalks). It is a single stage, it matches the SLD's own resident-facing,
# population-weighted-centroid construction, and boundary slivers carry no
# weight because they contain no block internal points. A tract with no
# residents falls back to block land area (ALAND20) so that every tract
# receives a value.
#
# -99999. EPA codes transit fields -99999 where a block group has no transit
# service within range. Within the Denver region, which has full GTFS coverage
# (RTD), that is a real zero: no jobs are reachable by transit. The primary
# construction therefore reads it as 0, so a tract's value is the access of its
# average resident, including residents without service. The alternative
# reading (missing: average only the served block groups) is carried as a
# sensitivity. For non-transit fields -99999 is always missing.
#
# SENSITIVITY CONSTRUCTIONS (all written to the variants file)
#   bridge   block_pop   (primary)  2020 blocks, POP20 weights
#            block_area             2020 blocks, ALAND20 weights
#            rel_area               Census 2010->2020 tract relationship file,
#                                   AREALAND_PART weights, all overlapping
#                                   parents (SLD first averaged to 2010 tracts
#                                   with TotPop weights)
#            rel_dominant           as rel_area, largest-overlap parent only
#   nodata   zero (primary) | missing
#
# Force a rebuild by deleting clean/p3_tract_sld.rds.
# SLD source: Sys.getenv("SLD_PATH"), else raw/, else download from EPA.
# ==============================================================================

source("50_p3_setup.R")
suppressPackageStartupMessages({ library(sf); library(tigris) })
options(tigris_use_cache = TRUE, tigris_class = "sf")

STATE        <- "08"
out_file     <- P3_SLD_FILE
variant_file <- file.path(DIR_CLEAN, "p3_tract_sld_variants.rds")
PRIMARY      <- "block_pop|zero"
NODATA       <- -9999          # anything at or below this is EPA's no-data code

if (file.exists(out_file) && file.exists(variant_file)) {
  message("  p3_tract_sld.rds and its variants exist -- delete to rebuild.")
} else {

## ---- 1. locate and read the SLD (Colorado rows only) -------------------------
CSV_URL <- paste0("https://edg.epa.gov/EPADataCommons/public/OA/",
                  "EPA_SmartLocationDatabase_V3_Jan_2021_Final.csv")
csv_dest <- file.path(DIR_RAW, "EPA_SmartLocationDatabase_V3.csv")
src <- path.expand(Sys.getenv("SLD_PATH", ""))
if (src == "") {
  if (!(file.exists(csv_dest) && file.size(csv_dest) > 1e6)) {
    message("Downloading the SLD (~200 MB) ...")
    utils::download.file(CSV_URL, csv_dest, mode = "wb", method = "libcurl")
  }
  src <- csv_dest
}
if (!file.exists(src) && !dir.exists(src)) stop("SLD not found at: ", src)

message("Reading ", basename(src), " ...")
sld <- if (grepl("[.]csv$", src, ignore.case = TRUE)) {
  as.data.frame(data.table::fread(src, showProgress = FALSE))
} else {
  sf::st_drop_geometry(sf::st_read(src, layer = sf::st_layers(src)$name[1],
                                   quiet = TRUE))
}

# column names in this file are inconsistently cased (D5AR but TotPop)
pick <- function(want) {
  hit <- names(sld)[tolower(names(sld)) == tolower(want)]
  if (!length(hit)) stop("SLD column not found: ", want)
  hit[1]
}
# GEOID10 is the 2010 TIGER block-group code (SLD data dictionary). A 12-digit
# id read as a number loses Colorado's leading zero; restore it via sprintf on
# a double, which is exact at this magnitude.
pad12 <- function(x) sprintf("%012.0f", suppressWarnings(as.numeric(as.character(x))))

sld <- tibble(
    bg10   = pad12(sld[[pick("GEOID10")]]),
    totpop = as.numeric(sld[[pick("TotPop")]]),
    d5ar   = as.numeric(sld[[pick("D5AR")]]),
    d5br   = as.numeric(sld[[pick("D5BR")]]),
    walk   = as.numeric(sld[[pick("NatWalkInd")]])) |>
  filter(substr(bg10, 1, 2) == STATE) |>
  mutate(totpop = ifelse(is.finite(totpop) & totpop > 0, totpop, 0),
         d5ar   = ifelse(d5ar <= NODATA, NA_real_, d5ar),
         walk   = ifelse(walk <= NODATA, NA_real_, walk),
         d5br_nodata = !is.na(d5br) & d5br <= NODATA)
stopifnot(nrow(sld) > 3000, !anyDuplicated(sld$bg10))
message("  ", nrow(sld), " Colorado block groups; ", sum(sld$d5br_nodata),
        " carry the transit no-data code")

## ---- 2. 2020 blocks -> 2010 block groups -------------------------------------
message("Loading 2020 blocks and 2010 block groups ...")
bl <- tigris::blocks(state = STATE, year = 2020, progress_bar = FALSE)
need <- c("GEOID20", "POP20", "ALAND20", "INTPTLAT20", "INTPTLON20")
if (!all(need %in% names(bl)))
  stop("TIGER 2020 blocks are missing: ",
       paste(setdiff(need, names(bl)), collapse = ", "))
blk <- tibble(block_id = as.character(bl$GEOID20),
              tract_id = substr(as.character(bl$GEOID20), 1, 11),
              pop      = as.numeric(bl$POP20),
              aland    = as.numeric(bl$ALAND20),
              lat      = as.numeric(bl$INTPTLAT20),
              lon      = as.numeric(bl$INTPTLON20))
pts <- sf::st_as_sf(blk, coords = c("lon", "lat"), crs = sf::st_crs(bl)) |>
  sf::st_transform(CRS_METERS)
rm(bl)

bg_shp <- file.path(DIR_RAW, "tl_2010_08_bg10.shp")
bg <- if (file.exists(bg_shp)) sf::st_read(bg_shp, quiet = TRUE) else
  tigris::block_groups(state = STATE, year = 2010, progress_bar = FALSE)
bg <- bg |> transmute(bg10 = as.character(GEOID10)) |> sf::st_transform(CRS_METERS)

hit <- sf::st_within(pts, bg)
n_hit <- lengths(hit)
if (any(n_hit > 1)) stop("A block internal point fell in more than one block group.")
blk$bg10 <- NA_character_
blk$bg10[n_hit == 1] <- bg$bg10[unlist(hit[n_hit == 1])]
# the few points exactly on the state outline take the nearest block group
if (any(n_hit == 0)) {
  miss <- which(n_hit == 0)
  blk$bg10[miss] <- bg$bg10[sf::st_nearest_feature(pts[miss, ], bg)]
  message("  ", length(miss), " blocks assigned by nearest block group")
}
stopifnot(!anyNA(blk$bg10))
message("  ", nrow(blk), " blocks assigned; total POP20 = ",
        format(sum(blk$pop), big.mark = ","))

## ---- 3. the constructions ----------------------------------------------------
wmean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) NA_real_ else sum(x[ok] * w[ok]) / sum(w[ok])
}
apply_nodata <- function(d, nodata)
  d |> mutate(d5br = ifelse(d5br_nodata,
                            if (nodata == "zero") 0 else NA_real_, d5br))

block_bridge <- function(nodata, weight) {
  b <- blk |> inner_join(apply_nodata(sld, nodata), by = "bg10")
  b$w <- if (weight == "pop") b$pop else b$aland
  out <- b |> group_by(tract_id) |>
    summarise(sld_D5AR = wmean(d5ar, w), sld_D5BR = wmean(d5br, w),
              sld_NatWalkInd = wmean(walk, w), sld_pop = sum(pop),
              # unpopulated tract: fall back to land area
              .a = wmean(d5ar, aland), .b = wmean(d5br, aland),
              .k = wmean(walk, aland), .groups = "drop")
  out |> mutate(sld_D5AR = coalesce(sld_D5AR, .a),
                sld_D5BR = coalesce(sld_D5BR, .b),
                sld_NatWalkInd = coalesce(sld_NatWalkInd, .k)) |>
    select(-.a, -.b, -.k)
}

rel_file <- file.path(DIR_RAW, "tract_rel_2010_2020_co.rds")
rel <- if (file.exists(rel_file)) readRDS(rel_file) else {
  message("Downloading the 2010->2020 tract relationship file ...")
  d <- readr::read_delim(paste0(
    "https://www2.census.gov/geo/docs/maps-data/data/rel2020/tract/",
    "tab20_tract20_tract10_st08.txt"), delim = "|", show_col_types = FALSE) |>
    transmute(tract20 = as.character(GEOID_TRACT_20),
              tract10 = as.character(GEOID_TRACT_10),
              w = as.numeric(AREALAND_PART)) |>
    filter(!is.na(tract20), !is.na(tract10), w > 0)
  saveRDS(d, rel_file); d
}
rel_bridge <- function(nodata, dominant) {
  t10 <- apply_nodata(sld, nodata) |>
    mutate(tract10 = substr(bg10, 1, 11),
           tw = ifelse(totpop > 0, totpop, 0)) |>
    group_by(tract10) |>
    summarise(a = coalesce(wmean(d5ar, tw), wmean(d5ar, rep(1, n()))),
              b = coalesce(wmean(d5br, tw), wmean(d5br, rep(1, n()))),
              k = coalesce(wmean(walk, tw), wmean(walk, rep(1, n()))),
              p = sum(totpop), .groups = "drop")
  r <- rel |> inner_join(t10, by = "tract10")
  if (dominant) r <- r |> group_by(tract20) |> slice_max(w, n = 1, with_ties = FALSE) |> ungroup()
  r |> group_by(tract_id = tract20) |>
    summarise(sld_D5AR = wmean(a, w), sld_D5BR = wmean(b, w),
              sld_NatWalkInd = wmean(k, w), sld_pop = sum(p), .groups = "drop")
}

grid <- tidyr::expand_grid(bridge = c("block_pop", "block_area",
                                      "rel_area", "rel_dominant"),
                           nodata = c("zero", "missing"))
variants <- purrr::pmap_dfr(grid, function(bridge, nodata) {
  t <- switch(bridge,
    block_pop    = block_bridge(nodata, "pop"),
    block_area   = block_bridge(nodata, "area"),
    rel_area     = rel_bridge(nodata, dominant = FALSE),
    rel_dominant = rel_bridge(nodata, dominant = TRUE))
  t |> mutate(construction = paste(bridge, nodata, sep = "|"),
              bridge = bridge, nodata = nodata, .before = 1)
})

## ---- 4. checks, write --------------------------------------------------------
primary <- variants |> filter(construction == PRIMARY) |>
  select(tract_id, sld_D5AR, sld_D5BR, sld_NatWalkInd, sld_pop)
cent <- readRDS(file.path(DIR_CLEAN, "tract_centroids_km.rds"))
if (any(variants$sld_D5AR < 0, na.rm = TRUE) || any(variants$sld_D5BR < 0, na.rm = TRUE))
  stop("Negative accessibility: the -99999 code leaked into an average.")
if (!all(cent$tract_id %in% primary$tract_id))
  stop(sum(!cent$tract_id %in% primary$tract_id),
       " 2020 tracts received no SLD value under the primary construction.")
if (anyNA(primary$sld_D5AR)) stop("Primary construction has NA auto access.")

saveRDS(primary,  out_file)
saveRDS(variants |> select(construction, bridge, nodata, tract_id,
                           sld_D5AR, sld_D5BR), variant_file)
message("  wrote p3_tract_sld.rds (", nrow(primary), " tracts, primary = ",
        PRIMARY, ") and ", dplyr::n_distinct(variants$construction),
        " constructions to p3_tract_sld_variants.rds")

write_diag(
  variants |> group_by(construction) |>
    summarise(n_tracts = n(),
              n_auto_missing = sum(is.na(sld_D5AR)),
              median_D5AR = round(median(sld_D5AR, na.rm = TRUE)),
              n_transit_positive = sum(sld_D5BR > 0, na.rm = TRUE),
              median_D5BR_positive = round(median(sld_D5BR[sld_D5BR > 0], na.rm = TRUE)),
              min_D5BR_positive = signif(min(sld_D5BR[sld_D5BR > 0], na.rm = TRUE), 3),
              .groups = "drop") |>
    mutate(primary = construction == PRIMARY),
  "53_sld_summary")
}

message("53_sld.R complete.")
