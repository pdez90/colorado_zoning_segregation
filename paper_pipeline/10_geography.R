# ==============================================================================
# 10_geography.R
# Block and tract geography for Colorado, projected to EPSG:5070, in KM.
#
# Produces:
#   clean/block_centroids_km.rds   block_id, tract_id, X_km, Y_km
#   clean/tract_centroids_km.rds   GEOID, tract_id, CBSA_Code, X_km, Y_km
#   clean/tract_aland_2020.rds     tract_id, aland_km2
#
# CBSA_Code comes from the Census delineation crosswalk so the Colorado scripts
# can filter to 19740 / 14500 / 24540 exactly as before.
# ==============================================================================

source("50_p3_setup.R")
suppressPackageStartupMessages(library(tigris))
options(tigris_use_cache = TRUE, tigris_class = "sf")

STATE <- "08"

## ---- tract geography ---------------------------------------------------------
tr_file <- file.path(DIR_CLEAN, "tract_centroids_km.rds")
al_file <- file.path(DIR_CLEAN, "tract_aland_2020.rds")

if (!file.exists(tr_file) || !file.exists(al_file)) {
  message("Downloading TIGER tracts (CO, 2024)...")
  tr <- tigris::tracts(state = STATE, year = 2024, progress_bar = FALSE)
  tr <- sf::st_transform(tr, CRS_METERS)
  xy <- suppressWarnings(sf::st_coordinates(sf::st_centroid(tr)))

  # county -> CBSA crosswalk from the Census delineation file
  cbsa <- tryCatch({
    d <- readr::read_csv(
      "https://www2.census.gov/programs-surveys/metro-micro/geographies/reference-files/2023/delineation-files/list1_2023.csv",
      skip = 2, show_col_types = FALSE)
    d |>
      transmute(county_fips = paste0(sprintf("%02d", as.integer(`FIPS State Code`)),
                                     sprintf("%03d", as.integer(`FIPS County Code`))),
                CBSA_Code = as.character(`CBSA Code`)) |>
      filter(!is.na(county_fips))
  }, error = function(e) {
    message("  delineation download failed (", conditionMessage(e),
            ") -- falling back to the hard-coded Denver-region crosswalk")
    tibble(county_fips = c("08001","08005","08014","08019","08031","08035",
                           "08039","08047","08059","08093",          # 19740
                           "08013",                                   # 14500
                           "08123"),                                  # 24540
           CBSA_Code   = c(rep("19740", 10), "14500", "24540"))
  })

  tract_centroids <- tibble(
      GEOID    = as.character(tr$GEOID),
      tract_id = as.character(tr$GEOID),
      county_fips = substr(as.character(tr$GEOID), 1, 5),
      X_km = xy[, 1] / 1000, Y_km = xy[, 2] / 1000,
      aland_km2 = as.numeric(tr$ALAND) / 1e6) |>
    left_join(cbsa, by = "county_fips")

  saveRDS(tract_centroids |> select(GEOID, tract_id, county_fips,
                                    CBSA_Code, X_km, Y_km), tr_file)
  saveRDS(tract_centroids |> select(tract_id, aland_km2), al_file)
  message("  wrote tract_centroids_km.rds (", nrow(tract_centroids), " tracts)")
  write_diag(tract_centroids |> count(CBSA_Code), "10_tracts_by_cbsa")
} else message("  tract geography exists")

## ---- the TIGER folder 61_co_zoning_tract.R reads ------------------------------
# 61 prefers a TIGER 2024 tract shapefile on disk and falls back to tigris'
# 2023 vintage. The zoning measures are built on the 2024 vintage, so write the
# 2024 layer where 60_co_setup.R's CO_TIGER_DIR points.
tig_dir <- path.expand("~/Downloads/LODES/TIGER2024_TRACT_UNZIPPED/tl_2024_08_tract")
tig_shp <- file.path(tig_dir, "tl_2024_08_tract.shp")
if (!file.exists(tig_shp)) {
  dir.create(tig_dir, showWarnings = FALSE, recursive = TRUE)
  tr24 <- tigris::tracts(state = STATE, year = 2024, progress_bar = FALSE)
  sf::st_write(tr24, tig_shp, quiet = TRUE, delete_dsn = TRUE)
  message("  wrote ", tig_shp)
} else message("  TIGER 2024 tract shapefile exists")

## ---- block geography ---------------------------------------------------------
bl_file <- file.path(DIR_CLEAN, "block_centroids_km.rds")
if (!file.exists(bl_file)) {
  message("Downloading TIGER blocks (CO, 2020) -- this is the slow one...")
  bl <- tigris::blocks(state = STATE, year = 2020, progress_bar = FALSE)
  bl <- sf::st_transform(bl, CRS_METERS)
  bxy <- suppressWarnings(sf::st_coordinates(sf::st_centroid(bl)))
  gid <- as.character(bl$GEOID20)
  block_centroids <- tibble(block_id = gid,
                            tract_id = substr(gid, 1, 11),
                            X_km = bxy[, 1] / 1000, Y_km = bxy[, 2] / 1000)
  saveRDS(block_centroids, bl_file)
  message("  wrote block_centroids_km.rds (", nrow(block_centroids), " blocks)")
} else message("  block geography exists")

message("10_geography.R complete.")
