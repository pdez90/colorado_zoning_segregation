# ==============================================================================
# 60_co_setup.R      [COLORADO ZONING case study, step 0]
# Shared configuration for the Denver-region zoning add-on to PAPER 3:
# does land-use zoning structure the tract-level link between residential
# segregation and the workplace-segregation exposure of a tract's residents?
# (See CO_ZONING_DESIGN.md in this folder.)
#
# This pipeline lives in ~/Downloads/LODES/Colorado and REUSES the Paper 2/3
# machinery + outputs (nothing national is recomputed):
#   * paper_pipeline/50_p3_setup.R (which sources 30 -> 00): config, seg fns
#   * output/paper2/p2_tract_segregation_panel.rds   (script 32)  <- REQUIRED
#   * clean/tract_centroids_km.rds                   (script 21/02)
#   * clean/p2_tract_lodes_cov_panel.rds, p2_tract_income_panel.rds (35)
#   * clean/tract_aland_2020.rds                     (36)
#   * raw/od_tract_{part}_{yr}_co.rds                (51 caches, reused if
#     present; otherwise 62 downloads Colorado itself -- ~1 min/year)
#
# NEW data this pipeline adds:
#   * ~/Wellbeing/Zoning/ALL_DenverMSA_10.3.23.shp -- harmonized zoning
#     districts for the Denver region, October 2023 snapshot.
#
# ZONING-DATA FACTS (verified against the DBF, Aug 2026):
#   * 100,837 polygons, 51 jurisdictions, EPSG:32613 (WGS84 / UTM 13N).
#   * GenZone2 classes (exact strings): Residential_Low / _Med / _MedHigh /
#     _High, MixedUseRes_Low / _Med / _MedHigh / _High, MixedUse_Conditional,
#     Commercial, Industrial, OpenSpace, Civic, MobileHome, Uncertain.
#   * ADU + GroupHome allowance flags are messy free text ("Yes", "YEs",
#     "Yes'", "NO", "N", "") -> normalized in 61; blanks stay NA.
#   * Coverage is the DENVER REGION, not Colorado: Denver-Aurora-Lakewood
#     (CBSA 19740) jurisdictions PLUS Boulder-county (14500) and southwest
#     Weld-county (24540: Frederick, Firestone, Dacono, Mead, ...)
#     jurisdictions. Rural fringe counties of 19740 (Elbert, Park, Clear
#     Creek, Gilpin) are NOT in the file -> tract coverage screen below.
#   * SNAPSHOT: one point in time (Oct 2023) -> cross-sectional design at
#     the 2023 anchor + "does 2023 zoning predict 2011->2023 trajectories"
#     (zoning changes slowly, but state clearly: no causal claims -- zoning
#     both shapes and follows who lives where).
# ==============================================================================

## ---- Source the Paper 3 pipeline (relative sources need its wd) --------------
PIPE_DIR <- path.expand("~/Downloads/LODES/paper_pipeline")
if (!dir.exists(PIPE_DIR))
  stop("paper_pipeline not found at ", PIPE_DIR, " -- edit PIPE_DIR in 60.")
owd <- setwd(PIPE_DIR)
source("50_p3_setup.R")     # -> 30_p2_setup.R -> 00_setup_and_functions.R
setwd(owd)

## ---- Guard: the 2018-revision corrections must be in force -------------------
# Two errors were found and fixed in the 2018 paper's revision (see the
# alignment audit): (1) the "White" group was silently ALL NON-BLACK workers
# because CR01 was never retained -- White must be CR01 (White alone) and
# Black CR02 (Black alone), never a complement; (2) distance decay was applied
# to unprojected DEGREE coordinates -- the corrected convention is EPSG:5070,
# exp(-beta * d_km) with beta = 0.5, 10-km cutoff. Everything downstream of
# this file inherits both corrections from the Paper 2 machinery; these
# assertions stop the pipeline if that machinery ever drifts.
stopifnot(
  identical(P2_MEASURES$whiteblack,   c("CR01", "CR02")),
  identical(P2_MEASURES$hisp_nonhisp, c("CT02", "CT01")),
  BETAS[["half"]] == 0.5,
  MAXDIST_KM      == 10,
  CRS_METERS      == 5070
)
# LODES covariates must come from the CORRECTED cache (v3: pct_white built
# from CR01 + raw group counts). The un-versioned file predates the fix.
P4_LODES_COV_FILE <- file.path(DIR_CLEAN, "p2_tract_lodes_cov_panel_v3.rds")
if (!file.exists(P4_LODES_COV_FILE))
  stop("Corrected covariate cache not found:\n  ", P4_LODES_COV_FILE,
       "\nRun the corrected paper_pipeline/35_p2_covariates.R first -- do ",
       "NOT fall back to p2_tract_lodes_cov_panel.rds (its 'White' group ",
       "is all non-Black workers).")

## ---- Colorado-case-study paths (everything new lands HERE) -------------------
DIR_CO       <- path.expand("~/Downloads/LODES/Colorado")
DIR_CO_CLEAN <- file.path(DIR_CO, "clean")
DIR_CO_OUT   <- file.path(DIR_CO, "output")
DIR_CO_MOD   <- file.path(DIR_CO_OUT, "models")
DIR_CO_FIG   <- file.path(DIR_CO_OUT, "figures")
DIR_CO_DIAG  <- file.path(DIR_CO, "diagnostics")
for (d in c(DIR_CO_CLEAN, DIR_CO_OUT, DIR_CO_MOD, DIR_CO_FIG, DIR_CO_DIAG))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

# QC tables/maps go to the CASE-STUDY diagnostics folder, not the national one
write_codiag <- function(df, name) {
  path <- file.path(DIR_CO_DIAG, paste0(name, ".csv"))
  write.csv(df, path, row.names = FALSE)
  message("  [diag] ", path)
  invisible(df)
}
co_diag_map <- function(geom_sf, values_df, id_col, value_col, name,
                        title = name) {
  df <- dplyr::left_join(geom_sf, values_df,
                         by = stats::setNames(id_col, names(geom_sf)[1]))
  p <- ggplot2::ggplot(df) +
    ggplot2::geom_sf(ggplot2::aes(fill = .data[[value_col]]), color = NA) +
    ggplot2::scale_fill_viridis_c(na.value = "red") +
    ggplot2::labs(title = title,
                  subtitle = "red = missing (investigate if unexpected)") +
    ggplot2::theme_void()
  path <- file.path(DIR_CO_DIAG, paste0(name, ".png"))
  ggplot2::ggsave(path, p, width = 10, height = 6, dpi = 200, bg = "white")
  message("  [diag] ", path)
  invisible(p)
}

# NA-safe percentage share (zero denominators -> NA, never Inf/NaN)
pct_of <- function(num, den) ifelse(den > 0, 100 * num / den, NA_real_)

## ---- Inputs ------------------------------------------------------------------
CO_ZONING_SHP <- path.expand("~/Wellbeing/data/Zoning/ALL_DenverMSA_10.3.23.shp")
CO_TIGER_DIR  <- path.expand("~/Downloads/LODES/TIGER2024_TRACT_UNZIPPED/tl_2024_08_tract")

## ---- Scope -------------------------------------------------------------------
CO_STATE_FIPS <- "08"
CO_CBSA_MAIN  <- "19740"                 # Denver-Aurora-Lakewood
CO_CBSA_EXTRA <- c("14500", "24540")     # Boulder; Greeley (SW Weld jurisd.)
CO_CBSA_KEEP  <- c(CO_CBSA_MAIN, CO_CBSA_EXTRA)

CO_ANCHOR_YEAR <- 2023L                  # zoning snapshot vintage
CO_YEARS       <- P3_YEARS               # 2011:2023 (trajectories + trends)

# A tract enters the MODELS only if zoning polygons cover >= this share of
# its land area (uncovered land = jurisdiction absent from the file, so its
# zoning is UNKNOWN, not zero). Sensitivity at 0.50 in 64.
CO_MIN_ZONED_COVER <- 0.80

## ---- GenZone2 -> analysis groups (exact strings; 61 stops on any new one) ----
CO_ZONE_GROUPS <- list(
  res_low     = "Residential_Low",
  res_midhigh = c("Residential_Med", "Residential_MedHigh",
                  "Residential_High"),
  mobile_home = "MobileHome",
  mixed_res   = c("MixedUseRes_Low", "MixedUseRes_Med",
                  "MixedUseRes_MedHigh", "MixedUseRes_High",
                  "MixedUse_Conditional"),
  commercial  = "Commercial",
  industrial  = "Industrial",
  open_civic  = c("OpenSpace", "Civic"),
  uncertain   = "Uncertain"
)
# groups counted as "residential land" (ADU / GroupHome denominators)
CO_RES_GROUPS <- c("res_low", "res_midhigh", "mobile_home", "mixed_res")
# groups entering the zoning-mix entropy (uncertain excluded)
CO_ENTROPY_GROUPS <- c("res_low", "res_midhigh", "mobile_home", "mixed_res",
                       "commercial", "industrial", "open_civic")

## ---- Checkpoint naming -------------------------------------------------------
co_od_cache <- function(part, yr)        # Colorado per-year-part raw OD cache
  file.path(DIR_CO_CLEAN, sprintf("co_od_tract_%s_%s.rds", part, yr))
co_wexp_year <- function(yr)             # per-year weighted exposures
  file.path(DIR_CO_CLEAN, sprintf("co_wexp_%s.rds", yr))

message(sprintf(
  "60_co_setup.R loaded | zoning snapshot %d | CBSAs %s | min cover %.0f%%",
  CO_ANCHOR_YEAR, paste(CO_CBSA_KEEP, collapse = "+"),
  100 * CO_MIN_ZONED_COVER))

## ---- major employment centers: ONE definition for every script ---------------
# Top 2% of tracts in the three-CBSA study region by CO_ANCHOR_YEAR workplace
# jobs (LODES WAC C000); tracts with no recorded jobs count as zero. This is the
# rule 68 uses to build dist_empctr_km, so every descriptive that mentions
# "major employment centers" (71, 78, 80) refers to the same set of tracts.
co_employment_centers <- function() {
  cent <- readRDS(file.path(DIR_CLEAN, "tract_centroids_km.rds")) |>
    dplyr::filter(CBSA_Code %in% CO_CBSA_KEEP) |>
    dplyr::transmute(tract_id = as.character(GEOID))
  wac <- readRDS(file.path(DIR_CO_CLEAN,
                           sprintf("co_wac_tract_%d.rds", CO_ANCHOR_YEAR))) |>
    dplyr::transmute(tract_id = as.character(tract_id), jobs = C000)
  d <- cent |> dplyr::left_join(wac, by = "tract_id") |>
    dplyr::mutate(jobs = ifelse(is.na(jobs), 0, jobs))
  d$tract_id[d$jobs >= stats::quantile(d$jobs, 0.98)]
}

## ---- manuscript figures: titles live in the captions ---------------------------
# The figures listed here are placed in the manuscript under a full caption, so
# their in-plot title, subtitle and caption are dropped at save time. Everything
# else about the figure (data, panels, legends, size) is untouched. Set
# CO_FIG_TITLES <- TRUE before sourcing a script to keep the titles (e.g. for
# slides). Panel titles inside multi-panel figures are never removed.
if (!exists("CO_FIG_TITLES")) CO_FIG_TITLES <- FALSE
CO_STRIP_TITLES <- c(
  "p4_figS2_adu.png", "p4_figS3_accessibility.png",
  "co_fig5_trends_by_tercile.png", "p4_figF2_exemplar_pair.png",
  "p4_fig_opportunity_sorting.png", "p4_fig_transit_marginal.png",
  "p4_fig_cervero_quadrants.png", "p4_fig_lowwage_dependence.png",
  "p4_fig_balance_matching.png", "p4_fig_decentralization.png",
  "p4_fig_excess_commuting.png", "p4_fig_workbased_menu.png")
ggsave <- function(filename, plot = ggplot2::last_plot(), ...) {
  if (!CO_FIG_TITLES && basename(filename) %in% CO_STRIP_TITLES) {
    plot <- if (inherits(plot, "patchwork"))
      plot + patchwork::plot_annotation(title = NULL, subtitle = NULL, caption = NULL)
    else
      plot + ggplot2::labs(title = NULL, subtitle = NULL, caption = NULL)
  }
  ggplot2::ggsave(filename, plot, ...)
}
