# ==============================================================================
# RUN_ALL.R
# Rebuilds every upstream input the Colorado zoning scripts expect.
#
#   cd ~/Downloads/LODES/paper_pipeline
#   Rscript RUN_ALL.R          # or source() it from an R session in that folder
#
# Rough timings on a laptop, everything cold:
#   10_geography     ~10 min   (TIGER blocks for Colorado is the slow download)
#   20_lodes_blocks  ~40 min   (13 years x RAC + WAC + OD main/aux)
#   32_segregation   ~25 min   (the index itself; the only CPU-bound step)
#   35_covariates    ~10 min   (ACS needs a Census API key)
#   53_sld           ~20 min   (~1 GB one-time download)
#
# Every step caches. Re-running skips whatever is already on disk, so an
# interrupted run picks up where it stopped. To force a step, delete its
# output and run again.
# ==============================================================================

setwd(path.expand("~/Downloads/LODES/paper_pipeline"))

STEPS <- c("10_geography.R",
           "20_lodes_blocks.R",
           "32_segregation_panel.R",
           "35_covariates.R",
           "53_sld.R")

for (s in STEPS) {
  message("\n", strrep("=", 78), "\n== ", s, "\n", strrep("=", 78))
  t0 <- Sys.time()
  ok <- tryCatch({ source(s, echo = FALSE); TRUE },
                 error = function(e) { message("!! ", s, " FAILED: ",
                                               conditionMessage(e)); FALSE })
  message("== ", s, if (ok) " ok " else " FAILED ",
          "(", round(difftime(Sys.time(), t0, units = "mins"), 1), " min)")
  if (!ok && s != "53_sld.R")
    stop("Stopping: ", s, " is required by everything after it.")
}

## ---- what should now exist ---------------------------------------------------
NEEDED <- c(
  file.path(DIR_CLEAN,  "tract_centroids_km.rds"),
  file.path(DIR_CLEAN,  "block_centroids_km.rds"),
  file.path(DIR_CLEAN,  "tract_aland_2020.rds"),
  file.path(DIR_CLEAN,  "p2_tract_lodes_cov_panel_v3.rds"),
  file.path(DIR_CLEAN,  "p2_tract_income_panel.rds"),
  file.path(DIR_P2_OUT, "p2_tract_segregation_panel.rds"))
OPTIONAL <- file.path(DIR_CLEAN, "p3_tract_sld.rds")

message("\n", strrep("=", 78))
for (f in NEEDED)
  message(if (file.exists(f)) "  ok      " else "  MISSING ", f)
message(if (file.exists(OPTIONAL)) "  ok      " else "  missing ", OPTIONAL,
        "  (network-accessibility results need this)")

if (all(file.exists(NEEDED)))
  message("\nUpstream rebuilt. Next, in ~/Downloads/LODES/Colorado:\n",
          "  for (s in sprintf('%d_*.R', 61:87)) ...  # or run 61,62,63,64 first\n",
          "then source('90_validate_rebuild.R') from paper_pipeline/ to check\n",
          "the rebuild against the committed results.")
