# ==============================================================================
# 50_p3_setup.R
# Paper-3 layer: commuting-flow conventions. This is the file 60_co_setup.R
# sources; it pulls 30 -> 00 behind it, so the Colorado scripts need no edits.
# ==============================================================================

if (!exists("P2_MEASURES")) source("30_p2_setup.R")

## ---- OD conventions ----------------------------------------------------------
P3_YEARS     <- P2_YEARS                 # 2011:2023
P3_OD_PARTS  <- c("main", "aux")         # aux = jobs held by out-of-state residents
# SA01-03 (age <=29 / 30-54 / 55+) are needed too: 66_co_group_flows.R splits
# flows by ten groups including age, and 67 differences the age groups.
# Leaving them out silently yielded 7 groups instead of 10 and made 67 fail on
# `mean_dist_km_age_u30`.
P3_OD_COLS   <- c("S000", "SA01", "SA02", "SA03",
                  "SE01", "SE02", "SE03",
                  "SI01", "SI02", "SI03")
P3_MIN_COMMUTERS <- 20L                  # tract enters models at >= 20 commuters

# Which segregation series the flow-weighting reads out of the panel.
P3_SEG_MEASURES <- names(P2_MEASURES)                    # 3 measures
P3_SEG_POWERS   <- c("quarter", "half", "one", "aspatial")

## ---- cache paths -------------------------------------------------------------
# The national pipeline wrote per-state OD caches here; the Colorado scripts
# look for these first and fall back to their own co_od_cache().
p3_od_cache <- function(part, yr, state = "co")
  file.path(DIR_RAW, sprintf("od_tract_%s_%s_%s.rds", part, yr, state))

P3_SLD_FILE <- file.path(DIR_CLEAN, "p3_tract_sld.rds")

message("50_p3_setup.R loaded | OD panel ", min(P3_YEARS), "-", max(P3_YEARS),
        " | parts: ", paste(P3_OD_PARTS, collapse = "+"))
