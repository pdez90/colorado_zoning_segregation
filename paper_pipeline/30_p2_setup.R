# ==============================================================================
# 30_p2_setup.R
# Paper-2 layer: which groups are measured, which years, and where the
# segregation panel lives. Sources 00 so a single source("50_p3_setup.R")
# from 60_co_setup.R pulls in the whole chain, as before.
# ==============================================================================

if (!exists("compute_spatial_D")) source("00_setup_and_functions.R")

## ---- group definitions -------------------------------------------------------
# ORDER MATTERS ONLY FOR READABILITY: the index is symmetric in the two groups.
# 60_co_setup.R asserts whiteblack == c("CR01","CR02") and
# hisp_nonhisp == c("CT02","CT01"); do not reorder those two.
P2_MEASURES <- list(
  whiteblack     = c("CR01", "CR02"),   # White alone; Black alone
  hisp_nonhisp   = c("CT02", "CT01"),   # Hispanic; non-Hispanic
  college_lesshs = c("CD04", "CD01")    # BA+; less than high school (age 30+)
)

# LODES margin columns we need from the block files.
# CR/CT/CD feed the segregation index; CE feeds the low-wage covariate;
# CNS05 (manufacturing) is WAC-only and feeds pct_manuf_wac in 64's BLOCK ZC.
# (The SI01-03 industry supergroups live in the OD files, not WAC; 66 and 71
# read those directly.)
P2_RAC_COLS <- c("C000", "CR01", "CR02", "CT01", "CT02",
                 "CD01", "CD04", "CE01", "CE02", "CE03")
P2_WAC_COLS <- c(P2_RAC_COLS, "CNS05")

P2_YEARS    <- 2011:2023
P2_PANDEMIC <- c(2020L, 2021L)          # flagged, never dropped
P2_OUT      <- DIR_P2_OUT

P2_SEG_PANEL <- file.path(DIR_P2_OUT, "p2_tract_segregation_panel.rds")

message("30_p2_setup.R loaded | panel ", min(P2_YEARS), "-", max(P2_YEARS),
        " | measures: ", paste(names(P2_MEASURES), collapse = ", "))
