# ==============================================================================
# 00_setup_and_functions.R
#
# Core conventions and the spatial segregation index. Everything downstream of
# 60_co_setup.R inherits the two corrections that file asserts:
#   (1) White = CR01 (White alone), Black = CR02 — never a non-Black complement.
#   (2) Distance decay on PROJECTED coordinates (EPSG:5070), exp(-beta * d_km),
#       beta = 0.5, 10-km truncation.
# The guards in 60 check exactly these constants, so do not edit them casually.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(tibble)
  library(sf); library(data.table)
  # The figure scripts (65, 70, 74, 75, 76) call ggplot(), theme_void() and
  # friends unqualified, so ggplot2 is attached here.
  # Leaving it out made all five fail with "could not find function ggplot".
  library(ggplot2); library(scales)
  # optional; attach only if present so a missing one is not fatal
  for (p in c("ggrepel", "patchwork"))
    if (requireNamespace(p, quietly = TRUE))
      library(p, character.only = TRUE)
})

## ---- corrected conventions (asserted by 60_co_setup.R) -----------------------
CRS_METERS  <- 5070L          # NAD83 / Conus Albers — equal area, metres
MAXDIST_KM  <- 10             # weights below ~1% of their value at zero
BETAS       <- list(quarter = 0.25, half = 0.5, one = 1.0)

## ---- Reardon & O'Sullivan (2004) spatial dissimilarity -----------------------
# D_j = sum_b tau_b * sum_m |pi~_bm - pi_m| / (2 * T_j * I_j)
#
#   b        blocks of tract j (the local environment is WITHIN-tract: this is
#            what makes the index a measure of within-tract compositional
#            sorting, and what keeps residential and workplace measures from
#            overlapping mechanically)
#   tau_b    total of the two groups at block b;  T_j = sum_b tau_b
#   pi_m     tract-wide proportion of group m
#   pi~_bm   distance-weighted local proportion around b, weights
#            w_bk = exp(-beta * d_bk) for d_bk <= maxdist_km, else 0 (w_bb = 1)
#   I_j      interaction index, sum_m pi_m (1 - pi_m) = 2*pi_1*pi_2 for 2 groups
#
# beta = NA (or maxdist_km = 0) gives the ASPATIAL index: each block is its own
# local environment, which reduces to the classical index of dissimilarity.
#
# Returns NA when the tract has no workers in the two groups, only one group
# present (I_j = 0), or fewer than two blocks containing either group.
compute_spatial_D <- function(xy, counts, beta, maxdist_km = MAXDIST_KM) {
  N <- as.matrix(counts)
  storage.mode(N) <- "double"
  N[!is.finite(N)] <- 0

  tau <- rowSums(N)
  T_j <- sum(tau)
  if (!is.finite(T_j) || T_j <= 0) return(NA_real_)

  # A tract with fewer than two blocks holding members of these two groups has
  # no internal geography for this measure: the local environment IS the tract
  # and the formula collapses to 0, which would read as "perfectly integrated"
  # when it means "not measurable here". Return NA instead.
  if (sum(tau > 0) < 2) return(NA_real_)

  pi_m <- colSums(N) / T_j
  I_j  <- sum(pi_m * (1 - pi_m))
  if (!is.finite(I_j) || I_j <= 0) return(NA_real_)   # one group only

  aspatial <- is.na(beta) || is.null(beta) || maxdist_km <= 0
  if (aspatial) {
    local <- N
  } else {
    xy <- as.matrix(xy)
    storage.mode(xy) <- "double"
    d <- as.matrix(stats::dist(xy))                   # km, projected
    W <- exp(-beta * d)
    W[d > maxdist_km] <- 0
    diag(W) <- 1                                      # d = 0 -> weight 1
    local <- W %*% N
  }
  loc_tot <- rowSums(local)
  ok <- loc_tot > 0
  if (!any(ok)) return(NA_real_)
  pit <- local[ok, , drop = FALSE] / loc_tot[ok]

  num <- sum(tau[ok] * rowSums(abs(sweep(pit, 2, pi_m, "-"))))
  num / (2 * T_j * I_j)
}

# Run compute_spatial_D for every unit in a long block table.
# pts: block rows with unit_id, X_km, Y_km and the two group-count columns.
run_spatial_D <- function(pts, cols, beta, maxdist_km = MAXDIST_KM) {
  pts <- pts[stats::complete.cases(pts[, c("X_km", "Y_km")]), , drop = FALSE]
  idx <- split(seq_len(nrow(pts)), pts$unit_id)
  v <- vapply(idx, function(i)
    compute_spatial_D(pts[i, c("X_km", "Y_km"), drop = FALSE],
                      pts[i, cols, drop = FALSE], beta, maxdist_km),
    numeric(1))
  tibble(unit_id = names(v), value = unname(v))
}

## ---- small shared helpers ----------------------------------------------------
# Data root. Override with the LODES_ROOT environment variable.
DIR_ROOT   <- path.expand(Sys.getenv("LODES_ROOT", "~/Downloads/LODES"))
DIR_RAW    <- file.path(DIR_ROOT, "raw")
DIR_CLEAN  <- file.path(DIR_ROOT, "clean")
DIR_OUT    <- file.path(DIR_ROOT, "output")
DIR_P2_OUT <- file.path(DIR_OUT, "paper2")
DIR_DIAG   <- file.path(DIR_ROOT, "diagnostics")
for (d in c(DIR_RAW, DIR_CLEAN, DIR_OUT, DIR_P2_OUT, DIR_DIAG))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

write_diag <- function(df, name) {
  p <- file.path(DIR_DIAG, paste0(name, ".csv"))
  utils::write.csv(df, p, row.names = FALSE)
  message("  [diag] ", p); invisible(df)
}

# checkpoint helper: run expr only if the cache is absent
cached <- function(path, expr, label = basename(path)) {
  if (file.exists(path)) { message("  ", label, " exists"); return(readRDS(path)) }
  v <- force(expr)
  saveRDS(v, path)
  message("  wrote ", label)
  v
}

message("00_setup_and_functions.R loaded | CRS ", CRS_METERS,
        " | beta half = ", BETAS[["half"]], " | maxdist ", MAXDIST_KM, " km")
