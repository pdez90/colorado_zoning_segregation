# ==============================================================================
# 32_segregation_panel.R
# THE keystone output: tract-year spatial dissimilarity, residential (RAC) and
# workplace (WAC) sides, for 3 measures x 4 decay settings, 2011-2023.
#
# Produces:
#   output/paper2/p2_tract_segregation_panel.rds
#     tract_id, year, d_{measure}_{rac|wac}_{quarter|half|one|aspatial}
#     (3 measures x 2 sides x 4 powers = 24 value columns)
#
# 60_co_setup.R asserts the conventions this file uses; 62 asserts every one
# of the 24 columns exists. Nothing downstream recomputes D except 72, which
# re-runs the same function at 5 km and 20 km cutoffs.
#
# SPEED NOTE. The honest implementation is compute_spatial_D() in 00, called
# once per tract x measure x power. That recomputes the same block distance
# matrix 12 times per tract. This script computes the distance matrix ONCE per
# tract and reuses it, then PROVES the shortcut agrees with compute_spatial_D
# to 1e-10 on a random sample of tracts before writing anything. If the check
# ever fails the script stops: 00's function is the definition, this is only
# an optimisation of it.
# ==============================================================================

source("30_p2_setup.R")

BLOCK_CENT <- file.path(DIR_CLEAN, "block_centroids_km.rds")
if (!file.exists(BLOCK_CENT))
  stop("Run 10_geography.R first -- ", BLOCK_CENT, " is missing.")
block_cent <- readRDS(BLOCK_CENT) |> select(block_id, X_km, Y_km)

# power label -> beta. NA = aspatial (each block is its own environment).
POWERS <- c(BETAS, list(aspatial = NA_real_))
MEAS   <- P2_MEASURES

## ---- one tract, all measures x all powers ------------------------------------
# Mirrors compute_spatial_D exactly; the only change is that the distance
# matrix and the weight matrices are built once and shared.
seg_one_tract <- function(xy, N_list, powers, maxdist_km = MAXDIST_KM) {
  n <- nrow(xy)

  out <- matrix(NA_real_, nrow = length(N_list), ncol = length(powers),
                dimnames = list(names(N_list), names(powers)))

  # single populated block -> no internal geography -> NA, matching
  # compute_spatial_D. This clause is the whole reason the fast path has to be
  # checked against the definition on DEGENERATE tracts and not just random
  # ones: single-block tracts are about 0.1% of the panel, so a 25-tract
  # random sample will essentially never contain one, and the two paths can
  # disagree on exactly the tracts that matter while the check reports 0.
  if (n < 2) return(out)

  spatial_needed <- any(!is.na(unlist(powers)))
  d <- if (spatial_needed && n > 1) as.matrix(stats::dist(as.matrix(xy))) else NULL

  W_list <- lapply(powers, function(b) {
    if (is.na(b) || is.null(d)) return(NULL)      # aspatial, or a lone block
    W <- exp(-b * d)
    W[d > maxdist_km] <- 0
    diag(W) <- 1                                   # d = 0 -> weight 1
    W
  })

  for (mi in seq_along(N_list)) {
    N <- N_list[[mi]]
    storage.mode(N) <- "double"
    N[!is.finite(N)] <- 0
    tau <- rowSums(N)
    T_j <- sum(tau)
    if (!is.finite(T_j) || T_j <= 0) next
    pi_m <- colSums(N) / T_j
    I_j  <- sum(pi_m * (1 - pi_m))
    if (!is.finite(I_j) || I_j <= 0) next          # one group only -> NA
    for (pj in seq_along(powers)) {
      W <- W_list[[pj]]
      local <- if (is.null(W)) N else W %*% N
      loc_tot <- rowSums(local)
      ok <- loc_tot > 0
      if (!any(ok)) next
      pit <- local[ok, , drop = FALSE] / loc_tot[ok]
      out[mi, pj] <- sum(tau[ok] * rowSums(abs(sweep(pit, 2, pi_m, "-")))) /
                     (2 * T_j * I_j)
    }
  }
  out
}

## ---- one side (rac/wac) of one year ------------------------------------------
seg_one_side <- function(blocks, side) {
  need <- unique(unlist(MEAS))
  miss <- setdiff(need, names(blocks))
  if (length(miss))
    stop("Block file is missing margin columns: ", paste(miss, collapse = ", "))

  pts <- blocks |>
    inner_join(block_cent, by = "block_id") |>
    filter(is.finite(X_km), is.finite(Y_km))

  idx <- split(seq_len(nrow(pts)), pts$tract_id)
  XY  <- as.matrix(pts[, c("X_km", "Y_km")])
  CNT <- lapply(MEAS, function(cols) as.matrix(pts[, cols]))

  vals <- vapply(idx, function(i) {
    as.vector(seg_one_tract(XY[i, , drop = FALSE],
                            lapply(CNT, function(M) M[i, , drop = FALSE]),
                            POWERS))
  }, numeric(length(MEAS) * length(POWERS)))

  nm <- as.vector(outer(names(MEAS), names(POWERS),
                        function(m, p) sprintf("d_%s_%s_%s", m, side, p)))
  out <- as_tibble(t(vals))
  names(out) <- nm
  out$tract_id <- names(idx)
  out |> relocate(tract_id)
}

## ---- the shortcut must agree with 00's definition ----------------------------
verify_against_definition <- function(blocks, n_tracts = 25, seed = 20260830) {
  set.seed(seed)
  pts <- blocks |> inner_join(block_cent, by = "block_id") |>
    filter(is.finite(X_km), is.finite(Y_km))

  # STRATIFY. A uniform random sample checks the easy middle of the
  # distribution and nothing else. The paths can only diverge on degenerate
  # tracts, so take every tract with fewer than three blocks, every tract
  # where one of the two groups is absent, and a random sample on top.
  prof <- pts |>
    group_by(tract_id) |>
    summarise(nb = n(),
              g1 = sum(.data[[MEAS$whiteblack[1]]], na.rm = TRUE),
              g2 = sum(.data[[MEAS$whiteblack[2]]], na.rm = TRUE),
              .groups = "drop")
  degenerate <- prof |> filter(nb < 3 | g1 == 0 | g2 == 0) |> pull(tract_id)
  ids <- unique(c(degenerate,
                  sample(prof$tract_id, min(n_tracts, nrow(prof)))))
  message("  verifying on ", length(ids), " tracts (",
          length(degenerate), " degenerate, ", n_tracts, " random)")
  worst <- 0
  for (tid in ids) {
    p <- pts[pts$tract_id == tid, ]
    fast <- seg_one_tract(as.matrix(p[, c("X_km", "Y_km")]),
                          lapply(MEAS, function(cols) as.matrix(p[, cols])),
                          POWERS)
    for (m in names(MEAS)) for (pw in names(POWERS)) {
      slow <- compute_spatial_D(p[, c("X_km", "Y_km")], p[, MEAS[[m]]],
                                POWERS[[pw]])
      a <- fast[m, pw]
      if (is.na(a) && is.na(slow)) next
      if (is.na(a) || is.na(slow))
        stop("Shortcut/definition disagree (NA) at tract ", tid, " ", m, " ", pw)
      worst <- max(worst, abs(a - slow))
    }
  }
  if (worst > 1e-10)
    stop("Shortcut differs from compute_spatial_D by ", worst, " -- STOP.")
  message("  verification vs compute_spatial_D: max |diff| = ",
          format(worst, scientific = TRUE))
  invisible(worst)
}

## ---- run ---------------------------------------------------------------------
dir.create(DIR_P2_OUT, showWarnings = FALSE, recursive = TRUE)
verified <- FALSE
panel <- list()

for (yr in P2_YEARS) {
  ck <- file.path(DIR_P2_OUT, sprintf("p2_seg_%d.rds", yr))
  if (file.exists(ck)) {
    message("  p2_seg_", yr, ".rds exists")
    panel[[as.character(yr)]] <- readRDS(ck)
    next
  }
  sides <- list()
  for (side in c("rac", "wac")) {
    f <- file.path(DIR_CLEAN, sprintf("blocks_%s_%s.rds", side, yr))
    if (!file.exists(f)) stop("Missing ", f, " -- run 20_lodes_blocks.R.")
    blocks <- readRDS(f)
    if (!verified) { verify_against_definition(blocks); verified <- TRUE }
    message("Segregation ", yr, " ", toupper(side), " (",
            format(nrow(blocks), big.mark = ","), " blocks) ...")
    sides[[side]] <- seg_one_side(blocks, side)
    rm(blocks); gc(verbose = FALSE)
  }
  yr_tbl <- full_join(sides$rac, sides$wac, by = "tract_id") |>
    mutate(year = as.integer(yr)) |>
    relocate(tract_id, year)
  saveRDS(yr_tbl, ck)
  message("  wrote p2_seg_", yr, ".rds (", nrow(yr_tbl), " tracts)")
  panel[[as.character(yr)]] <- yr_tbl
}

seg_panel <- bind_rows(panel) |> arrange(tract_id, year)

# the 24 columns 62 asserts
want <- with(expand.grid(m = names(MEAS), s = c("rac", "wac"),
                         p = names(POWERS), stringsAsFactors = FALSE),
             sprintf("d_%s_%s_%s", m, s, p))
stopifnot(all(want %in% names(seg_panel)))

saveRDS(seg_panel, P2_SEG_PANEL)
message("Wrote ", P2_SEG_PANEL, " (", nrow(seg_panel), " tract-years, ",
        n_distinct(seg_panel$tract_id), " tracts)")

write_diag(
  seg_panel |> group_by(year) |>
    summarise(n_tracts = n(),
              mean_d_wb_rac_half = mean(d_whiteblack_rac_half, na.rm = TRUE),
              mean_d_wb_wac_half = mean(d_whiteblack_wac_half, na.rm = TRUE),
              pctna_rac = round(100 * mean(is.na(d_whiteblack_rac_half)), 2),
              pctna_wac = round(100 * mean(is.na(d_whiteblack_wac_half)), 2)),
  "32_segregation_panel_by_year")

message("32_segregation_panel.R complete.")
