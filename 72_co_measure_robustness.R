# ==============================================================================
# 72_co_measure_robustness.R      [PAPER 4, step 12]
# Sensitivity to the GENERATED SPATIAL MEASURES: is the headline
# interaction an artifact of the kernel/cutoff choices?
#
#  PART A  Decay-parameter sensitivity (CHEAP -- no recomputation).
#          2023 is an anchor year, so the Paper 2 panel already carries D at
#          beta = 0.25 ("quarter"), 0.5 ("half"), 1.0 ("one") and the
#          aspatial index. Rebuild wexp under each, refit the naive and
#          preferred (position-adjusted) specifications.
#          Expectation: the preferred estimate is similar across the three
#          spatial decays; only the aspatial index -- which discards the spatial structure the
#          construct is about -- is materially weaker.
#
#  PART B  Max-distance cutoff sensitivity (LOCAL RECOMPUTATION, ~minutes).
#          The panel is computed with a 10-km hard cutoff only, so 5-km and
#          20-km variants require re-running the block-level segregation for
#          the Denver-region tracts (2023, whiteblack, beta = 0.5, RAC+WAC).
#          Uses 31's block files + block centroids; ~1,400 CO tracts, so
#          this is minutes, not the overnight national job.
#
# Output: output/models/p4_measure_sensitivity.csv
#         clean/co_seg_maxdist_{5,20}_2023.rds
# ==============================================================================

source("60_co_setup.R")
library(fixest)

seg   <- readRDS(file.path(DIR_P2_OUT, "p2_tract_segregation_panel.rds"))
s23   <- seg |> filter(year == CO_ANCHOR_YEAR); rm(seg)
cent  <- readRDS(file.path(DIR_CLEAN, "tract_centroids_km.rds")) |>
  transmute(tract_id = as.character(GEOID), CBSA_Code)
acc   <- readRDS(file.path(DIR_CO_CLEAN, "co_accessibility_2023.rds"))
dat   <- readRDS(file.path(DIR_CO_OUT, "co_analysis_home_panel.rds"))

read_co_od23 <- function() {
  map(P3_OD_PARTS, function(part) {
    for (f in c(p3_od_cache(part, CO_ANCHOR_YEAR, "co"),
                co_od_cache(part, CO_ANCHOR_YEAR)))
      if (file.exists(f)) return(readRDS(f))
    NULL
  }) |> compact() |> bind_rows() |>
    group_by(w_tract, h_tract) |>
    summarise(S000 = sum(S000, na.rm = TRUE), .groups = "drop") |>
    inner_join(cent, by = c("h_tract" = "tract_id")) |>
    rename(cbsa_h = CBSA_Code) |>
    inner_join(cent, by = c("w_tract" = "tract_id")) |>
    rename(cbsa_w = CBSA_Code) |>
    filter(cbsa_h %in% CO_CBSA_KEEP, cbsa_h == cbsa_w)
}
od <- read_co_od23()

wexp_from <- function(dcol_by_tract) {          # named vector: tract -> D_wac
  od |>
    mutate(d = dcol_by_tract[w_tract]) |>
    filter(!is.na(d), S000 > 0) |>
    group_by(tract_id = h_tract) |>
    summarise(wexp = sum(d * S000) / sum(S000), .groups = "drop")
}

zscore <- function(x) { x[!is.finite(x)] <- NA_real_; as.numeric(scale(x)) }
base_frame <- dat |>
  filter(year == CO_ANCHOR_YEAR, in_scope, denver_msa,
         n_commuters >= P3_MIN_COMMUTERS) |>
  left_join(acc, by = "tract_id") |>
  mutate(across(c(pct_reslow_of_res, pct_black_rac, pct_lowincome_rac,
                  log_worker_density_rac, income_percapita_k,
                  income_percapita_k_sq, dist_cbd_km, dist_empctr_km),
                zscore, .names = "z_{.col}"))
COVS <- paste(c("z_pct_black_rac", "z_pct_lowincome_rac",
                "z_log_worker_density_rac", "z_income_percapita_k",
                "z_income_percapita_k_sq"), collapse = " + ")

fit_pair <- function(frame, label) {
  frame <- frame |> mutate(z_y = zscore(wexp), z_x = zscore(res))
  POS <- "z_dist_cbd_km + z_dist_empctr_km + z_x:z_dist_cbd_km + z_x:z_dist_empctr_km"
  map_df(c(naive = sprintf(
             "z_y ~ z_x * z_pct_reslow_of_res + %s | county_fips", COVS),
           preferred = sprintf(
             "z_y ~ z_x * z_pct_reslow_of_res + %s + %s | county_fips",
             COVS, POS)),
         .id = "spec", function(fml) {
    fit <- tryCatch(feols(as.formula(fml), data = frame,
                          cluster = ~jurisd_main), error = function(e) NULL)
    if (is.null(fit)) return(tibble())
    ct <- summary(fit)$coeftable
    i <- grep("^z_x:z_pct_reslow_of_res$", rownames(ct))
    tibble(variant = label, estimate = ct[i, 1], std.error = ct[i, 2],
           p.value = ct[i, 4], n_obs = fit$nobs)
  })
}

## ---- PART A: decay parameters ------------------------------------------------
res <- list()
for (pw in c("quarter", "half", "one", "aspatial")) {
  d_wac <- setNames(s23[[sprintf("d_whiteblack_wac_%s", pw)]], s23$tract_id)
  d_rac <- setNames(s23[[sprintf("d_whiteblack_rac_%s", pw)]], s23$tract_id)
  fr <- base_frame |>
    inner_join(wexp_from(d_wac), by = "tract_id") |>
    mutate(res = d_rac[tract_id])
  res[[length(res) + 1]] <- fit_pair(fr, sprintf("decay_%s", pw))
}

## ---- PART B: maxdist cutoffs (block-level recompute, CO only) ----------------
block_cent <- readRDS(file.path(DIR_CLEAN, "block_centroids_km.rds"))
run_seg_maxdist <- function(pts, cols, beta, md) {  # 00's runner + explicit md
  split_idx <- split(seq_len(nrow(pts)), pts$unit_id)
  v <- vapply(split_idx, function(idx)
    compute_spatial_D(as.matrix(pts[idx, c("X_km", "Y_km")]),
                      pts[idx, cols, drop = FALSE], beta,
                      maxdist_km = md), numeric(1))
  tibble(tract_id = names(v), value = unname(v))
}
region_tracts <- cent |> filter(CBSA_Code %in% CO_CBSA_KEEP) |> pull(tract_id)

for (md in c(5, 20)) {
  ck <- file.path(DIR_CO_CLEAN, sprintf("co_seg_maxdist_%d_2023.rds", md))
  if (file.exists(ck)) { message(basename(ck), " exists"); next }
  out <- list()
  for (ds in c("rac", "wac")) {
    blocks <- readRDS(file.path(DIR_CLEAN,
                                sprintf("blocks_%s_%s.rds", ds, CO_ANCHOR_YEAR)))
    pts <- blocks |>
      filter(tract_id %in% region_tracts) |>
      inner_join(block_cent |> select(block_id, X_km, Y_km), by = "block_id") |>
      rename(unit_id = tract_id)
    message(sprintf("maxdist %d km | %s | %s tracts", md, ds,
                    n_distinct(pts$unit_id)))
    out[[ds]] <- run_seg_maxdist(pts, P2_MEASURES$whiteblack, BETAS[["half"]],
                                 md) |>
      rename(!!paste0("d_", ds) := value)
    rm(blocks, pts); gc()
  }
  saveRDS(full_join(out$rac, out$wac, by = "tract_id"), ck)
}
for (md in c(5, 20)) {
  sm <- readRDS(file.path(DIR_CO_CLEAN,
                          sprintf("co_seg_maxdist_%d_2023.rds", md)))
  fr <- base_frame |>
    inner_join(wexp_from(setNames(sm$d_wac, sm$tract_id)), by = "tract_id") |>
    mutate(res = setNames(sm$d_rac, sm$tract_id)[tract_id])
  res[[length(res) + 1]] <- fit_pair(fr, sprintf("maxdist_%dkm", md))
}

all_res <- bind_rows(res)
write.csv(all_res, file.path(DIR_CO_MOD, "p4_measure_sensitivity.csv"),
          row.names = FALSE)
print(as.data.frame(all_res))
# Read: the PREFERRED estimates across decay_quarter/half/one and
# maxdist_5/10(=main)/20 should sit in a tight band; that stability -- not
# any single p-value -- is what answers the generated-measure concern.
message("72 complete.")
