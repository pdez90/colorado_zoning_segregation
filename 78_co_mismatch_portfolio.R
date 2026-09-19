# ==============================================================================
# 78_co_mismatch_portfolio.R      [PAPER 4, step 18 -- SI extensions 1 + 3]
# Two SI analyses from the labor-market-geography review:
#
#  PART A (ext. 1) -- JOBS-WORKERS MISMATCH. Do neighborhoods contain jobs that
#    match the workers who live there? For each tract, compare the composition
#    of resident workers (RAC) with the composition of local jobs (WAC):
#        mismatch = 0.5 * sum_g | share_g^RAC - share_g^WAC |
#    (total-variation distance; 0 = identical composition, 1 = disjoint) over
#    (a) the three earnings bands CE01-03 and (b) the twenty NAICS sectors
#    CNS01-20. Defined only where both sides have mass: tracts with >= 100
#    resident workers AND >= 100 jobs (n reported). Then: does exclusionary
#    zoning predict larger mismatch? Naive and position-adjusted rungs, same
#    conventions as the main models (county FE, jurisdiction clusters).
#
#  PART B (ext. 3) -- WHAT IS DIFFERENT ABOUT THE DESTINATIONS? Exclusionary
#    tracts reach as MANY destinations as anyone (90.4/101.7/89.1); this
#    characterizes what those destinations ARE. Per origin tract, flow-weighted
#    destination attributes: share of flows to major employment centers (top
#    2% of tracts by jobs, 68's definition), mean destination distance from
#    the CBD, mean destination low-wage share (CE01/C000), mean destination
#    industry entropy (CNS01-20, normalized Shannon), mean destination scale
#    (log10 jobs), and concentration in the primary destination jurisdiction.
#    Summarized by exclusionary tercile, commuter-weighted.
#
# Descriptive + decomposition language only, as everywhere in this paper.
# Output: output/models/p4_mismatch_tercile.csv
#         output/models/p4_mismatch_models.csv
#         output/models/p4_dest_portfolio_tercile.csv
# ==============================================================================

source("60_co_setup.R")
library(fixest)

## ---- shared frames -----------------------------------------------------------
dat <- readRDS(file.path(DIR_CO_OUT, "co_analysis_home_panel.rds"))
acc <- readRDS(file.path(DIR_CO_CLEAN, "co_accessibility_2023.rds"))
zscore <- function(x) { x[!is.finite(x)] <- NA_real_; as.numeric(scale(x)) }
CE_COLS  <- sprintf("CE%02d", 1:3)
CNS_COLS <- sprintf("CNS%02d", 1:20)
norm_entropy <- function(m) {
  tot <- rowSums(m)
  p <- m / pmax(tot, .Machine$double.eps)
  h <- -rowSums(ifelse(p > 0, p * log(p), 0))
  ifelse(tot > 0, h / log(ncol(m)), NA_real_)
}

rac <- readRDS(file.path(DIR_CO_CLEAN,
                         sprintf("co_rac_tract_%d.rds", CO_ANCHOR_YEAR)))
wac <- readRDS(file.path(DIR_CO_CLEAN,
                         sprintf("co_wac_tract_%d.rds", CO_ANCHOR_YEAR)))
stopifnot(all(c(CE_COLS, CNS_COLS) %in% names(rac)),
          all(c(CE_COLS, CNS_COLS) %in% names(wac)))

xs <- dat |>
  filter(year == CO_ANCHOR_YEAR, in_scope, denver_msa,
         n_commuters >= P3_MIN_COMMUTERS) |>
  left_join(acc, by = "tract_id") |>
  mutate(tercile = ntile(pct_reslow_of_res, 3))

## ---- PART A: composition mismatch --------------------------------------------
tv_dist <- function(a, b) {           # rows = tracts; a, b share matrices
  0.5 * rowSums(abs(a - b))
}
shares <- function(df, cols) {
  m <- as.matrix(df[, cols]); m[is.na(m)] <- 0
  m / pmax(rowSums(m), 1)
}
mm <- xs |>
  select(tract_id) |>
  left_join(rac |> select(tract_id, workers = C000, all_of(c(CE_COLS, CNS_COLS))),
            by = "tract_id") |>
  rename_with(~ paste0("rac_", .x), all_of(c(CE_COLS, CNS_COLS))) |>
  left_join(wac |> select(tract_id, jobs = C000, all_of(c(CE_COLS, CNS_COLS))),
            by = "tract_id") |>
  rename_with(~ paste0("wac_", .x), all_of(c(CE_COLS, CNS_COLS))) |>
  filter(workers >= 100, jobs >= 100)
mm$mismatch_earn <- tv_dist(shares(mm, paste0("rac_", CE_COLS)),
                            shares(mm, paste0("wac_", CE_COLS)))
mm$mismatch_ind  <- tv_dist(shares(mm, paste0("rac_", CNS_COLS)),
                            shares(mm, paste0("wac_", CNS_COLS)))
message(sprintf("Mismatch defined for %d tracts (>=100 workers AND >=100 jobs)",
                nrow(mm)))

ma <- xs |>
  inner_join(mm |> select(tract_id, mismatch_earn, mismatch_ind),
             by = "tract_id") |>
  mutate(across(c(mismatch_earn, mismatch_ind, pct_reslow_of_res,
                  pct_black_rac, pct_lowincome_rac, log_worker_density_rac,
                  income_percapita_k, income_percapita_k_sq,
                  dist_cbd_km, dist_empctr_km),
                zscore, .names = "z_{.col}"))
terc_a <- ma |>
  group_by(tercile) |>
  summarise(n = n(),
            mismatch_earnings = mean(mismatch_earn),
            mismatch_industry = mean(mismatch_ind), .groups = "drop")
write.csv(terc_a, file.path(DIR_CO_MOD, "p4_mismatch_tercile.csv"),
          row.names = FALSE)
print(as.data.frame(terc_a), digits = 3)

COVS <- paste(c("z_pct_black_rac", "z_pct_lowincome_rac",
                "z_log_worker_density_rac", "z_income_percapita_k",
                "z_income_percapita_k_sq"), collapse = " + ")
POS  <- "z_dist_cbd_km + z_dist_empctr_km"
tidy1 <- function(fit, id, keep) {
  ct <- as.data.frame(summary(fit)$coeftable)
  ct$term <- rownames(ct); rownames(ct) <- NULL
  names(ct)[1:4] <- c("estimate", "std.error", "statistic", "p.value")
  ct |> filter(term %in% keep) |> mutate(model_id = id, n_obs = fit$nobs)
}
res <- list()
for (y in c("mismatch_earn", "mismatch_ind")) {
  res[[length(res) + 1]] <- tidy1(feols(as.formula(sprintf(
    "z_%s ~ z_pct_reslow_of_res + %s | county_fips", y, COVS)),
    data = ma, cluster = ~jurisd_main),
    sprintf("MM_%s_naive", y), "z_pct_reslow_of_res")
  res[[length(res) + 1]] <- tidy1(feols(as.formula(sprintf(
    "z_%s ~ z_pct_reslow_of_res + %s + %s | county_fips", y, COVS, POS)),
    data = ma, cluster = ~jurisd_main),
    sprintf("MM_%s_position", y), "z_pct_reslow_of_res")
}
mm_models <- bind_rows(res)
write.csv(mm_models, file.path(DIR_CO_MOD, "p4_mismatch_models.csv"),
          row.names = FALSE)
print(as.data.frame(mm_models), digits = 3)

## ---- PART B: destination-portfolio composition -------------------------------
cent <- readRDS(file.path(DIR_CLEAN, "tract_centroids_km.rds")) |>
  transmute(tract_id = as.character(GEOID), X_km, Y_km)
od <- map(P3_OD_PARTS, function(part) {
  for (f in c(p3_od_cache(part, CO_ANCHOR_YEAR, "co"),
              co_od_cache(part, CO_ANCHOR_YEAR)))
    if (file.exists(f)) return(readRDS(f))
  NULL
}) |> compact() |> bind_rows() |>
  group_by(w_tract, h_tract) |>
  summarise(S000 = sum(S000, na.rm = TRUE), .groups = "drop")

dest <- wac |>
  transmute(tract_id, jobs = C000,
            lowwage_share = ifelse(C000 > 0, CE01 / C000, NA_real_)) |>
  mutate(ind_entropy = norm_entropy({
    m <- as.matrix(wac[, CNS_COLS]); m[is.na(m)] <- 0; m })) |>
  left_join(acc |> select(tract_id, dist_cbd_km), by = "tract_id") |>
  mutate(empctr = tract_id %in% co_employment_centers())  # single definition, 60
message(sum(dest$empctr, na.rm = TRUE), " employment-center tracts (regional top 2%)")
jmap <- xs |> select(tract_id, jurisd_main)

# same frame as the workplace-exposure measure (Eq. 2): destinations inside the
# home tract's MSA, which for this sample is the Denver MSA
msa_tracts <- readRDS(file.path(DIR_CLEAN, "tract_centroids_km.rds")) |>
  filter(CBSA_Code == CO_CBSA_MAIN) |> pull(GEOID) |> as.character()
odd <- od |>
  filter(w_tract %in% msa_tracts) |>
  inner_join(xs |> select(tract_id, tercile), by = c(h_tract = "tract_id")) |>
  inner_join(dest, by = c(w_tract = "tract_id")) |>
  left_join(jmap |> rename(w_jur = jurisd_main), by = c(w_tract = "tract_id"))
port <- odd |>
  group_by(h_tract, tercile) |>
  summarise(
    flows = sum(S000),
    pct_to_empctr = 100 * sum(S000[empctr], na.rm = TRUE) / sum(S000),
    dest_cbd_km = weighted.mean(dist_cbd_km, S000, na.rm = TRUE),
    dest_lowwage_share = 100 * weighted.mean(lowwage_share, S000, na.rm = TRUE),
    dest_ind_entropy = weighted.mean(ind_entropy, S000, na.rm = TRUE),
    dest_log10_jobs = weighted.mean(log10(pmax(jobs, 1)), S000, na.rm = TRUE),
    pct_primary_jur = {
      ok <- !is.na(w_jur)
      if (any(ok)) 100 * max(tapply(S000[ok], w_jur[ok], sum)) / sum(S000)
      else NA_real_ }, .groups = "drop")
terc_b <- port |>
  group_by(tercile) |>
  summarise(n = n(),
            across(c(pct_to_empctr, dest_cbd_km, dest_lowwage_share,
                     dest_ind_entropy, dest_log10_jobs, pct_primary_jur),
                   ~ weighted.mean(.x, flows, na.rm = TRUE)), .groups = "drop")
write.csv(terc_b, file.path(DIR_CO_MOD, "p4_dest_portfolio_tercile.csv"),
          row.names = FALSE)
print(as.data.frame(terc_b), digits = 3)
message("78 complete.")
