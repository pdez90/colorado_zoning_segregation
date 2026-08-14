# ==============================================================================
# 85_co_jhbalance_literature.R   [PAPER 4, step 25 -- SI: speaking to the
#                                 jobs-housing balance literature]
#
#  PART A -- JOBS-HOUSING BALANCE ON ITS OWN TERMS (Cervero 1989; Peng 1997;
#    Zheng et al. 2021). Two methodological claims from that literature are
#    tested here against OUR outcome (workplace-location segregation exposure):
#      (i)  SCALE. Peng argues the jobs-housing ratio must be measured in a
#           floating catchment rather than within fixed boundaries, and Zheng
#           et al. restate this as a modifiable-areal-unit problem. We compute
#           the ratio four ways -- own tract, 5 km catchment, 10 km catchment,
#           and jurisdiction -- and ask which, if any, predicts exposure.
#      (ii) NONLINEARITY. Peng finds VMT responds to the ratio only outside
#           roughly 1.2-2.8 (jobs per HOUSEHOLD) and is flat between. Our ratio
#           is jobs per resident WORKER, so his thresholds do not transfer
#           directly (at ~1.5 workers per household his band is roughly 0.8-1.9
#           in our units); rather than import them, we test for curvature
#           empirically with a cubic in the ratio and a joint Wald test on the
#           quadratic and cubic terms, alongside a binned display of the shape.
#    A horse race against metropolitan position asks the question a reviewer
#    from this literature will ask directly: why not simply use jobs-housing
#    balance as the predictor?
#
#  PART B -- THE LADDER ON THE LITERATURE'S OWN OUTCOMES (Giuliano & Small
#    1993). Our specification ladder explains workplace-location segregation.
#    Here the identical naive -> position-adjusted comparison is run with the
#    dependent variables this literature actually uses -- mean commute
#    distance and destination portfolio breadth -- for both the exclusionary
#    zoning measure and the jobs-housing ratio. If those associations attenuate
#    as ours does, the decomposition is not peculiar to our outcome. NOTE the
#    specification differs from the main text by design: the travel literature
#    estimates a MAIN EFFECT of land use on commuting, so these are main
#    effects, not the residential-segregation interactions of the headline.
#
# Descriptive/decomposition language throughout; no causal claims.
# Output: output/models/p4_jhbalance_scales.csv
#         output/models/p4_jhbalance_bins.csv
#         output/models/p4_literature_ladder.csv
# ==============================================================================

source("60_co_setup.R")
library(fixest)

zscore <- function(x) { x[!is.finite(x)] <- NA_real_; as.numeric(scale(x)) }
winz <- function(x, p = 0.01) {
  q <- quantile(x, c(p, 1 - p), na.rm = TRUE)
  pmin(pmax(x, q[1]), q[2])
}
tidy1 <- function(fit, id, keep) {
  ct <- as.data.frame(summary(fit)$coeftable)
  ct$term <- rownames(ct); rownames(ct) <- NULL
  names(ct)[1:4] <- c("estimate", "std.error", "statistic", "p.value")
  ct |> filter(term %in% keep) |> mutate(model_id = id, n_obs = fit$nobs)
}

## ---- frames ------------------------------------------------------------------
dat  <- readRDS(file.path(DIR_CO_OUT, "co_analysis_home_panel.rds"))
acc  <- readRDS(file.path(DIR_CO_CLEAN, "co_accessibility_2023.rds"))
cent <- readRDS(file.path(DIR_CLEAN, "tract_centroids_km.rds")) |>
  transmute(tract_id = as.character(GEOID), CBSA_Code, X_km, Y_km) |>
  filter(CBSA_Code == CO_CBSA_MAIN)
rac <- readRDS(file.path(DIR_CO_CLEAN,
                         sprintf("co_rac_tract_%d.rds", CO_ANCHOR_YEAR))) |>
  transmute(tract_id, workers = C000)
wac <- readRDS(file.path(DIR_CO_CLEAN,
                         sprintf("co_wac_tract_%d.rds", CO_ANCHOR_YEAR))) |>
  transmute(tract_id, jobs = C000)

# regional frame for catchments: ALL Denver-MSA tracts, so a catchment is not
# truncated at the edge of the zoning sample
reg <- cent |>
  left_join(rac, by = "tract_id") |>
  left_join(wac, by = "tract_id") |>
  mutate(across(c(workers, jobs), ~ replace_na(.x, 0)))
D <- as.matrix(dist(reg[, c("X_km", "Y_km")]))
catch_ratio <- function(r) {
  W <- (D <= r)
  as.numeric(W %*% reg$jobs) / pmax(as.numeric(W %*% reg$workers), 1)
}
reg$jhr_5km  <- catch_ratio(5)
reg$jhr_10km <- catch_ratio(10)
reg$jhr_tract <- reg$jobs / pmax(reg$workers, 1)
message(sprintf("Catchment ratios built for %d MSA tracts", nrow(reg)))

xs <- dat |>
  filter(year == CO_ANCHOR_YEAR, in_scope, denver_msa,
         n_commuters >= P3_MIN_COMMUTERS) |>
  left_join(acc, by = "tract_id") |>
  left_join(reg |> select(tract_id, jhr_tract, jhr_5km, jhr_10km),
            by = "tract_id")
# jurisdiction-scale ratio (frame tracts, as in 77/79)
jur <- xs |>
  left_join(rac, by = "tract_id") |> left_join(wac, by = "tract_id") |>
  mutate(across(c(workers, jobs), ~ replace_na(.x, 0))) |>
  group_by(jurisd_main) |>
  summarise(jhr_jurisd = sum(jobs) / pmax(sum(workers), 1), .groups = "drop")
xs <- xs |> left_join(jur, by = "jurisd_main")

SCALES <- c(tract = "jhr_tract", catch5 = "jhr_5km",
            catch10 = "jhr_10km", jurisdiction = "jhr_jurisd")
xs <- xs |>
  mutate(across(all_of(unname(SCALES)), winz)) |>
  mutate(across(c(wexp_whiteblack_wac_half, d_whiteblack_rac_half,
                  mean_dist_km, eff_n_dest, pct_reslow_of_res,
                  pct_black_rac, pct_lowincome_rac, log_worker_density_rac,
                  income_percapita_k, income_percapita_k_sq,
                  dist_cbd_km, dist_empctr_km, unname(SCALES)),
                zscore, .names = "z_{.col}"))

COVS <- paste(c("z_pct_black_rac", "z_pct_lowincome_rac",
                "z_log_worker_density_rac", "z_income_percapita_k",
                "z_income_percapita_k_sq"), collapse = " + ")
POS  <- "z_dist_cbd_km + z_dist_empctr_km"
Y    <- "z_wexp_whiteblack_wac_half"

## =============================================================================
## PART A -- balance measured at four scales; shape; horse race
## =============================================================================
# how close is each ratio to being a position measure? (our thesis, their terms)
cors <- sapply(SCALES, function(v)
  cor(xs[[paste0("z_", v)]], xs$z_dist_cbd_km, use = "complete.obs"))
message("corr(jobs-housing ratio, CBD distance): ",
        paste(sprintf("%s %.2f", names(cors), cors), collapse = " | "))

res <- list()
for (nm in names(SCALES)) {
  v <- paste0("z_", SCALES[[nm]])
  # (1) linear association with our outcome
  f1 <- feols(as.formula(sprintf("%s ~ %s + %s | county_fips", Y, v, COVS)),
              data = xs, cluster = ~jurisd_main)
  res[[length(res) + 1]] <- tidy1(f1, sprintf("A1_%s_linear", nm), v)
  # (2) curvature: cubic, joint Wald on the quadratic + cubic terms.
  #     Polynomial columns are precomputed rather than written as I(x^2) in
  #     the formula, so the coefficient names are predictable and the Wald
  #     test can index the vcov reliably.
  xs2 <- xs
  xs2$poly2 <- xs2[[v]]^2
  xs2$poly3 <- xs2[[v]]^3
  f2 <- feols(as.formula(sprintf("%s ~ %s + poly2 + poly3 + %s | county_fips",
                                 Y, v, COVS)),
              data = xs2, cluster = ~jurisd_main)
  nl <- intersect(c("poly2", "poly3"), names(coef(f2)))
  if (length(nl)) {
    b <- coef(f2)[nl]; V <- vcov(f2)[nl, nl, drop = FALSE]
    Wst <- tryCatch(as.numeric(t(b) %*% solve(V) %*% b),
                    error = function(e) NA_real_)
    pnl <- if (is.na(Wst)) NA_real_ else
      pchisq(Wst, df = length(nl), lower.tail = FALSE)
  } else {
    Wst <- NA_real_; pnl <- NA_real_
    message("  (", nm, ": polynomial terms dropped -- no curvature test)")
  }
  res[[length(res) + 1]] <- tibble(
    estimate = Wst, std.error = NA_real_, statistic = Wst, p.value = pnl,
    term = sprintf("joint: quadratic + cubic (%d df)", length(nl)),
    model_id = sprintf("A2_%s_nonlinearity", nm), n_obs = f2$nobs)
  # (3) horse race against metropolitan position
  f3 <- feols(as.formula(sprintf("%s ~ %s + %s + %s | county_fips",
                                 Y, v, POS, COVS)),
              data = xs, cluster = ~jurisd_main)
  res[[length(res) + 1]] <- tidy1(f3, sprintf("A3_%s_vs_position", nm),
                                  c(v, "z_dist_cbd_km", "z_dist_empctr_km"))
}
partA <- bind_rows(res) |>
  mutate(corr_with_cbd_dist = cors[sub("^A[0-9]_([a-z0-9]+)_.*$", "\\1", model_id)])
write.csv(partA, file.path(DIR_CO_MOD, "p4_jhbalance_scales.csv"),
          row.names = FALSE)
print(as.data.frame(partA |> select(model_id, term, estimate, std.error, p.value)),
      digits = 3)

# the shape, displayed rather than modeled (Peng's question, our outcome)
bins <- xs |>
  filter(!is.na(jhr_10km)) |>
  mutate(bin = ntile(jhr_10km, 5)) |>
  group_by(bin) |>
  summarise(n = n(),
            jh_ratio_10km = mean(jhr_10km),
            wexp = mean(wexp_whiteblack_wac_half, na.rm = TRUE),
            res_seg = mean(d_whiteblack_rac_half, na.rm = TRUE),
            commute_km = mean(mean_dist_km, na.rm = TRUE),
            dist_cbd_km = mean(dist_cbd_km, na.rm = TRUE), .groups = "drop")
write.csv(bins, file.path(DIR_CO_MOD, "p4_jhbalance_bins.csv"),
          row.names = FALSE)
print(as.data.frame(bins), digits = 3)

## =============================================================================
## PART B -- the ladder on the literature's own outcomes (main effects)
## =============================================================================
OUTCOMES <- c(commute_distance = "z_mean_dist_km",
              destination_breadth = "z_eff_n_dest",
              workplace_segregation = Y)
PREDS <- c(exclusionary = "z_pct_reslow_of_res",
           jobs_housing_ratio = "z_jhr_10km")
res2 <- list()
for (on in names(OUTCOMES)) for (pn in names(PREDS)) {
  y <- OUTCOMES[[on]]; x <- PREDS[[pn]]
  f_n <- feols(as.formula(sprintf("%s ~ %s + %s | county_fips", y, x, COVS)),
               data = xs, cluster = ~jurisd_main)
  f_p <- feols(as.formula(sprintf("%s ~ %s + %s + %s | county_fips",
                                  y, x, POS, COVS)),
               data = xs, cluster = ~jurisd_main)
  res2[[length(res2) + 1]] <- tidy1(f_n, sprintf("B_%s_%s_naive", on, pn), x)
  res2[[length(res2) + 1]] <- tidy1(f_p, sprintf("B_%s_%s_position", on, pn), x)
}
partB <- bind_rows(res2)
partB$pct_of_naive <- NA_real_
for (on in names(OUTCOMES)) for (pn in names(PREDS)) {
  i_n <- partB$model_id == sprintf("B_%s_%s_naive", on, pn)
  i_p <- partB$model_id == sprintf("B_%s_%s_position", on, pn)
  if (any(i_n) && any(i_p))
    partB$pct_of_naive[i_p] <- 100 * partB$estimate[i_p] / partB$estimate[i_n]
}
write.csv(partB, file.path(DIR_CO_MOD, "p4_literature_ladder.csv"),
          row.names = FALSE)
print(as.data.frame(partB |> select(model_id, estimate, std.error, p.value,
                                    pct_of_naive)), digits = 3)
message("85 complete.")
