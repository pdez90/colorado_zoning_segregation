# ==============================================================================
# 73_co_spatial_inference.R      [PAPER 4, step 13]
# Spatial inference for the preferred specification. Every major variable in
# this paper is spatial, so residual spatial autocorrelation is expected;
# the question is whether it changes the inference.
#
#   (a) Moran's I on the residuals of the preferred model, computed manually
#       (no spdep dependency) under two weight schemes: k = 8 nearest
#       neighbours and inverse distance within 10 km; analytic z under
#       randomization.
#   (b) Conley spatial-HAC standard errors (Bartlett kernel) at 10 / 20 /
#       50 km cutoffs for every coefficient of the preferred model, computed
#       manually on the projected km coordinates: the textbook
#       sandwich X'KX with K_ij = max(1 - d_ij / cutoff, 0), no finite-sample
#       correction. Written out so the kernel and distances are explicit.
#
# Residual spatial autocorrelation is reported alongside the Conley SEs;
# Conley p-values use the same t reference (G - 1 df) as the clustered row.
#
# Output: output/models/p4_spatial_inference.csv
# ==============================================================================

source("60_co_setup.R")
library(fixest)

acc  <- readRDS(file.path(DIR_CO_CLEAN, "co_accessibility_2023.rds"))
dat  <- readRDS(file.path(DIR_CO_OUT, "co_analysis_home_panel.rds"))
cent <- readRDS(file.path(DIR_CLEAN, "tract_centroids_km.rds")) |>
  transmute(tract_id = as.character(GEOID), X_km, Y_km)

zscore <- function(x) { x[!is.finite(x)] <- NA_real_; as.numeric(scale(x)) }
xs <- dat |>
  filter(year == CO_ANCHOR_YEAR, in_scope, denver_msa,
         n_commuters >= P3_MIN_COMMUTERS) |>
  left_join(acc,  by = "tract_id") |>
  left_join(cent, by = "tract_id") |>
  mutate(across(c(wexp_whiteblack_wac_half, d_whiteblack_rac_half,
                  pct_reslow_of_res, pct_black_rac, pct_lowincome_rac,
                  log_worker_density_rac, income_percapita_k,
                  income_percapita_k_sq, dist_cbd_km, dist_empctr_km),
                zscore, .names = "z_{.col}"))

X   <- "z_d_whiteblack_rac_half"
COVS <- paste(c("z_pct_black_rac", "z_pct_lowincome_rac",
                "z_log_worker_density_rac", "z_income_percapita_k",
                "z_income_percapita_k_sq"), collapse = " + ")
POS <- sprintf(paste("z_dist_cbd_km + z_dist_empctr_km +",
                     "%s:z_dist_cbd_km + %s:z_dist_empctr_km"), X, X)
fml <- sprintf(
  "z_wexp_whiteblack_wac_half ~ %s * z_pct_reslow_of_res + %s + %s | county_fips",
  X, COVS, POS)
fit <- feols(as.formula(fml), data = xs, cluster = ~jurisd_main)
used <- obs(fit)
d    <- xs[used, ]
e    <- resid(fit)
n    <- length(e)
D    <- as.matrix(dist(d[, c("X_km", "Y_km")]))

## ---- (a) Moran's I -----------------------------------------------------------
morans <- function(W) {
  diag(W) <- 0
  rs <- rowSums(W); W <- W / ifelse(rs > 0, rs, 1)
  z <- e - mean(e); S0 <- sum(W)
  I  <- (n / S0) * sum(z * (W %*% z)) / sum(z^2)
  EI <- -1 / (n - 1)
  S1 <- 0.5 * sum((W + t(W))^2); S2 <- sum((rowSums(W) + colSums(W))^2)
  b2 <- n * sum(z^4) / (sum(z^2)^2)
  A  <- n * ((n^2 - 3 * n + 3) * S1 - n * S2 + 3 * S0^2)
  B  <- b2 * ((n^2 - n) * S1 - 2 * n * S2 + 6 * S0^2)
  VI <- (A - B) / ((n - 1) * (n - 2) * (n - 3) * S0^2) - EI^2
  zz <- (I - EI) / sqrt(VI)
  c(I = I, z = zz, p = 2 * pnorm(-abs(zz)))
}
W8 <- matrix(0, n, n)
for (i in seq_len(n)) W8[i, order(D[i, ])[2:9]] <- 1
Wd <- ifelse(D > 0 & D <= 10, 1 / pmax(D, 0.1), 0)
m1 <- morans(W8); m2 <- morans(Wd)

## ---- (b) Conley spatial-HAC SEs ---------------------------------------------
# rebuild the exact design matrix (county FE as dummies, same rows)
mm <- model.matrix(
  ~ z_d_whiteblack_rac_half * z_pct_reslow_of_res + z_pct_black_rac +
    z_pct_lowincome_rac + z_log_worker_density_rac + z_income_percapita_k +
    z_income_percapita_k_sq + z_dist_cbd_km + z_dist_empctr_km +
    z_d_whiteblack_rac_half:z_dist_cbd_km +
    z_d_whiteblack_rac_half:z_dist_empctr_km + factor(county_fips), data = d)
bols <- solve(crossprod(mm), crossprod(mm, d$z_wexp_whiteblack_wac_half))
eo   <- d$z_wexp_whiteblack_wac_half - mm %*% bols   # OLS residuals
XtXi <- solve(crossprod(mm))
conley_se <- function(cutoff) {
  K <- pmax(1 - D / cutoff, 0); diag(K) <- 1
  Xe <- mm * as.numeric(eo)
  V  <- XtXi %*% (t(Xe) %*% K %*% Xe) %*% XtXi
  sqrt(diag(V))
}
tt <- "z_d_whiteblack_rac_half:z_pct_reslow_of_res"
b  <- bols[tt, 1]
n_clu <- dplyr::n_distinct(d$jurisd_main)
rows <- list(
  tibble(stat = "morans_I_knn8", value = m1["I"], z = m1["z"], p = m1["p"]),
  tibble(stat = "morans_I_invdist10km", value = m2["I"], z = m2["z"],
         p = m2["p"]),
  tibble(stat = "interaction_jurisdiction_cluster",
         value = coef(fit)[tt], z = NA,
         p = summary(fit)$coeftable[tt, 4]))
for (ck in c(10, 20, 50)) {
  se <- conley_se(ck)[tt]
  rows[[length(rows) + 1]] <- tibble(
    stat = sprintf("interaction_conley_%dkm", ck), value = b,
    std.error = se, z = b / se,
    # same reference distribution as the clustered row above (t, G - 1 df),
    # so the two are comparable; the asymptotic-normal p is kept alongside
    p = 2 * pt(-abs(b / se), df = n_clu - 1),
    p_normal = 2 * pnorm(-abs(b / se)))
}
# NOTE: Conley SEs need not increase monotonically with bandwidth -- spatial
# covariance can contain positive and negative components, so a wide kernel
# can average them down. The 10-20 km entries are the informative ones.

## ---- (c) cross-check the manual estimator against fixest::vcov_conley --------
# The manual estimator uses a Bartlett kernel on PROJECTED km distances with
# no finite-sample correction; fixest::vcov_conley is a different
# implementation (spherical distances, its own kernel and adjustments) and
# returns a SMALLER standard error here. The manuscript reports the manual,
# more conservative estimate; the fixest value is saved alongside it for
# transparency. Older fixest versions lack vcov_conley, hence the tryCatch.
chk <- tryCatch({
  ll <- readRDS(file.path(DIR_CLEAN, "tract_centroids_km.rds")) |>
    transmute(tract_id = as.character(GEOID), X_km, Y_km) |>
    filter(tract_id %in% d$tract_id)
  pts <- st_as_sf(ll |> mutate(x_m = X_km * 1000, y_m = Y_km * 1000),
                  coords = c("x_m", "y_m"), crs = CRS_METERS) |>
    st_transform(4326)
  co <- st_coordinates(pts)
  dl <- d |> left_join(tibble(tract_id = ll$tract_id,
                              lon = co[, 1], lat = co[, 2]), by = "tract_id")
  fit_c <- feols(as.formula(fml), data = dl,
                 vcov = conley(20, distance = "spherical"))
  se_fx <- summary(fit_c)$coeftable[tt, 2]
  se_mn <- conley_se(20)[tt]
  message(sprintf(
    "Conley cross-check @20km: manual SE = %.4f | fixest SE = %.4f | ratio %.3f",
    se_mn, se_fx, se_mn / se_fx))
  tibble(stat = "conley20_fixest_crosscheck_se", value = se_fx,
         std.error = se_mn, z = NA, p = NA)
}, error = function(e) {
  message("fixest vcov_conley cross-check skipped: ", conditionMessage(e))
  NULL
})
if (!is.null(chk)) rows[[length(rows) + 1]] <- chk
out <- bind_rows(rows)
write.csv(out, file.path(DIR_CO_MOD, "p4_spatial_inference.csv"),
          row.names = FALSE)
print(as.data.frame(out))
message("73 complete.")
