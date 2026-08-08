# ==============================================================================
# 69_co_robustness.R      [PAPER 4, step 9]
# The three remaining referee-proofing analyses for the JAPA paper.
#
#  BLOCK L  (levels ladder) Does zoning predict the LEVEL of residential
#           segregation once you condition step by step? Sequential buildup:
#             M0 raw -> M1 +demographics -> M2 +income -> M3 +density
#             -> M4 +metropolitan position
#           Verified: pct_reslow_of_res goes +0.175 (p=.029) -> +0.124 ->
#           +0.108 -> -0.049 (ns) -> -0.053 (ns). DENSITY is what absorbs it;
#           position adds nothing further. Report the whole path in the paper
#           rather than only the final column: the raw association exists and
#           replicates the zoning-segregation literature, and the point is
#           what happens to it, not that it was never there.
#           CAVEAT to state in the text: density is plausibly a MEDIATOR of
#           zoning (zoning caps density), so M3-M4 may over-control; M0-M2 and
#           M3-M4 bracket the truth.
#
#  BLOCK A  (ADU validity) What is pct_adu_res actually measuring in 2023?
#           Verified: NOT central early adopters. Correlation with distance to
#           CBD is only +0.09, and ADU-permissive tracts sit FARTHER out
#           (18.2 km vs 15.5 km). It is a jurisdictional patchwork -- Lakewood
#           99.6%, Douglas Co. uninc. 98.3%, JeffCo 96.3%, Aurora 93.8% vs
#           Centennial 3.7%, Westminster 9.0%, DENVER only 53.7%. So the
#           variable marks suburban large-lot jurisdictions where an accessory
#           unit is trivially accommodated, NOT ADU liberalization in the sense
#           the reform literature means. Demoted to a measurement caveat.
#           NOTE the timing: Colorado's HB24-1152 did not require ADU
#           permission until June 30 2025, so the Oct 2023 snapshot is a
#           genuine PRE-MANDATE baseline -- a real before/after becomes
#           possible once a post-2025 zoning vintage exists.
#
#  BLOCK T  (effect translation) Everything above is in SD units, which no
#           practitioner can act on. Converts the preferred (rung B) estimates
#           into native units, per-tercile contrasts, and simple slopes.
#
# Needs: 68 (accessibility), 63's home panel.
# Output: output/models/p4_levels_ladder.csv, p4_adu_validity.csv,
#         p4_effect_translation.csv, p4_tercile_descriptives.csv
# ==============================================================================

source("60_co_setup.R")
library(fixest)

acc <- readRDS(file.path(DIR_CO_CLEAN, "co_accessibility_2023.rds"))
dat <- readRDS(file.path(DIR_CO_OUT, "co_analysis_home_panel.rds"))
zscore <- function(x) { x[!is.finite(x)] <- NA_real_; as.numeric(scale(x)) }

MODEL_VARS <- c("wexp_whiteblack_wac_half", "d_whiteblack_rac_half",
                "pct_res_low", "pct_reslow_of_res", "pct_adu_res",
                "zoning_entropy", "pct_job_zone", "pct_black_rac",
                "pct_lowincome_rac", "log_worker_density_rac",
                "income_percapita_k", "income_percapita_k_sq",
                "dist_cbd_km", "dist_empctr_km", "log_jobs_grav",
                "jobs_housing_ratio", "mean_dist_km", "eff_n_dest")
xs <- dat |>
  filter(year == CO_ANCHOR_YEAR, in_scope, denver_msa,
         n_commuters >= P3_MIN_COMMUTERS) |>
  left_join(acc, by = "tract_id") |>
  mutate(across(any_of(MODEL_VARS), zscore, .names = "z_{.col}"))

tidy_fx <- function(fit, model_id, keep) {
  if (is.null(fit)) return(NULL)
  ct <- as.data.frame(summary(fit)$coeftable)
  ct$term <- rownames(ct); rownames(ct) <- NULL
  names(ct)[1:4] <- c("estimate", "std.error", "statistic", "p.value")
  ct |> filter(term %in% keep) |>
    mutate(model_id = model_id, n_obs = fit$nobs)
}
safe <- function(fml, model_id, keep, d = xs)
  tidy_fx(tryCatch(feols(as.formula(fml), data = d, cluster = ~jurisd_main),
                   error = function(e) NULL), model_id, keep)

## =============================================================================
## BLOCK L -- levels ladder
## =============================================================================
COV_FULL <- paste(c("z_pct_black_rac", "z_pct_lowincome_rac",
                    "z_log_worker_density_rac", "z_income_percapita_k",
                    "z_income_percapita_k_sq"), collapse = " + ")
ladder <- list(
  M0_raw      = "",
  M1_demog    = "z_pct_black_rac + z_pct_lowincome_rac",
  M2_income   = "z_pct_black_rac + z_pct_lowincome_rac + z_income_percapita_k + z_income_percapita_k_sq",
  M3_density  = COV_FULL,
  M4_position = paste(COV_FULL, "+ z_dist_cbd_km + z_dist_empctr_km"))

lev <- list()
for (zv in c("pct_reslow_of_res", "pct_res_low", "pct_job_zone",
             "zoning_entropy", "pct_adu_res"))
  for (s in names(ladder)) {
    rhs <- if (nzchar(ladder[[s]])) paste("z_", zv, " + ", ladder[[s]], sep = "")
           else paste0("z_", zv)
    lev[[length(lev) + 1]] <- safe(
      sprintf("z_d_whiteblack_rac_half ~ %s | county_fips", rhs),
      sprintf("%s__%s", zv, s), paste0("z_", zv))
  }
levels_ladder <- bind_rows(lev) |>
  mutate(zoning = sub("__.*", "", model_id), spec = sub(".*__", "", model_id))
write.csv(levels_ladder, file.path(DIR_CO_MOD, "p4_levels_ladder.csv"),
          row.names = FALSE)
message("Levels ladder:")
print(levels_ladder |> select(zoning, spec, estimate, p.value) |>
        pivot_wider(names_from = spec, values_from = c(estimate, p.value)) |>
        as.data.frame())

## =============================================================================
## BLOCK A -- what is the ADU variable measuring?
## =============================================================================
adu_cor <- xs |>
  summarise(across(c(dist_cbd_km, log_jobs_grav, log_worker_density_rac,
                     pct_reslow_of_res, income_percapita_k),
                   ~ cor(.x, pct_adu_res, use = "complete.obs"))) |>
  pivot_longer(everything(), names_to = "variable",
               values_to = "cor_with_pct_adu_res")
adu_jur <- xs |>
  group_by(jurisd_main) |>
  summarise(n_tracts = n(), mean_adu = mean(pct_adu_res, na.rm = TRUE),
            mean_dist_cbd = mean(dist_cbd_km, na.rm = TRUE), .groups = "drop") |>
  filter(n_tracts >= 8) |> arrange(desc(mean_adu))
adu_split <- xs |>
  mutate(adu_permissive = pct_adu_res >= 90) |>
  group_by(adu_permissive) |>
  summarise(n = n(), mean_dist_cbd = mean(dist_cbd_km, na.rm = TRUE),
            mean_reslow = mean(pct_reslow_of_res, na.rm = TRUE),
            .groups = "drop")
write.csv(adu_cor, file.path(DIR_CO_MOD, "p4_adu_validity.csv"),
          row.names = FALSE)
write_codiag(adu_jur, "69_adu_by_jurisdiction")
write_codiag(adu_split, "69_adu_permissive_split")
print(as.data.frame(adu_jur)); print(as.data.frame(adu_split))

## =============================================================================
## BLOCK T -- effect translation (preferred rung-B specification)
## =============================================================================
X   <- "z_d_whiteblack_rac_half"
POS <- sprintf(paste("z_dist_cbd_km + z_dist_empctr_km +",
                     "%s:z_dist_cbd_km + %s:z_dist_empctr_km"), X, X)
sd_w <- sd(xs$wexp_whiteblack_wac_half, na.rm = TRUE)
sd_r <- sd(xs$d_whiteblack_rac_half,    na.rm = TRUE)
mean_w <- mean(xs$wexp_whiteblack_wac_half, na.rm = TRUE)

trans <- map_df(c("pct_reslow_of_res", "zoning_entropy", "pct_job_zone",
                  "pct_adu_res"), function(zv) {
  fit <- feols(as.formula(sprintf(
    "z_wexp_whiteblack_wac_half ~ %s * z_%s + %s + %s | county_fips",
    X, zv, COV_FULL, POS)), data = xs, cluster = ~jurisd_main)
  b_x <- coef(fit)[X]; b_i <- coef(fit)[sprintf("%s:z_%s", X, zv)]
  tibble(zoning = zv,
         interaction_sd = b_i,
         slope_lo_sd = b_x - b_i, slope_hi_sd = b_x + b_i,
         # native: change in wexp per 1-unit change in residential D
         slope_lo_native = (b_x - b_i) * sd_w / sd_r,
         slope_hi_native = (b_x + b_i) * sd_w / sd_r,
         # +1 SD of residential segregation, as % of mean wexp
         pct_of_mean_lo = 100 * (b_x - b_i) * sd_w / mean_w,
         pct_of_mean_hi = 100 * (b_x + b_i) * sd_w / mean_w)
})
write.csv(trans, file.path(DIR_CO_MOD, "p4_effect_translation.csv"),
          row.names = FALSE)
message(sprintf("SD wexp=%.4f | SD res_seg=%.4f | mean wexp=%.4f",
                sd_w, sd_r, mean_w))
print(as.data.frame(trans))

## ---- native-unit tercile descriptives (the practitioner-legible table) -------
terc <- xs |>
  mutate(tercile = ntile(pct_reslow_of_res, 3)) |>
  group_by(tercile) |>
  summarise(n_tracts = n(),
            across(c(pct_reslow_of_res, d_whiteblack_rac_half,
                     wexp_whiteblack_wac_half, jobs_housing_ratio,
                     mean_dist_km, eff_n_dest, dist_cbd_km, pct_black_rac,
                     income_percapita_k),
                   ~ mean(.x, na.rm = TRUE)), .groups = "drop") |>
  mutate(across(where(is.numeric), ~ round(.x, 3)))
write.csv(terc, file.path(DIR_CO_MOD, "p4_tercile_descriptives.csv"),
          row.names = FALSE)
print(as.data.frame(terc))

## ---- headline coupling correlation -------------------------------------------
r_couple <- cor(xs$d_whiteblack_rac_half, xs$wexp_whiteblack_wac_half,
                use = "complete.obs")
message(sprintf("Tract-level coupling r (Denver %d, in-scope frame): %.3f",
                CO_ANCHOR_YEAR, r_couple))
# NOTE: 62's diagnostic reports a HIGHER r because it covers all regional
# tracts in the 3 CBSAs, not this Denver in-scope model frame. Quote whichever
# matches the sample you are describing -- do not mix them.
write_codiag(tibble(sample = c("Denver in-scope model frame"),
                    year = CO_ANCHOR_YEAR, n = nrow(xs), r = r_couple),
             "69_coupling_correlation")
message("69 complete.")
