# ==============================================================================
# 88_co_manuscript_statistics.R
# Small descriptive statistics quoted in the text that no other script writes
# to a file. Each is computed on the same 2023 Denver model frame as 68.
#
#   tract geometry   circle-equivalent tract diameter: median and the share of
#                    tracts wider than 5, 10 and 20 km (how often the index's
#                    distance truncation can bind)
#   tract size       correlation of log land area with exclusionary zoning and
#                    with CBD distance
#   levels ladder    the county-fixed-effects-only association of exclusionary
#                    zoning with residential segregation, refit on the
#                    complete-case sample the covariate-adjusted rows use
#   retention        correlation of residential restrictiveness with resident
#                    retention across jurisdictions (from 86's jurisdiction file)
#
# Output: output/models/p4_manuscript_statistics.csv   (statistic, value, n)
# ==============================================================================

source("60_co_setup.R")
library(fixest)

acc <- readRDS(file.path(DIR_CO_CLEAN, "co_accessibility_2023.rds"))
dat <- readRDS(file.path(DIR_CO_OUT, "co_analysis_home_panel.rds"))
zscore <- function(x) { x[!is.finite(x)] <- NA_real_; as.numeric(scale(x)) }
xs <- dat |>
  filter(year == CO_ANCHOR_YEAR, in_scope, denver_msa,
         n_commuters >= P3_MIN_COMMUTERS) |>
  left_join(acc |> select(tract_id, dist_cbd_km), by = "tract_id")
stopifnot("aland_km2" %in% names(xs))

out <- list()
add <- function(statistic, value, n = NA_integer_)
  out[[length(out) + 1]] <<- tibble(statistic = statistic, value = value, n = n)

## tract geometry
diam <- 2 * sqrt(xs$aland_km2 / pi)
add("tract_diameter_km_median", median(diam, na.rm = TRUE), sum(!is.na(diam)))
for (k in c(5, 10, 20))
  add(sprintf("pct_tracts_diameter_gt_%dkm", k),
      100 * mean(diam > k, na.rm = TRUE), sum(!is.na(diam)))

## tract size
la <- log(xs$aland_km2)
add("cor_log_area_exclusionary", cor(la, xs$pct_reslow_of_res, use = "complete.obs"),
    sum(is.finite(la) & is.finite(xs$pct_reslow_of_res)))
add("cor_log_area_dist_cbd",     cor(la, xs$dist_cbd_km,       use = "complete.obs"),
    sum(is.finite(la) & is.finite(xs$dist_cbd_km)))

## levels ladder, baseline on the complete-case sample
covs <- c("pct_black_rac", "pct_lowincome_rac", "log_worker_density_rac",
          "income_percapita_k")
cc <- xs |> filter(if_all(all_of(c("d_whiteblack_rac_half", "pct_reslow_of_res",
                                   covs)), ~ is.finite(.x))) |>
  mutate(z_y = zscore(d_whiteblack_rac_half), z_x = zscore(pct_reslow_of_res))
fit <- feols(z_y ~ z_x | county_fips, data = cc, cluster = ~jurisd_main)
ct <- summary(fit)$coeftable
add("levels_baseline_common_sample_estimate", ct["z_x", 1], fit$nobs)
add("levels_baseline_common_sample_se",       ct["z_x", 2], fit$nobs)
add("levels_baseline_common_sample_p",        ct["z_x", 4], fit$nobs)

## retention vs restrictiveness across jurisdictions
jf <- file.path(DIR_CO_MOD, "p4_cervero_jurisdiction.csv")
if (file.exists(jf)) {
  j <- read.csv(jf)
  add("cor_restrictiveness_retention_jurisdictions",
      cor(j$mean_pct_reslow, j$pct_residents_retained, use = "complete.obs"),
      nrow(j))
}

out <- bind_rows(out)
write.csv(out, file.path(DIR_CO_MOD, "p4_manuscript_statistics.csv"), row.names = FALSE)
print(as.data.frame(out |> mutate(value = round(value, 4))))
message("88 complete.")
