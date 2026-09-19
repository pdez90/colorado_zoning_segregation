# ==============================================================================
# 86_co_cervero_test.R      [PAPER 4, step 26 -- SI: testing Cervero's claim]
#
# Cervero (1989) attributes jobs-housing imbalance to fiscal and exclusionary
# zoning in EMPLOYMENT-RICH jurisdictions: places that host jobs restrict the
# housing that would let those workers live nearby, so the workforce commutes
# in. His cases are Bay Area employment hubs -- Palo Alto, Santa Clara,
# Sunnyvale -- with jobs-housing ratios above 2.5 alongside restrictive
# residential zoning. The mechanism requires a CONJUNCTION: job-richness AND
# housing restriction in the same jurisdiction.
#
# We have what Cervero did not: a harmonized regulatory map, so the conjunction
# can be measured directly rather than inferred from housing prices. This
# script asks whether the Denver region exhibits it.
#
#   PART A  Tract level: does exclusionary zoning predict the jobs-housing
#           ratio, before and after metropolitan position?
#   PART B  Jurisdiction level (Cervero's own unit): is job-richness
#           associated with residential restrictiveness, and among job-rich
#           jurisdictions does restrictiveness track worker importing?
#   PART C  The conjunction as a 2x2: how much of the region's employment and
#           housing sits in each quadrant of job-richness x restrictiveness?
#
# Interpretation discipline: a negative or absent conjunction is NOT evidence
# against Cervero's Bay Area findings. It is evidence that the regional
# configuration differs -- which is itself the point, since it locates where
# exclusion sits relative to employment in this region.
#
# Output: output/models/p4_cervero_tract.csv
#         output/models/p4_cervero_jurisdiction.csv
#         output/figures/p4_fig_cervero_quadrants.png
# ==============================================================================

source("60_co_setup.R")
library(fixest)
library(ggrepel)

zscore <- function(x) { x[!is.finite(x)] <- NA_real_; as.numeric(scale(x)) }
tidy1 <- function(fit, id, keep) {
  ct <- as.data.frame(summary(fit)$coeftable)
  ct$term <- rownames(ct); rownames(ct) <- NULL
  names(ct)[1:4] <- c("estimate", "std.error", "statistic", "p.value")
  ct |> filter(term %in% keep) |> mutate(model_id = id, n_obs = fit$nobs)
}

dat <- readRDS(file.path(DIR_CO_OUT, "co_analysis_home_panel.rds"))
acc <- readRDS(file.path(DIR_CO_CLEAN, "co_accessibility_2023.rds"))
rac <- readRDS(file.path(DIR_CO_CLEAN,
                         sprintf("co_rac_tract_%d.rds", CO_ANCHOR_YEAR))) |>
  transmute(tract_id, workers = C000)
wac <- readRDS(file.path(DIR_CO_CLEAN,
                         sprintf("co_wac_tract_%d.rds", CO_ANCHOR_YEAR))) |>
  transmute(tract_id, jobs = C000)

xs <- dat |>
  filter(year == CO_ANCHOR_YEAR, in_scope, denver_msa,
         n_commuters >= P3_MIN_COMMUTERS) |>
  left_join(acc, by = "tract_id") |>
  left_join(rac, by = "tract_id") |> left_join(wac, by = "tract_id") |>
  mutate(across(c(workers, jobs), ~ replace_na(.x, 0)),
         jhr = jobs / pmax(workers, 1),
         log_jhr = log(pmax(jhr, 0.01))) |>
  mutate(across(c(log_jhr, pct_reslow_of_res, pct_job_zone,
                  pct_black_rac, pct_lowincome_rac, log_worker_density_rac,
                  income_percapita_k, income_percapita_k_sq,
                  dist_cbd_km, dist_empctr_km), zscore, .names = "z_{.col}"))

COVS <- paste(c("z_pct_black_rac", "z_pct_lowincome_rac",
                "z_income_percapita_k", "z_income_percapita_k_sq"),
              collapse = " + ")   # NOTE: worker density deliberately omitted --
# it is close to mechanically related to the jobs-housing ratio's denominator
POS <- "z_dist_cbd_km + z_dist_empctr_km"

## ---- PART A: tract level -----------------------------------------------------
res <- list()
for (zv in c("z_pct_reslow_of_res", "z_pct_job_zone")) {
  res[[length(res)+1]] <- tidy1(feols(as.formula(sprintf(
    "z_log_jhr ~ %s + %s | county_fips", zv, COVS)),
    data = xs, cluster = ~jurisd_main), paste0("A_", zv, "_naive"), zv)
  res[[length(res)+1]] <- tidy1(feols(as.formula(sprintf(
    "z_log_jhr ~ %s + %s + %s | county_fips", zv, POS, COVS)),
    data = xs, cluster = ~jurisd_main), paste0("A_", zv, "_position"), zv)
}
partA <- bind_rows(res)
write.csv(partA, file.path(DIR_CO_MOD, "p4_cervero_tract.csv"), row.names = FALSE)
print(as.data.frame(partA |> select(model_id, term, estimate, std.error, p.value)),
      digits = 3)

## ---- PART B: jurisdiction level (Cervero's unit) -----------------------------
jf <- file.path(DIR_CO_MOD, "p4_selfcontainment_jurisdiction.csv")
stopifnot(file.exists(jf))
jur <- read.csv(jf) |>
  filter(workers_housed >= 2000) |>
  mutate(pct_jobs_imported = 100 - pct_jobs_filled_locally)

cp <- cor(jur$jobs_per_worker, jur$mean_pct_reslow)
cs <- cor(jur$jobs_per_worker, jur$mean_pct_reslow, method = "spearman")
message(sprintf(
  "CONJUNCTION TEST (n=%d jurisdictions): corr(jobs per worker, %% residential land low-density) = %.2f (Pearson), %.2f (Spearman)",
  nrow(jur), cp, cs))
jr <- jur |> filter(jobs_per_worker >= 1)
message(sprintf(
  "Among the %d job-rich jurisdictions: corr(restrictiveness, %% jobs imported) = %.2f; import share ranges %.0f-%.0f%%",
  nrow(jr), cor(jr$mean_pct_reslow, jr$pct_jobs_imported),
  min(jr$pct_jobs_imported), max(jr$pct_jobs_imported)))

fit <- lm(log(jobs_per_worker) ~ mean_pct_reslow, data = jur,
          weights = workers_housed)
message(sprintf("  weighted slope: %.4f log-points per pp res-low (p = %.3f)",
                coef(fit)[2], summary(fit)$coefficients[2, 4]))
co_write_stats("86_cervero_jurisdiction",
  n_jurisdictions = nrow(jur), cor_pearson = cp, cor_spearman = cs,
  n_job_rich = nrow(jr),
  cor_restrictive_imported_jobrich = cor(jr$mean_pct_reslow, jr$pct_jobs_imported),
  weighted_slope = coef(fit)[[2]], weighted_slope_p = summary(fit)$coefficients[2, 4])

## ---- PART C: the conjunction as a 2x2 ----------------------------------------
med <- median(jur$mean_pct_reslow)
jur <- jur |>
  mutate(quadrant = paste(ifelse(jobs_per_worker >= 1, "job-rich", "job-poor"),
                          ifelse(mean_pct_reslow >= med, "restrictive",
                                 "permissive"), sep = ", "))
quad <- jur |>
  group_by(quadrant) |>
  summarise(n = n(),
            mean_pct_reslow = weighted.mean(mean_pct_reslow, workers_housed),
            mean_jobs_per_worker = sum(jobs_hosted) / sum(workers_housed),
            pct_of_region_jobs = 100 * sum(jobs_hosted) / sum(jur$jobs_hosted),
            pct_of_region_workers = 100 * sum(workers_housed) /
                                          sum(jur$workers_housed),
            mean_pct_jobs_imported = weighted.mean(pct_jobs_imported,
                                                   jobs_hosted),
            .groups = "drop")
write.csv(jur |> select(jurisd_main, jobs_per_worker, mean_pct_reslow,
                        pct_jobs_imported, pct_residents_retained, quadrant),
          file.path(DIR_CO_MOD, "p4_cervero_jurisdiction.csv"), row.names = FALSE)
print(as.data.frame(quad), digits = 3)
write.csv(quad, file.path(DIR_CO_MOD, "p4_cervero_quadrants.csv"), row.names = FALSE)

pC <- ggplot(jur, aes(mean_pct_reslow, jobs_per_worker)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey55", linewidth = .4) +
  geom_vline(xintercept = med, linetype = 2, colour = "grey55", linewidth = .4) +
  geom_point(aes(size = jobs_hosted), colour = "#08519c", alpha = .75) +
  geom_text_repel(aes(label = jurisd_main), size = 2.4, max.overlaps = 14,
                  seed = 4) +
  scale_y_log10() +
  scale_size_area(max_size = 10, name = "jobs hosted", labels = scales::comma) +
  annotate("text", x = max(jur$mean_pct_reslow), y = max(jur$jobs_per_worker),
           hjust = 1, vjust = 1, size = 3, fontface = "italic", colour = "grey35",
           label = "Cervero's archetype:\njob-rich and housing-restrictive") +
  labs(title = "Where does exclusionary zoning sit relative to employment?",
       subtitle = paste("Denver-region jurisdictions with at least 2,000",
                        "resident workers. Cervero's archetype is the",
                        "upper-right quadrant:\nfive jurisdictions holding 13%",
                        "of regional jobs but housing 8% of its workers, and",
                        "importing 93% of their workforce."),
       x = "% of residential land zoned low-density (jurisdiction mean)",
       y = "local jobs per resident worker (log scale)") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")
ggsave(file.path(DIR_CO_FIG, "p4_fig_cervero_quadrants.png"), pC,
       width = 8.0, height = 6.6, dpi = 350, bg = "white")
message("86 complete.")
