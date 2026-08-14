# ==============================================================================
# 87_co_workbased.R      [PAPER 4, step 27 -- SI: the work-based direction]
#
# The paper so far is HOME-BASED: from where people live, how segregated are
# the workplaces they reach? The jobs-housing literature insists on the mirror
# question -- Zheng et al. (2021) estimate home-based AND work-based
# jobs-housing balance and find the built environment relates differently to
# each -- but poses it for commuting distance, never for segregation. This
# script supplies the work-based direction for our outcome, reusing the
# construct script 62 already computes and the paper has not used:
#
#     wres_i = flow-weighted mean RESIDENTIAL segregation of the neighborhoods
#              that supply workers to workplace tract i
#
#   PART A  The ladder, mirrored. Does the zoning AROUND A WORKPLACE moderate
#           the link between a workplace's own segregation and the segregation
#           of the neighborhoods it draws from -- and does metropolitan
#           position absorb it, as it does on the home side? The residential
#           restrictiveness of a work tract is the Cervero-relevant measure
#           here: it asks whether workplaces embedded in restrictively zoned
#           residential land draw from different neighborhoods.
#
#   PART B  The menu, mirrored. For each workplace, compare the segregation of
#           the ACCESSIBLE HOUSING (gravity-weighted over resident workers)
#           with the segregation of its REALIZED labor shed. On the home side
#           92% of the variance in realized exposure tracked the accessible
#           set; this asks whether the same holds looking outward from jobs.
#
#   PART C  The Cervero link at workplace scale: among job-rich work tracts,
#           does surrounding residential restrictiveness track longer
#           in-commutes and more segregated origins?
#
# Output: output/models/p4_workbased_ladder.csv
#         output/models/p4_workbased_menu.csv
#         output/figures/p4_fig_workbased_menu.png
# ==============================================================================

source("60_co_setup.R")
library(fixest)

zscore <- function(x) { x[!is.finite(x)] <- NA_real_; as.numeric(scale(x)) }
tidy1 <- function(fit, id, keep) {
  ct <- as.data.frame(summary(fit)$coeftable)
  ct$term <- rownames(ct); rownames(ct) <- NULL
  names(ct)[1:4] <- c("estimate", "std.error", "statistic", "p.value")
  ct |> filter(term %in% keep) |> mutate(model_id = id, n_obs = fit$nobs)
}

datw <- readRDS(file.path(DIR_CO_OUT, "co_analysis_work_panel.rds"))
acc  <- readRDS(file.path(DIR_CO_CLEAN, "co_accessibility_2023.rds"))
cent <- readRDS(file.path(DIR_CLEAN, "tract_centroids_km.rds")) |>
  transmute(tract_id = as.character(GEOID), CBSA_Code, X_km, Y_km) |>
  filter(CBSA_Code == CO_CBSA_MAIN)
seg23 <- readRDS(file.path(DIR_P2_OUT, "p2_tract_segregation_panel.rds")) |>
  filter(year == CO_ANCHOR_YEAR) |>
  transmute(tract_id, d_rac = d_whiteblack_rac_half)
rac <- readRDS(file.path(DIR_CO_CLEAN,
                         sprintf("co_rac_tract_%d.rds", CO_ANCHOR_YEAR))) |>
  transmute(tract_id, workers = C000)

xw <- datw |>
  filter(year == CO_ANCHOR_YEAR, in_scope, denver_msa,
         n_workers_od >= P3_MIN_COMMUTERS) |>
  left_join(acc, by = "tract_id")
stopifnot("wres_whiteblack_rac_half" %in% names(xw))

WCOV <- c("pct_black_wac", "pct_lowincome_wac", "log_worker_density_wac")
WCOV <- intersect(WCOV, names(xw))
xw <- xw |>
  mutate(across(any_of(c("wres_whiteblack_rac_half", "d_whiteblack_wac_half",
                         "mean_dist_km_work", "pct_reslow_of_res",
                         "pct_job_zone", "zoning_entropy", WCOV,
                         "dist_cbd_km", "dist_empctr_km")),
                zscore, .names = "z_{.col}"))
COVS <- paste(paste0("z_", WCOV), collapse = " + ")
Y <- "z_wres_whiteblack_rac_half"
X <- "z_d_whiteblack_wac_half"
POS <- sprintf(paste("z_dist_cbd_km + z_dist_empctr_km +",
                     "%s:z_dist_cbd_km + %s:z_dist_empctr_km"), X, X)

## ---- PART A: the ladder, mirrored --------------------------------------------
res <- list()
for (zv in c("pct_reslow_of_res", "pct_job_zone", "zoning_entropy")) {
  v <- paste0("z_", zv); keep <- sprintf("%s:%s", X, v)
  res[[length(res)+1]] <- tidy1(feols(as.formula(sprintf(
    "%s ~ %s * %s + %s | county_fips", Y, X, v, COVS)),
    data = xw, cluster = ~jurisd_main), sprintf("W_%s_naive", zv), keep)
  res[[length(res)+1]] <- tidy1(feols(as.formula(sprintf(
    "%s ~ %s * %s + %s + %s | county_fips", Y, X, v, COVS, POS)),
    data = xw, cluster = ~jurisd_main), sprintf("W_%s_position", zv), keep)
}
ladder <- bind_rows(res)
ladder$pct_of_naive <- NA_real_
for (zv in c("pct_reslow_of_res", "pct_job_zone", "zoning_entropy")) {
  i_n <- ladder$model_id == sprintf("W_%s_naive", zv)
  i_p <- ladder$model_id == sprintf("W_%s_position", zv)
  if (any(i_n) && any(i_p))
    ladder$pct_of_naive[i_p] <- 100 * ladder$estimate[i_p] / ladder$estimate[i_n]
}
write.csv(ladder, file.path(DIR_CO_MOD, "p4_workbased_ladder.csv"),
          row.names = FALSE)
print(as.data.frame(ladder |> select(model_id, estimate, std.error, p.value,
                                     pct_of_naive)), digits = 3)

## ---- PART B: the menu, mirrored ----------------------------------------------
reg <- cent |>
  inner_join(rac, by = "tract_id") |>
  inner_join(seg23, by = "tract_id") |>
  filter(!is.na(d_rac))
D  <- as.matrix(dist(as.matrix(reg[, c("X_km", "Y_km")])))
W  <- exp(-0.10 * D) * matrix(reg$workers, nrow(reg), nrow(reg), byrow = TRUE)
accR <- as.numeric((W %*% reg$d_rac) / rowSums(W))
opp <- tibble(tract_id = reg$tract_id, accR = accR)

xb <- xw |>
  left_join(opp, by = "tract_id") |>
  mutate(sorting_gap = wres_whiteblack_rac_half - accR) |>
  filter(!is.na(accR), !is.na(sorting_gap))
message(sprintf(
  "WORK-BASED var(wres) split: accessible-housing %.0f%% | sorting %.0f%% | r(accR,gap) = %.2f",
  100 * cov(xb$wres_whiteblack_rac_half, xb$accR) /
        var(xb$wres_whiteblack_rac_half),
  100 * cov(xb$wres_whiteblack_rac_half, xb$sorting_gap) /
        var(xb$wres_whiteblack_rac_half),
  cor(xb$accR, xb$sorting_gap)))

menu <- xb |>
  mutate(dec = ntile(d_whiteblack_wac_half, 10)) |>
  group_by(dec) |>
  summarise(n = n(),
            realized_labor_shed = mean(wres_whiteblack_rac_half, na.rm = TRUE),
            accessible_housing = mean(accR, na.rm = TRUE),
            gap = mean(sorting_gap, na.rm = TRUE), .groups = "drop")
write.csv(menu, file.path(DIR_CO_MOD, "p4_workbased_menu.csv"),
          row.names = FALSE)
print(as.data.frame(menu), digits = 3)

pM <- ggplot(xb, aes(accR, wres_whiteblack_rac_half)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = .4,
              colour = "grey50") +
  geom_point(size = 1.1, alpha = .55, colour = "#08519c") +
  labs(title = "Looking outward from the workplace: accessible housing vs realized labor shed",
       subtitle = paste("Each point is a workplace tract, 2023. Horizontal:",
                        "segregation of the housing reachable from it",
                        "(gravity-weighted over\nresident workers). Vertical:",
                        "segregation of the neighborhoods that actually supply",
                        "its workers. Euclidean-gravity benchmark."),
       x = "segregation of the ACCESSIBLE housing",
       y = "segregation of the REALIZED labor shed (wres)") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank())
ggsave(file.path(DIR_CO_FIG, "p4_fig_workbased_menu.png"), pM,
       width = 7.4, height = 6.2, dpi = 350, bg = "white")

## ---- PART C: the Cervero link at workplace scale -----------------------------
if ("z_mean_dist_km_work" %in% names(xw)) {
  for (y in c("z_mean_dist_km_work", Y)) {
    f1 <- feols(as.formula(sprintf("%s ~ z_pct_reslow_of_res + %s | county_fips",
                                   y, COVS)), data = xw, cluster = ~jurisd_main)
    f2 <- feols(as.formula(sprintf(
      "%s ~ z_pct_reslow_of_res + z_dist_cbd_km + z_dist_empctr_km + %s | county_fips",
      y, COVS)), data = xw, cluster = ~jurisd_main)
    message(sprintf("  %s ~ work-tract restrictiveness: naive %.3f (p=%.3f) | position %.3f (p=%.3f)",
      y, coef(f1)["z_pct_reslow_of_res"],
      summary(f1)$coeftable["z_pct_reslow_of_res", 4],
      coef(f2)["z_pct_reslow_of_res"],
      summary(f2)$coeftable["z_pct_reslow_of_res", 4]))
  }
}
message("87 complete.")
