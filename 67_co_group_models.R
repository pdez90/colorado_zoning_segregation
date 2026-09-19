# ==============================================================================
# 67_co_group_models.R      [COLORADO ZONING case study, step 7]
# Who travels farther, who lands in mixed vs monoculture workhoods, and does
# ZONING structure the class gaps? (fixest, 55/64 conventions)
#
#  Table G1  2023 group profile: distance, >24km share, destination
#            portfolio, wexp, destination earnings/education/industry mix --
#            by earnings class, industry supergroup, age band; overall and
#            by exclusionary-zoning tercile. THE "who travels farther" table.
#  BLOCK GA  Within-tract class gaps ~ zoning (2023 cross-section, county FE,
#            jurisdiction-clustered): gap_dist = mean_dist(SE01)-(SE03);
#            gap_wexp; gap_portfolio. Does exclusionary zoning widen the gap
#            between the low-wage workers who SERVE a place and the high-wage
#            workers who LIVE its lifestyle? ("drive till you qualify" /
#            service-worker exile.)
#  BLOCK GB  Destination class-mixing: dest_earn/edu/ind_entropy (all
#            workers) ~ res_seg x zoning -- the mixing analogue of ZB: do
#            residents of exclusionary-zoned segregated tracts commute into
#            class-homogeneous workhoods?
#  BLOCK GC  Group-specific ZB: wexp_g ~ res_seg x pct_reslow_of_res fitted
#            separately per group -- is the zoning amplification of the
#            res->workhood link borne by low-earnings workers?
#  Sensitivities: min group commuters (suppression noise floors), no-density
#            covariates, industry-supergroup gaps, age gaps.
#
# Education note (stated wherever used): education has NO OD flows; it
# enters only (a) as destination-mix entropy (a WAC margin attribute --
# assumption-free) and (b) via home-tract CD weights in descriptive rows
# (shared-destination assumption, as in Paper 3's race handling).
#
# Needs: 66 outputs + 63's home panel (zoning + covariates).
# Output: output/models/co_group_table_g1.csv
#         output/models/co_group_model_coefficients.csv + figures
# ==============================================================================

source("60_co_setup.R")
library(fixest)

gp   <- readRDS(file.path(DIR_CO_CLEAN, "co_group_flows_panel.rds"))
dat  <- readRDS(file.path(DIR_CO_OUT, "co_analysis_home_panel.rds"))

CO_MIN_GROUP_COMMUTERS <- 10   # suppression-noise floor for group cells

zscore <- function(x) { x[!is.finite(x)] <- NA_real_; as.numeric(scale(x)) }
tidy_fx <- function(fit, model_id, keep = "z_|:") {
  if (is.null(fit)) return(NULL)
  ct <- as.data.frame(summary(fit)$coeftable)
  ct$term <- rownames(ct); rownames(ct) <- NULL
  names(ct)[1:4] <- c("estimate", "std.error", "statistic", "p.value")
  ct |> filter(grepl(keep, term)) |>
    mutate(conf.low  = estimate - 1.96 * std.error,
           conf.high = estimate + 1.96 * std.error,
           model_id = model_id, n_obs = fit$nobs)
}
safe_feols <- function(fml, data, cluster, model_id, ...) {
  fit <- tryCatch(feols(as.formula(fml), data = data, cluster = cluster, ...),
                  error = function(e) { message("  FAILED ", model_id, ": ",
                                                conditionMessage(e)); NULL })
  tidy_fx(fit, model_id)
}
res <- list()
add <- function(x) if (!is.null(x)) res[[length(res) + 1]] <<- x

## ---- zoning + covariate frame (2023 anchor, as 64) ---------------------------
zx <- dat |>
  filter(year == CO_ANCHOR_YEAR, in_scope, denver_msa,
         n_commuters >= P3_MIN_COMMUTERS) |>   # the 656-tract analysis frame
  select(tract_id, CBSAFP, county_fips, jurisd_main,
         pct_res_low, pct_reslow_of_res, pct_adu_res, zoning_entropy,
         pct_job_zone, d_whiteblack_rac_half,
         pct_black_rac, pct_lowincome_rac, log_worker_density_rac,
         income_percapita_k, income_percapita_k_sq) |>
  mutate(reslow_tercile = ntile(pct_reslow_of_res, 3))

## =============================================================================
## Table G1 -- who travels farther / into what, 2023
## =============================================================================
g23 <- gp |>
  filter(year == CO_ANCHOR_YEAR, n_commuters >= CO_MIN_GROUP_COMMUTERS) |>
  inner_join(zx |> select(tract_id, reslow_tercile), by = "tract_id")

g1_overall <- g23 |>
  group_by(group) |>
  summarise(
    tercile = "all", n_tracts = n(), workers = sum(n_commuters),
    across(c(mean_dist_km, pct_gt24km, eff_n_dest, wexp_whiteblack,
             dest_earn_entropy, dest_edu_entropy, dest_ind_entropy),
           ~ weighted.mean(.x, n_commuters, na.rm = TRUE)),
    .groups = "drop")
g1_tercile <- g23 |>
  group_by(group, tercile = as.character(reslow_tercile)) |>
  summarise(
    n_tracts = n(), workers = sum(n_commuters),
    across(c(mean_dist_km, pct_gt24km, eff_n_dest, wexp_whiteblack,
             dest_earn_entropy, dest_edu_entropy, dest_ind_entropy),
           ~ weighted.mean(.x, n_commuters, na.rm = TRUE)),
    .groups = "drop")
# education rows (SHARED-DESTINATION ASSUMPTION -- education has no OD
# flows; each home tract's all-worker measures are weighted by its RAC
# education counts, i.e. we assume destinations don't differ by education
# within a home tract; same handling as Paper 3's race rows)
rac_wt <- readRDS(file.path(DIR_CO_CLEAN, "co_rac_weights_panel.rds"))
g_all23 <- g23 |> filter(group == "all") |>
  inner_join(rac_wt |> filter(year == CO_ANCHOR_YEAR) |>
               select(tract_id, CD01, CD04), by = "tract_id")
g1_edu <- map(c(edu_less_hs = "CD01", edu_ba_plus = "CD04"),
              function(col) {
  w <- g_all23[[col]]
  g_all23 |>
    summarise(
      tercile = "all", n_tracts = n(), workers = sum(w),
      across(c(mean_dist_km, pct_gt24km, eff_n_dest, wexp_whiteblack,
               dest_earn_entropy, dest_edu_entropy, dest_ind_entropy),
             ~ weighted.mean(.x, w, na.rm = TRUE)))
}) |> bind_rows(.id = "group") |>
  mutate(group = paste0(group, "_shared_dest"))

g1 <- bind_rows(g1_overall, g1_tercile, g1_edu) |>
  mutate(across(where(is.numeric), ~ round(.x, 3)))
write.csv(g1, file.path(DIR_CO_MOD, "co_group_table_g1.csv"),
          row.names = FALSE)
message("Table G1 written (", nrow(g1), " rows)")

## =============================================================================
## Build the wide gap frame (2023)
## =============================================================================
gw <- g23 |>
  select(tract_id, group, n_commuters, mean_dist_km, pct_gt24km,
         eff_n_dest, wexp_whiteblack, dest_earn_entropy, dest_edu_entropy,
         dest_ind_entropy) |>
  pivot_wider(names_from = group,
              values_from = c(n_commuters, mean_dist_km, pct_gt24km,
                              eff_n_dest, wexp_whiteblack, dest_earn_entropy,
                              dest_edu_entropy, dest_ind_entropy)) |>
  inner_join(zx, by = "tract_id") |>
  mutate(
    gap_dist_earn   = mean_dist_km_earn_low - mean_dist_km_earn_high,
    gap_wexp_earn   = wexp_whiteblack_earn_low - wexp_whiteblack_earn_high,
    gap_dest_earn   = eff_n_dest_earn_low - eff_n_dest_earn_high,
    gap_dist_age    = mean_dist_km_age_u30 - mean_dist_km_age_55p,
    gap_dist_goods_serv = mean_dist_km_ind_goods - mean_dist_km_ind_services)

MODEL_VARS <- c("gap_dist_earn", "gap_wexp_earn", "gap_dest_earn",
                "gap_dist_age", "gap_dist_goods_serv",
                "dest_earn_entropy_all", "dest_edu_entropy_all",
                "dest_ind_entropy_all",
                "mean_dist_km_earn_low", "mean_dist_km_earn_high",
                "wexp_whiteblack_earn_low", "wexp_whiteblack_earn_high",
                "pct_res_low", "pct_reslow_of_res", "pct_adu_res",
                "zoning_entropy", "pct_job_zone", "d_whiteblack_rac_half",
                "pct_black_rac", "pct_lowincome_rac",
                "log_worker_density_rac", "income_percapita_k",
                "income_percapita_k_sq")
gz <- gw |> mutate(across(any_of(MODEL_VARS), zscore, .names = "z_{.col}"))

COVS <- paste(c("z_pct_black_rac", "z_pct_lowincome_rac",
                "z_log_worker_density_rac", "z_income_percapita_k",
                "z_income_percapita_k_sq"), collapse = " + ")
COVS_ND <- paste(c("z_pct_black_rac", "z_pct_lowincome_rac",
                   "z_income_percapita_k", "z_income_percapita_k_sq"),
                 collapse = " + ")
ZMODS <- c("pct_reslow_of_res", "pct_adu_res", "zoning_entropy",
           "pct_job_zone")

message("Gap frame: ", nrow(gz), " tracts")

## =============================================================================
## BLOCK GA -- class gaps ~ zoning
## =============================================================================
for (out in c("z_gap_dist_earn", "z_gap_wexp_earn", "z_gap_dest_earn",
              "z_gap_dist_age", "z_gap_dist_goods_serv"))
  for (zv in ZMODS) {
    add(safe_feols(sprintf("%s ~ z_%s + %s | county_fips", out, zv, COVS),
                   gz, ~jurisd_main,
                   sprintf("GA_%s_%s", sub("^z_gap_", "", out), zv)))
    add(safe_feols(sprintf("%s ~ z_%s + %s | county_fips", out, zv, COVS_ND),
                   gz, ~jurisd_main,
                   sprintf("GA_nodens_%s_%s", sub("^z_gap_", "", out), zv)))
  }
# levels alongside gaps (is the gap from the low end rising or the high end
# falling?)
for (out in c("z_mean_dist_km_earn_low", "z_mean_dist_km_earn_high"))
  add(safe_feols(sprintf(
    "%s ~ z_pct_reslow_of_res + %s | county_fips", out, COVS),
    gz, ~jurisd_main, sprintf("GA_level_%s_reslow", sub("^z_", "", out))))

## =============================================================================
## BLOCK GB -- destination class-mixing ~ res_seg x zoning (mixing analogue
## of ZB)
## =============================================================================
for (out in c("z_dest_earn_entropy_all", "z_dest_edu_entropy_all",
              "z_dest_ind_entropy_all"))
  for (zv in ZMODS) {
    add(safe_feols(sprintf("%s ~ z_%s + %s | county_fips", out, zv, COVS),
                   gz, ~jurisd_main,
                   sprintf("GB_%s_%s", gsub("^z_dest_|_all$", "", out), zv)))
    add(safe_feols(sprintf(
      "%s ~ z_d_whiteblack_rac_half * z_%s + %s | county_fips",
      out, zv, COVS), gz, ~jurisd_main,
      sprintf("GB_interact_%s_%s",
              sub("^z_dest_", "", sub("_all$", "", out)), zv)))
  }

## =============================================================================
## BLOCK GC -- group-specific ZB: who bears the zoning amplification?
## =============================================================================
for (g in c("earn_low", "earn_mid", "earn_high", "ind_goods",
            "ind_tradetrans", "ind_services")) {
  ycol <- sprintf("wexp_whiteblack_%s", g)
  if (!ycol %in% names(gw)) next
  d_g <- gw |>
    filter(.data[[sprintf("n_commuters_%s", g)]] >= CO_MIN_GROUP_COMMUTERS) |>
    mutate(z_y = zscore(.data[[ycol]]),
           across(any_of(c("d_whiteblack_rac_half", "pct_reslow_of_res",
                           "pct_black_rac", "pct_lowincome_rac",
                           "log_worker_density_rac", "income_percapita_k",
                           "income_percapita_k_sq")),
                  zscore, .names = "z_{.col}"))
  add(safe_feols(paste(
    "z_y ~ z_d_whiteblack_rac_half * z_pct_reslow_of_res +", COVS,
    "| county_fips"), d_g, ~jurisd_main, sprintf("GC_zb_%s", g)))
}

## ---- write tables ------------------------------------------------------------
all_res <- bind_rows(res)
write.csv(all_res,
          file.path(DIR_CO_MOD, "co_group_model_coefficients.csv"),
          row.names = FALSE)
message("Fitted ", n_distinct(all_res$model_id), " group models.")

## ---- figures -----------------------------------------------------------------
theme_co <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())
GROUP_LABELS <- c(
  all = "All workers", earn_low = "Earnings: low (<=$1,250/mo)",
  earn_mid = "Earnings: middle", earn_high = "Earnings: high (>$3,333/mo)",
  ind_goods = "Industry: goods-producing",
  ind_tradetrans = "Industry: trade/transport/utilities",
  ind_services = "Industry: services",
  age_u30 = "Age: under 30", age_30_54 = "Age: 30-54", age_55p = "Age: 55+")

# Fig G1: distance by group x exclusionary tercile (dumbbell-ish dot plot)
fg1 <- g1_tercile |>
  mutate(glab = GROUP_LABELS[group],
         tercile = factor(tercile, 1:3, c("least", "middle", "most")))
pG1 <- ggplot(fg1, aes(mean_dist_km, reorder(glab, mean_dist_km),
                       color = tercile)) +
  geom_line(aes(group = glab), color = "grey80", linewidth = 1.5) +
  geom_point(size = 2.4) +
  scale_color_manual(values = c("#1b6ca8", "#888888", "#c8442c"),
                     name = "Exclusionary-zoning\ntercile") +
  labs(title = "Who commutes farther, and where is it worst? (2023)",
       subtitle = "Commuter-weighted mean distance by LODES flow group and home-tract zoning",
       x = "mean commute distance (km)", y = NULL) + theme_co
ggsave(file.path(DIR_CO_FIG, "co_figG1_distance_by_group.png"), pG1,
       width = 9, height = 5.5, dpi = 300, bg = "white")

# Fig G2: earnings-gap coefficient forest (GA, with + without density)
fb <- all_res |>
  filter(grepl("^GA_(nodens_)?dist_earn_", model_id),
         grepl("^z_pct_|^z_zoning", term))
if (nrow(fb) > 0) {
  pG2 <- fb |>
    mutate(spec = ifelse(grepl("nodens", model_id), "no density",
                         "with density"),
           mod = sub("^GA_(nodens_)?dist_earn_", "", model_id)) |>
    ggplot(aes(estimate, reorder(mod, estimate), color = spec)) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = .3) +
    geom_pointrange(aes(xmin = conf.low, xmax = conf.high),
                    position = position_dodge(width = .4)) +
    scale_color_manual(values = c("#c8442c", "#1b6ca8"), name = NULL) +
    labs(title = "Does zoning widen the low- vs high-earnings commute gap?",
         subtitle = "Outcome: mean_dist(SE01) - mean_dist(SE03), SD units; jurisdiction-clustered CI",
         x = "coefficient on zoning measure", y = NULL) + theme_co
  ggsave(file.path(DIR_CO_FIG, "co_figG2_earnings_gap_forest.png"), pG2,
         width = 8, height = 4.5, dpi = 300, bg = "white")
}

# Fig G3: destination class-mixing vs res seg by zoning tercile
fg3 <- gw |>
  mutate(reslow_tercile = factor(reslow_tercile, 1:3,
                                 c("least", "middle", "most")),
         bin = ntile(d_whiteblack_rac_half, 20)) |>
  group_by(reslow_tercile, bin) |>
  summarise(res = mean(d_whiteblack_rac_half, na.rm = TRUE),
            mix = mean(dest_earn_entropy_all, na.rm = TRUE),
            .groups = "drop")
pG3 <- ggplot(fg3, aes(res, mix, color = reslow_tercile)) +
  geom_point(size = 1.2) +
  geom_smooth(method = "lm", se = FALSE, linewidth = .6) +
  scale_color_manual(values = c("#1b6ca8", "#888888", "#c8442c"),
                     name = "Exclusionary-zoning\ntercile") +
  labs(title = "Residential segregation and the class mix of destination workhoods (2023)",
       subtitle = "Flow-weighted earnings entropy of destination tracts; binned means",
       x = "residential D (RAC, beta=0.5)",
       y = "destination earnings-mix entropy") + theme_co
ggsave(file.path(DIR_CO_FIG, "co_figG3_class_mixing.png"), pG3,
       width = 8, height = 5.5, dpi = 300, bg = "white")

## ---- DIAGNOSTICS -------------------------------------------------------------
write_codiag(
  all_res |>
    mutate(block = sub("_.*", "", model_id),
           sig = p.value < 0.05,
           direction = ifelse(estimate > 0, "pos", "neg")) |>
    count(block, direction, sig),
  "67_group_model_sign_summary")
write_codiag(
  g23 |> count(group, name = "n_tract_cells"),
  "67_group_cell_counts")
# Eyeball: earn_low cells << earn_high cells is EXPECTED (suppression +
# smaller group); if any group has < ~300 tracts the gap models for it are
# underpowered -- read CIs accordingly.
message("67 complete. Tables in ", DIR_CO_MOD)
