# ==============================================================================
# 75_co_network_access.R      [PAPER 4, step 15]
# ONE carefully designed network-accessibility analysis (not a transportation-
# variable expansion): does how the transportation system connects
# neighborhoods to the labor market -- not merely being peripheral -- track
# the residential -> workplace-location segregation coupling?
#
#  PART A (runs now; no new data)
#    Opportunity-set vs realized-destination decomposition. For each origin
#    tract j, compute the workplace-location segregation of its ACCESSIBLE
#    labor market:
#        accD_j = sum_i jobs_i * exp(-0.10 d_ji) * D_wac_i
#                 / sum_i jobs_i * exp(-0.10 d_ji)
#    (impedance-weighted composition of what is reachable), alongside the
#    REALIZED wexp_j, and the sorting gap = wexp - accD. Two stories become
#    distinguishable: an OPPORTUNITY-STRUCTURE story (the reachable labor
#    market is itself segregated: accD high, gap ~ 0) vs a SORTING-WITHIN-
#    OPPORTUNITIES story (integrated work is reachable but actual jobs are
#    disproportionately elsewhere: gap > 0). The near-equal effective
#    destination counts across terciles already say destination QUANTITY is not the issue; this asks about composition.
#    CAVEAT (state it): the impedance here is Euclidean-gravity, a
#    geographic, not network, opportunity set.
#
#  PART B (network measures; requires clean/p3_tract_sld.rds from
#  paper_pipeline/53_sld.R). The accessibility ladder, extending 68's:
#      A  zoning only                     (68, A_total)
#      B  geographic position             (68, B_plus_regional_position)
#      C  NETWORK accessibility instead:  wexp ~ res_seg x log(D5AR) [auto
#         45-min network jobs] and x log(D5BR) [transit 45-min jobs]
#      D  position + network accessibility together: does connectivity
#         explain variation AMONG neighborhoods at comparable positions?
#    Plus the interaction the exemplar-pair figure invites: does job
#    accessibility itself moderate res_seg -> wexp (negative = accessibility
#    loosens the coupling), separately for auto and transit access.
#    The SLD table (paper_pipeline/53_sld.R) is required.
#
# LANGUAGE DISCIPLINE (as agreed for position): infrastructure, zoning,
# residential sorting and employment location co-evolved. All of this is
# mechanism-consistent decomposition of a cross-sectional association --
# never mediation, never causal pathway identification.
#
# Output: output/models/p4_network_access_ladder.csv
#         output/models/p4_opportunity_vs_realized.csv
#         output/models/p4_access_marginal_effects.csv
#         output/figures/p4_fig_opportunity_sorting.png
#         output/figures/p4_fig_transit_marginal.png   (Figure S10)
# ==============================================================================

source("60_co_setup.R")
library(fixest)

## ---- frames ------------------------------------------------------------------
acc  <- readRDS(file.path(DIR_CO_CLEAN, "co_accessibility_2023.rds"))
dat  <- readRDS(file.path(DIR_CO_OUT, "co_analysis_home_panel.rds"))
cent <- readRDS(file.path(DIR_CLEAN, "tract_centroids_km.rds")) |>
  transmute(tract_id = as.character(GEOID), CBSA_Code, X_km, Y_km)
seg23 <- readRDS(file.path(DIR_P2_OUT, "p2_tract_segregation_panel.rds")) |>
  filter(year == CO_ANCHOR_YEAR) |>
  transmute(tract_id, d_wac = d_whiteblack_wac_half)
wac <- readRDS(file.path(DIR_CO_CLEAN,
                         sprintf("co_wac_tract_%d.rds", CO_ANCHOR_YEAR))) |>
  transmute(tract_id = as.character(tract_id), jobs = C000)

## ---- PART A: accessible-opportunity composition ------------------------------
reg <- cent |> filter(CBSA_Code == CO_CBSA_MAIN) |>
  inner_join(wac,  by = "tract_id") |>
  inner_join(seg23, by = "tract_id") |>
  filter(!is.na(d_wac))
XY <- as.matrix(reg[, c("X_km", "Y_km")])
D  <- as.matrix(dist(XY))
W  <- exp(-0.10 * D) * matrix(reg$jobs, nrow(reg), nrow(reg), byrow = TRUE)
accD <- as.numeric((W %*% reg$d_wac) / rowSums(W))
opp <- tibble(tract_id = reg$tract_id, accD = accD)

zscore <- function(x) { x[!is.finite(x)] <- NA_real_; as.numeric(scale(x)) }
xs <- dat |>
  filter(year == CO_ANCHOR_YEAR, in_scope, denver_msa,
         n_commuters >= P3_MIN_COMMUTERS) |>
  left_join(acc, by = "tract_id") |>
  left_join(opp, by = "tract_id") |>
  mutate(sorting_gap = wexp_whiteblack_wac_half - accD,
         tercile = ntile(pct_reslow_of_res, 3)) |>
  mutate(across(c(wexp_whiteblack_wac_half, d_whiteblack_rac_half,
                  pct_reslow_of_res, accD, sorting_gap,
                  pct_black_rac, pct_lowincome_rac, log_worker_density_rac,
                  income_percapita_k, income_percapita_k_sq,
                  dist_cbd_km, dist_empctr_km, log_jobs_grav),
                zscore, .names = "z_{.col}"))

# decomposition descriptives: is high exposure an opportunity-structure fact
# (accD elevated) or a sorting fact (gap elevated)?
desc <- xs |>
  mutate(res_decile = ntile(d_whiteblack_rac_half, 10)) |>
  group_by(res_decile) |>
  summarise(n = n(),
            realized = mean(wexp_whiteblack_wac_half, na.rm = TRUE),
            accessible = mean(accD, na.rm = TRUE),
            gap = mean(sorting_gap, na.rm = TRUE), .groups = "drop")
by_terc <- xs |>
  group_by(tercile) |>
  summarise(n = n(),
            realized = mean(wexp_whiteblack_wac_half, na.rm = TRUE),
            accessible = mean(accD, na.rm = TRUE),
            gap = mean(sorting_gap, na.rm = TRUE), .groups = "drop")
write.csv(bind_rows(desc |> mutate(split = "res_seg_decile") |>
                      rename(group = res_decile),
                    by_terc |> mutate(split = "exclusionary_tercile") |>
                      rename(group = tercile)),
          file.path(DIR_CO_MOD, "p4_opportunity_vs_realized.csv"),
          row.names = FALSE)
print(as.data.frame(by_terc))

COVS <- paste(c("z_pct_black_rac", "z_pct_lowincome_rac",
                "z_log_worker_density_rac", "z_income_percapita_k",
                "z_income_percapita_k_sq"), collapse = " + ")
X <- "z_d_whiteblack_rac_half"
tidy1 <- function(fit, id, keep) {
  ct <- as.data.frame(summary(fit)$coeftable)
  ct$term <- rownames(ct); rownames(ct) <- NULL
  names(ct)[1:4] <- c("estimate", "std.error", "statistic", "p.value")
  ct |> filter(term %in% keep) |> mutate(model_id = id, n_obs = fit$nobs)
}
safe <- function(fml, id, keep, d = xs)
  tryCatch(tidy1(feols(as.formula(fml), data = d, cluster = ~jurisd_main),
                 id, keep), error = function(e) NULL)
res <- list()
# does within-opportunity sorting rise with residential segregation?
res[[1]] <- safe(sprintf("z_sorting_gap ~ %s + %s | county_fips", X, COVS),
                 "A_gap_on_resseg", X)
# and with exclusionary zoning?
res[[2]] <- safe(sprintf(
  "z_sorting_gap ~ z_pct_reslow_of_res + %s | county_fips", COVS),
  "A_gap_on_zoning", "z_pct_reslow_of_res")
# variance split of realized exposure between opportunity and sorting
# wexp = accD + gap, so var(wexp) = cov(wexp, accD) + cov(wexp, gap): a
# covariance allocation, exact by construction; the two parts are correlated.
vd <- xs |> filter(!is.na(accD), !is.na(sorting_gap))
write.csv(tibble(
  n = nrow(vd),
  pct_var_accessible = 100 * cov(vd$wexp_whiteblack_wac_half, vd$accD) /
    var(vd$wexp_whiteblack_wac_half),
  pct_var_sorting_gap = 100 * cov(vd$wexp_whiteblack_wac_half, vd$sorting_gap) /
    var(vd$wexp_whiteblack_wac_half),
  cor_accessible_gap = cor(vd$accD, vd$sorting_gap)),
  file.path(DIR_CO_MOD, "p4_variance_allocation.csv"), row.names = FALSE)
message(sprintf(
  "var(wexp) split: cov with accD %.0f%% | cov with gap %.0f%% | r(accD,gap)=%.2f",
  100 * cov(vd$wexp_whiteblack_wac_half, vd$accD) /
    var(vd$wexp_whiteblack_wac_half),
  100 * cov(vd$wexp_whiteblack_wac_half, vd$sorting_gap) /
    var(vd$wexp_whiteblack_wac_half),
  cor(vd$accD, vd$sorting_gap)))

## ---- PART B: the accessibility ladder ----------------------------------------
sld_file <- file.path(DIR_CLEAN, "p3_tract_sld.rds")
if (file.exists(sld_file) || "sld_D5AR" %in% names(xs)) {
  # 63 joins the SLD when the file exists, so the panel may already carry the
  # columns; join only if it does not (a second join would create .x/.y pairs)
  if (!"sld_D5AR" %in% names(xs))
    xs <- xs |> left_join(readRDS(sld_file), by = "tract_id")
  xs <- xs |>
    mutate(log_d5ar = ifelse(sld_D5AR > 0, log(sld_D5AR), NA_real_),
           log_d5br = ifelse(sld_D5BR > 0, log(sld_D5BR), NA_real_),
           z_log_d5ar = zscore(log_d5ar), z_log_d5br = zscore(log_d5br))
  NET <- c(auto_45min = "z_log_d5ar", transit_45min = "z_log_d5br")
  message("SLD network measures found: auto (D5AR) + transit (D5BR).")
} else {
  stop("p3_tract_sld.rds not found -- run paper_pipeline/53_sld.R first.")
}
POS <- sprintf(paste("z_dist_cbd_km + z_dist_empctr_km +",
                     "%s:z_dist_cbd_km + %s:z_dist_empctr_km"), X, X)
for (nm in names(NET)) {
  v <- NET[[nm]]
  keep <- c(sprintf("%s:%s", X, "z_pct_reslow_of_res"),
            sprintf("%s:%s", X, v), v)
  # C: network accessibility replaces geographic position
  res[[length(res) + 1]] <- safe(sprintf(
    "z_wexp_whiteblack_wac_half ~ %s * z_pct_reslow_of_res + %s + %s + %s:%s | county_fips",
    X, COVS, v, X, v), sprintf("C_network_%s", nm), keep)
  # D: position AND network accessibility together
  res[[length(res) + 1]] <- safe(sprintf(
    "z_wexp_whiteblack_wac_half ~ %s * z_pct_reslow_of_res + %s + %s + %s + %s:%s | county_fips",
    X, COVS, POS, v, X, v), sprintf("D_position_plus_%s", nm), keep)
  # the exemplar-pair question, on its own: does accessibility moderate the
  # coupling? (negative interaction = accessibility loosens it)
  res[[length(res) + 1]] <- safe(sprintf(
    "z_wexp_whiteblack_wac_half ~ %s * %s + %s | county_fips", X, v, COVS),
    sprintf("ACCESSxSEG_%s", nm), c(X, sprintf("%s:%s", X, v), v))
}
all_res <- bind_rows(res)
write.csv(all_res, file.path(DIR_CO_MOD, "p4_network_access_ladder.csv"),
          row.names = FALSE)
print(as.data.frame(all_res |> filter(grepl(":", term))))

## ---- figure: opportunity vs realized -----------------------------------------
pF <- ggplot(xs |> filter(!is.na(accD)),
             aes(accD, wexp_whiteblack_wac_half,
                 color = factor(tercile, 1:3, c("least", "middle", "most")))) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = .4,
              color = "grey50") +
  geom_point(size = 1.1, alpha = .6) +
  scale_color_manual(values = c("#c6dbef", "#6baed6", "#08519c"),
                     name = "exclusionary tercile") +
  labs(title = "Accessible opportunity vs realized destinations",
       subtitle = paste(
         "Above the 45-degree line: workers sort into MORE segregated",
         "locations than their reachable labor market;",
         "\nbelow: less. Euclidean-gravity opportunity set (see 75 header)."),
       x = "workplace-location segregation of the ACCESSIBLE labor market",
       y = "realized workplace-location exposure (wexp)") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")
ggsave(file.path(DIR_CO_FIG, "p4_fig_opportunity_sorting.png"), pF,
       width = 7.6, height = 6.4, dpi = 350, bg = "white")

## ---- Figure S10: marginal effect of residential segregation across access -----
# The ACCESSxSEG interaction as a picture: the slope of res_seg -> wexp
# evaluated across the 5th-95th percentiles of network access, with
# delta-method 95% CIs from the cluster-robust (jurisdiction) vcov. The figure
# shows the specification AS FITTED; 75b_co_access_influence.R reports how
# sensitive these interactions are to single tracts and to the functional form
# of access. A moderation pattern in a cross-sectional decomposition, not a
# causal effect of transport service.
if (all(c("z_log_d5ar", "z_log_d5br") %in% names(xs))) {
  me_grid <- function(fit, xterm, vterm, vseq) {
    b <- coef(fit); V <- vcov(fit)   # vcov inherits ~jurisd_main clustering
    it <- intersect(c(paste0(xterm, ":", vterm), paste0(vterm, ":", xterm)),
                    names(b))
    stopifnot(length(it) == 1)
    slope <- b[[xterm]] + b[[it]] * vseq
    se <- sqrt(V[xterm, xterm] + vseq^2 * V[it, it] + 2 * vseq * V[xterm, it])
    tibble(access_z = vseq, slope = slope, se = se,
           lo = slope - 1.96 * se, hi = slope + 1.96 * se)
  }
  MODE_LAB <- c(auto_45min    = "auto (45-min network jobs)",
                transit_45min = "transit (45-min network jobs)")
  me <- list(); rugs <- list()
  for (nm in names(NET)) {
    v <- NET[[nm]]
    fit <- tryCatch(feols(as.formula(sprintf(
      "z_wexp_whiteblack_wac_half ~ %s * %s + %s | county_fips", X, v, COVS)),
      data = xs, cluster = ~jurisd_main), error = function(e) NULL)
    if (is.null(fit)) next
    vv <- xs[[v]][is.finite(xs[[v]])]
    vseq <- seq(quantile(vv, .05), quantile(vv, .95), length.out = 41)
    me[[nm]] <- me_grid(fit, X, v, vseq) |>
      mutate(mode = MODE_LAB[[nm]], n_obs = fit$nobs)
    # rug only over the plotted grid: a handful of extreme low-access
    # outliers otherwise stretch the shared axis and empty the panels
    rugs[[nm]] <- tibble(access_z = vv[vv >= min(vseq) & vv <= max(vseq)],
                         mode = MODE_LAB[[nm]])
  }
  me <- bind_rows(me)
  write.csv(me, file.path(DIR_CO_MOD, "p4_access_marginal_effects.csv"),
            row.names = FALSE)
  pM <- ggplot(me, aes(access_z, slope)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = .4,
               color = "grey50") +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#6baed6", alpha = .25) +
    geom_line(color = "#08519c", linewidth = .8) +
    geom_rug(data = bind_rows(rugs), aes(access_z), inherit.aes = FALSE,
             sides = "b", alpha = .15, length = unit(0.02, "npc")) +
    facet_wrap(~mode, scales = "free_x") +
    labs(title = "How network access reshapes the coupling",
         subtitle = paste(
           "Slope of residential segregation on workplace-location exposure",
           "(both in SD), across the access distribution;",
           "\n95% delta-method CIs, jurisdiction-clustered. Transit sample",
           "restricted to tracts with recorded transit access."),
         x = "network job accessibility (SD of log 45-min jobs)",
         y = "marginal effect of residential segregation (SD)") +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"))
  ggsave(file.path(DIR_CO_FIG, "p4_fig_transit_marginal.png"), pM,
         width = 8.6, height = 4.4, dpi = 350, bg = "white")
  message("Figure S10 (marginal effects by access) written.")
} else {
  message("SLD columns absent -- skipping Figure S10 (marginal effects).")
}
message("75 complete.")
