# ==============================================================================
# 75b_co_access_influence.R
# Are the interaction estimates carried by the sample or by a few tracts?
#
# Both network-accessibility measures (SLD D5AR auto, D5BR transit) enter 75 as
# z-scored logs. Job-accessibility counts have a floor near zero, so on the log
# scale a handful of fringe tracts sit many standard deviations below the mean;
# residential segregation is itself strongly right-skewed. An interaction
# between two such variables can be set by one or two observations. This script
# reports, for every network-accessibility model in 75 and for the headline
# zoning model in 68:
#
#   PART 1  the estimate under each SLD construction written by 53
#           (bridge x no-data rule) and each functional form of access:
#             log         z(log x), x > 0                   [75's specification]
#             log_wins    as log, winsorized at 2.5 / 97.5 %
#             log_trim    as log, tracts below 1,000 reachable jobs excluded
#             rank        z(rank of x), x > 0
#             log1p       z(log(1 + x)), all tracts (zero-access tracts kept)
#   PART 2  leave-one-out: the range of the interaction over all single-tract
#           deletions, and the ten most influential tracts
#   PART 3  the same leave-one-out for the headline zoning moderation (68's
#           A_total and B_plus_regional_position rungs)
#
# Reporting rule used in the manuscript: an interaction is described as a
# finding only if its sign and significance hold in PART 1 across forms under
# the primary construction AND no single deletion in PART 2 removes it.
#
# Output: output/models/p4_access_influence_forms.csv
#         output/models/p4_access_influence_loo.csv
#         output/models/p4_access_influence_top_tracts.csv
#         output/models/p4_headline_influence_loo.csv
# ==============================================================================

source("60_co_setup.R")
library(fixest)

PRIMARY <- "block_pop|zero"
acc <- readRDS(file.path(DIR_CO_CLEAN, "co_accessibility_2023.rds"))
dat <- readRDS(file.path(DIR_CO_OUT, "co_analysis_home_panel.rds"))
var <- readRDS(file.path(DIR_CLEAN, "p3_tract_sld_variants.rds"))

zscore <- function(x) { x[!is.finite(x)] <- NA_real_; as.numeric(scale(x)) }
# 75 fits on the SLD columns that 63 wrote into the panel; this script fits on
# the variants file. Stop if the two have drifted apart (53 rebuilt without 63).
if ("sld_D5AR" %in% names(dat)) {
  chk <- dat |> filter(year == CO_ANCHOR_YEAR) |>
    select(tract_id, a = sld_D5AR, b = sld_D5BR) |>
    inner_join(var |> filter(construction == PRIMARY), by = "tract_id")
  if (!isTRUE(all.equal(chk$a, chk$sld_D5AR)) || !isTRUE(all.equal(chk$b, chk$sld_D5BR)))
    stop("Panel SLD columns differ from the primary construction in ",
         "p3_tract_sld_variants.rds -- re-run 63 after 53.")
}
xs0 <- dat |>
  filter(year == CO_ANCHOR_YEAR, in_scope, denver_msa,
         n_commuters >= P3_MIN_COMMUTERS) |>
  select(-any_of(c("sld_D5AR", "sld_D5BR", "sld_NatWalkInd", "sld_pop",
                   "log_sld_D5BR"))) |>
  left_join(acc, by = "tract_id") |>
  mutate(across(c(wexp_whiteblack_wac_half, d_whiteblack_rac_half,
                  pct_reslow_of_res, pct_black_rac, pct_lowincome_rac,
                  log_worker_density_rac, income_percapita_k,
                  income_percapita_k_sq, dist_cbd_km, dist_empctr_km),
                zscore, .names = "z_{.col}"))
message("frame: ", nrow(xs0), " tracts")

X    <- "z_d_whiteblack_rac_half"
Y    <- "z_wexp_whiteblack_wac_half"
ZON  <- "z_pct_reslow_of_res"
COVS <- paste(c("z_pct_black_rac", "z_pct_lowincome_rac",
                "z_log_worker_density_rac", "z_income_percapita_k",
                "z_income_percapita_k_sq"), collapse = " + ")
POS  <- sprintf(paste("z_dist_cbd_km + z_dist_empctr_km +",
                      "%s:z_dist_cbd_km + %s:z_dist_empctr_km"), X, X)
# the three 75 specifications, with `z_acc` standing for the access measure
SPECS <- c(
  C_network      = sprintf("%s ~ %s * %s + %s + z_acc + %s:z_acc | county_fips",
                           Y, X, ZON, COVS, X),
  D_position_net = sprintf("%s ~ %s * %s + %s + %s + z_acc + %s:z_acc | county_fips",
                           Y, X, ZON, COVS, POS, X),
  ACCESSxSEG     = sprintf("%s ~ %s * z_acc + %s | county_fips", Y, X, COVS))
IT <- paste0(X, ":z_acc")

forms <- function(x) {
  lg <- ifelse(is.finite(x) & x > 0, log(x), NA_real_)
  q  <- quantile(lg, c(.025, .975), na.rm = TRUE)
  list(log      = lg,
       log_wins = pmin(pmax(lg, q[1]), q[2]),
       log_trim = ifelse(x >= 1000, lg, NA_real_),
       rank     = ifelse(is.na(lg), NA_real_, rank(lg, na.last = "keep")),
       log1p    = ifelse(is.finite(x), log1p(x), NA_real_))
}
# The clustered covariance is computed at estimation and read straight from the
# fitted object; summary() would re-evaluate the data argument lazily.
fit1 <- function(fml, d) {
  d <- d[!is.na(d$jurisd_main), ]
  d$jurisd_main <- as.character(d$jurisd_main)
  feols(as.formula(fml), data = d, vcov = ~jurisd_main, notes = FALSE)
}
grab <- function(fit, term) {
  ct <- fit$coeftable
  tibble(estimate = ct[term, 1], std.error = ct[term, 2], p.value = ct[term, 4],
         n_obs = fit$nobs)
}

## ---- PART 1: constructions x forms -------------------------------------------
out1 <- list()
for (cn in unique(var$construction)) {
  xs <- xs0 |> left_join(var |> filter(construction == cn) |>
                           select(tract_id, sld_D5AR, sld_D5BR), by = "tract_id")
  for (mode in c("auto", "transit")) {
    f <- forms(if (mode == "auto") xs$sld_D5AR else xs$sld_D5BR)
    for (fn in names(f)) {
      d <- xs |> mutate(z_acc = zscore(f[[fn]]))
      for (sp in names(SPECS)) {
        r <- tryCatch(grab(fit1(SPECS[[sp]], d), IT), error = function(e) NULL)
        if (!is.null(r)) out1[[length(out1) + 1]] <- r |>
          mutate(construction = cn, primary = cn == PRIMARY, mode = mode,
                 form = fn, spec = sp, .before = 1)
      }
    }
  }
}
out1 <- bind_rows(out1)
write.csv(out1, file.path(DIR_CO_MOD, "p4_access_influence_forms.csv"),
          row.names = FALSE)
message("\nPRIMARY construction, interaction res.seg x access:")
print(as.data.frame(out1 |> filter(primary) |>
  transmute(mode, spec, form, estimate = round(estimate, 4),
            p = round(p.value, 3), n_obs)))

## ---- PART 2: leave-one-out, primary construction, 75's log form --------------
# Refit with each tract deleted in turn, keeping the jurisdiction-clustered
# p-value of every refit, so the reporting rule in the header can be applied
# to ALL deletions and not only to the one that moves the estimate furthest.
loo <- function(fml, d, term) {
  vars <- all.vars(as.formula(fml))
  d <- d[stats::complete.cases(d[, c(vars, "jurisd_main")]), ]
  full <- grab(fit1(fml, d), term)
  ep <- vapply(seq_len(nrow(d)), function(i) {
    g <- grab(fit1(fml, d[-i, ]), term)
    c(g$estimate, g$p.value)
  }, numeric(2))
  list(full = full, d = d, est = ep[1, ], p = ep[2, ])
}
xsP <- xs0 |> left_join(var |> filter(construction == PRIMARY) |>
                          select(tract_id, sld_D5AR, sld_D5BR), by = "tract_id")
out2 <- list(); top <- list()
for (mode in c("auto", "transit")) {
  x <- if (mode == "auto") xsP$sld_D5AR else xsP$sld_D5BR
  d0 <- xsP |> mutate(acc_raw = x, z_acc = zscore(forms(x)$log))
  for (sp in names(SPECS)) {
    L <- loo(SPECS[[sp]], d0, IT)
    # the single deletion that moves the estimate furthest (either direction)
    worst <- which.max(abs(L$est - L$full$estimate))
    refit <- grab(fit1(SPECS[[sp]], L$d[-worst, ]), IT)
    out2[[length(out2) + 1]] <- tibble(
      mode = mode, spec = sp, n_obs = L$full$n_obs,
      estimate = L$full$estimate, p.value = L$full$p.value,
      loo_min = min(L$est), loo_max = max(L$est),
      sign_changes_under_one_deletion = any(sign(L$est) != sign(L$full$estimate)),
      loo_max_p = max(L$p), n_deletions_p_ge_05 = sum(L$p >= .05),
      most_influential_tract = L$d$tract_id[worst],
      estimate_without_it = refit$estimate, p_without_it = refit$p.value)
    o <- order(-abs(L$est - L$full$estimate))[1:10]
    top[[length(top) + 1]] <- L$d[o, ] |>
      transmute(mode = mode, spec = sp, tract_id, jurisd_main,
                access_raw = acc_raw, z_access = z_acc,
                z_res_seg = .data[[X]], z_wexp = .data[[Y]], dist_cbd_km,
                estimate_if_deleted = L$est[o],
                shift = L$est[o] - L$full$estimate)
  }
}
out2 <- bind_rows(out2); top <- bind_rows(top)
write.csv(out2, file.path(DIR_CO_MOD, "p4_access_influence_loo.csv"), row.names = FALSE)
write.csv(top,  file.path(DIR_CO_MOD, "p4_access_influence_top_tracts.csv"), row.names = FALSE)
message("\nLeave-one-out, primary construction:")
print(as.data.frame(out2 |> mutate(across(where(is.numeric), \(v) round(v, 4)))))

## ---- PART 3: the headline zoning moderation ----------------------------------
HEAD <- c(A_total = sprintf("%s ~ %s * %s + %s | county_fips", Y, X, ZON, COVS),
          B_plus_regional_position =
            sprintf("%s ~ %s * %s + %s + %s | county_fips", Y, X, ZON, COVS, POS))
ITZ <- paste0(X, ":", ZON)
out3 <- list()
for (sp in names(HEAD)) {
  L <- loo(HEAD[[sp]], xs0, ITZ)
  worst <- which.max(abs(L$est - L$full$estimate))
  refit <- grab(fit1(HEAD[[sp]], L$d[-worst, ]), ITZ)
  # residential segregation is right-skewed: also refit with it winsorized / ranked
  alt <- function(v) { d <- xs0; d[[X]] <- zscore(v); grab(fit1(HEAD[[sp]], d), ITZ) }
  rs <- xs0$d_whiteblack_rac_half
  w99 <- alt(pmin(rs, quantile(rs, .99, na.rm = TRUE)))
  rk  <- alt(rank(rs, na.last = "keep"))
  out3[[sp]] <- tibble(
    spec = sp, n_obs = L$full$n_obs, estimate = L$full$estimate,
    p.value = L$full$p.value, loo_min = min(L$est), loo_max = max(L$est),
    loo_max_p = max(L$p), n_deletions_p_ge_05 = sum(L$p >= .05),
    most_influential_tract = L$d$tract_id[worst],
    estimate_without_it = refit$estimate, p_without_it = refit$p.value,
    est_resseg_wins99 = w99$estimate, p_resseg_wins99 = w99$p.value,
    est_resseg_rank = rk$estimate, p_resseg_rank = rk$p.value)
}
out3 <- bind_rows(out3)
write.csv(out3, file.path(DIR_CO_MOD, "p4_headline_influence_loo.csv"), row.names = FALSE)
message("\nHeadline zoning moderation, influence:")
print(as.data.frame(out3 |> mutate(across(where(is.numeric), \(v) round(v, 4)))))

message("75b_co_access_influence.R complete.")
