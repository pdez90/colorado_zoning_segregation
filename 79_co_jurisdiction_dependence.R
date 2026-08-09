# ==============================================================================
# 79_co_jurisdiction_dependence.R   [PAPER 4, step 19 -- SI extensions 4 + 5 + 6]
# THE KEY QUESTION (Priyanka): do municipalities that restrict housing
# nevertheless depend on workers housed elsewhere?
#
#  Ext. 5 -- jurisdictional boundary accounting. Per jurisdiction: flow-based
#    worker imports/exports, net exchange, the earnings composition of
#    IMPORTED vs RESIDENT-RETAINED workers, and the mean commute distance of
#    each. Does local land-use autonomy separate the benefits of employment
#    from the responsibility of housing the workers who provide it?
#  Ext. 6 -- housing exclusion x low-wage labor dependence. Among jurisdictions
#    hosting substantial low-wage employment, does residential restrictiveness
#    predict the share of low-wage jobs filled by in-commuters? Job-weighted
#    correlation + regression; headline aggregates for exclusionary (>=80%
#    res-low) jurisdictions. Figure: p4_fig_lowwage_dependence.png.
#  Ext. 4 -- balance vs matching. Conventional jobs-housing balance
#    (jobs / resident workers) against actual flow retention (% of residents
#    working in the jurisdiction): two tracts/jurisdictions at "balance" can
#    differ wholly in matching. Typology (quadrants at balance = 1 and
#    retention = regional median) + figure: p4_fig_balance_matching.png.
#
# Conventions as 77: jurisdiction = tract's jurisd_main (dominant zoning
# jurisdiction; approximation); flows between zoning-frame tracts; SE01 =
# jobs <= $1,250/month. Descriptive; no causal identification claimed.
# Output: output/models/p4_jurisdiction_dependence.csv
#         output/models/p4_balance_matching_typology.csv
#         output/figures/p4_fig_lowwage_dependence.png
#         output/figures/p4_fig_balance_matching.png
# ==============================================================================

source("60_co_setup.R")
library(ggrepel)

## ---- frames (as 77) ----------------------------------------------------------
dat <- readRDS(file.path(DIR_CO_OUT, "co_analysis_home_panel.rds"))
xs <- dat |>
  filter(year == CO_ANCHOR_YEAR, in_scope, denver_msa,
         n_commuters >= P3_MIN_COMMUTERS) |>
  select(tract_id, jurisd_main, pct_reslow_of_res)
jmap <- xs |> select(tract_id, jurisd_main)
cent <- readRDS(file.path(DIR_CLEAN, "tract_centroids_km.rds")) |>
  transmute(tract_id = as.character(GEOID), X_km, Y_km)
rac <- readRDS(file.path(DIR_CO_CLEAN,
                         sprintf("co_rac_tract_%d.rds", CO_ANCHOR_YEAR))) |>
  transmute(tract_id, workers_res = C000, res_lowwage = CE01)
wac <- readRDS(file.path(DIR_CO_CLEAN,
                         sprintf("co_wac_tract_%d.rds", CO_ANCHOR_YEAR))) |>
  transmute(tract_id, jobs = C000, jobs_lowwage = CE01)

od <- map(P3_OD_PARTS, function(part) {
  for (f in c(p3_od_cache(part, CO_ANCHOR_YEAR, "co"),
              co_od_cache(part, CO_ANCHOR_YEAR)))
    if (file.exists(f)) return(readRDS(f))
  NULL
}) |> compact() |> bind_rows() |>
  group_by(w_tract, h_tract) |>
  summarise(S000 = sum(S000, na.rm = TRUE),
            SE01 = sum(SE01, na.rm = TRUE), .groups = "drop") |>
  inner_join(jmap |> rename(h_jur = jurisd_main), by = c(h_tract = "tract_id")) |>
  inner_join(jmap |> rename(w_jur = jurisd_main), by = c(w_tract = "tract_id")) |>
  left_join(cent |> rename(xh = X_km, yh = Y_km), by = c(h_tract = "tract_id")) |>
  left_join(cent |> rename(xw = X_km, yw = Y_km), by = c(w_tract = "tract_id")) |>
  mutate(same_jur = h_jur == w_jur,
         dist_km = sqrt((xh - xw)^2 + (yh - yw)^2))

## ---- ext. 5: per-jurisdiction boundary accounting ----------------------------
work_side <- od |>
  group_by(w_jur) |>
  summarise(
    inflow_total    = sum(S000),
    inflow_imported = sum(S000[!same_jur]),
    low_total       = sum(SE01),
    low_imported    = sum(SE01[!same_jur]),
    lowwage_share_imported = 100 * sum(SE01[!same_jur]) /
                                   pmax(sum(S000[!same_jur]), 1),
    lowwage_share_retained = 100 * sum(SE01[same_jur]) /
                                   pmax(sum(S000[same_jur]), 1),
    dist_imported = weighted.mean(dist_km[!same_jur], S000[!same_jur],
                                  na.rm = TRUE),
    dist_retained = weighted.mean(dist_km[same_jur], S000[same_jur],
                                  na.rm = TRUE), .groups = "drop") |>
  rename(jurisd_main = w_jur)
home_side <- od |>
  group_by(h_jur) |>
  summarise(outflow_total = sum(S000),
            retained      = sum(S000[same_jur]), .groups = "drop") |>
  rename(jurisd_main = h_jur)
stocks <- xs |>
  left_join(rac, by = "tract_id") |>
  left_join(wac, by = "tract_id") |>
  group_by(jurisd_main) |>
  summarise(n_tracts = n(),
            workers_housed = sum(workers_res, na.rm = TRUE),
            res_lowwage = sum(res_lowwage, na.rm = TRUE),
            jobs_hosted = sum(jobs, na.rm = TRUE),
            lowwage_jobs = sum(jobs_lowwage, na.rm = TRUE),
            mean_pct_reslow = mean(pct_reslow_of_res, na.rm = TRUE),
            .groups = "drop")

jur <- stocks |>
  left_join(work_side, by = "jurisd_main") |>
  left_join(home_side, by = "jurisd_main") |>
  mutate(
    jobs_per_worker = jobs_hosted / pmax(workers_housed, 1),
    pct_residents_retained = 100 * retained / pmax(outflow_total, 1),
    pct_jobs_imported = 100 * inflow_imported / pmax(inflow_total, 1),
    pct_lowwage_jobs_imported = 100 * low_imported / pmax(low_total, 1),
    net_flow_import = inflow_total - outflow_total,
    lowwage_share_resworkers = 100 * res_lowwage / pmax(workers_housed, 1))
write.csv(jur, file.path(DIR_CO_MOD, "p4_jurisdiction_dependence.csv"),
          row.names = FALSE)

## ---- ext. 6: the key test ----------------------------------------------------
key <- jur |> filter(lowwage_jobs >= 500)
wcor <- with(key, {
  w <- lowwage_jobs
  cx <- weighted.mean(mean_pct_reslow, w); cy <- weighted.mean(pct_lowwage_jobs_imported, w)
  sum(w * (mean_pct_reslow - cx) * (pct_lowwage_jobs_imported - cy)) /
    sqrt(sum(w * (mean_pct_reslow - cx)^2) * sum(w * (pct_lowwage_jobs_imported - cy)^2))
})
fit <- lm(pct_lowwage_jobs_imported ~ mean_pct_reslow, data = key,
          weights = lowwage_jobs)
message(sprintf(
  "KEY TEST (n=%d jurisdictions with >=500 low-wage jobs): weighted r = %.2f; slope = %.2f pp imported per pp res-low (p = %.4f)",
  nrow(key), wcor, coef(fit)[2], summary(fit)$coefficients[2, 4]))
excl <- jur |> filter(mean_pct_reslow >= 80)
message(sprintf(
  "Exclusionary jurisdictions (>=80%% res-low, n=%d): host %s low-wage jobs; %.1f%% filled by in-commuters (region overall: %.1f%%). They hold %.1f%% of region jobs but house %.1f%% of region workers.",
  nrow(excl), format(sum(excl$lowwage_jobs), big.mark = ","),
  100 * sum(excl$low_imported) / sum(excl$low_total),
  100 * sum(jur$low_imported, na.rm = TRUE) / sum(jur$low_total, na.rm = TRUE),
  100 * sum(excl$jobs_hosted) / sum(jur$jobs_hosted),
  100 * sum(excl$workers_housed) / sum(jur$workers_housed)))
message(sprintf(
  "Imported vs retained workforce, exclusionary jurisdictions: low-wage share %.1f%% vs %.1f%%; mean commute %.1f vs %.1f km.",
  with(excl, 100 * sum(low_imported) / pmax(sum(inflow_imported), 1)),
  with(excl, 100 * (sum(low_total) - sum(low_imported)) /
         pmax(sum(inflow_total) - sum(inflow_imported), 1)),
  with(excl, weighted.mean(dist_imported, inflow_imported, na.rm = TRUE)),
  with(excl, weighted.mean(dist_retained, inflow_total - inflow_imported,
                           na.rm = TRUE))))

pK <- ggplot(key, aes(mean_pct_reslow, pct_lowwage_jobs_imported)) +
  geom_smooth(aes(weight = lowwage_jobs), method = "lm", formula = y ~ x,
              color = "#08519c", fill = "#c6dbef", linewidth = .7) +
  geom_point(aes(size = lowwage_jobs), alpha = .65, color = "#08519c") +
  geom_text_repel(aes(label = jurisd_main), size = 2.4, max.overlaps = 12,
                  seed = 4) +
  scale_size_area(max_size = 9, name = "low-wage jobs hosted",
                  labels = scales::comma) +
  labs(title = "Housing exclusion and low-wage labor dependence",
       subtitle = paste("Jurisdictions hosting >= 500 low-wage jobs;",
                        "fit is low-wage-job-weighted least squares."),
       x = "% of residential land zoned low-density (jurisdiction mean)",
       y = "% of local low-wage jobs filled by in-commuters") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")
ggsave(file.path(DIR_CO_FIG, "p4_fig_lowwage_dependence.png"), pK,
       width = 7.8, height = 6.2, dpi = 350, bg = "white")

## ---- ext. 4: balance vs matching ---------------------------------------------
med_ret <- median(jur$pct_residents_retained, na.rm = TRUE)
jur <- jur |>
  mutate(typology = case_when(
    jobs_per_worker >= 1 & pct_residents_retained >= med_ret ~
      "employment center, self-supplying",
    jobs_per_worker >= 1 & pct_residents_retained <  med_ret ~
      "employment center, worker-importing",
    jobs_per_worker <  1 & pct_residents_retained >= med_ret ~
      "residential, locally employed",
    TRUE ~ "bedroom exporter"))
typ <- jur |>
  group_by(typology) |>
  summarise(n = n(),
            mean_pct_reslow = weighted.mean(mean_pct_reslow, workers_housed),
            mean_lowwage_import = weighted.mean(pct_lowwage_jobs_imported,
                                                lowwage_jobs, na.rm = TRUE),
            workers_housed = sum(workers_housed),
            jobs_hosted = sum(jobs_hosted),
            .groups = "drop")
write.csv(typ, file.path(DIR_CO_MOD, "p4_balance_matching_typology.csv"),
          row.names = FALSE)
print(as.data.frame(typ), digits = 3)

pT <- ggplot(jur |> filter(workers_housed >= 2000),
             aes(jobs_per_worker, pct_residents_retained)) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey55", linewidth = .4) +
  geom_hline(yintercept = med_ret, linetype = 2, color = "grey55",
             linewidth = .4) +
  geom_point(aes(size = workers_housed, color = mean_pct_reslow), alpha = .8) +
  geom_text_repel(aes(label = jurisd_main), size = 2.4, max.overlaps = 14,
                  seed = 4) +
  scale_x_log10() +
  scale_color_gradient(low = "#c6dbef", high = "#08519c",
                       name = "% res. land low-density") +
  scale_size_area(max_size = 10, name = "resident workers",
                  labels = scales::comma) +
  labs(title = "Jobs-housing balance is not matching",
       subtitle = paste("Quadrant lines: balance = 1 job/resident worker;",
                        "retention = regional median. Jurisdictions >= 2,000",
                        "resident workers."),
       x = "local jobs per resident worker (log scale)",
       y = "% of employed residents working within the jurisdiction") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        legend.box = "vertical")
ggsave(file.path(DIR_CO_FIG, "p4_fig_balance_matching.png"), pT,
       width = 8.2, height = 6.8, dpi = 350, bg = "white")
message("79 complete.")
