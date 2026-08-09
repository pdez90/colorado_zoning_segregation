# ==============================================================================
# 77_co_selfcontainment.R      [PAPER 4, step 17]
# Jobs-workers dependence and commuting self-containment: does "job-poor
# peripheral" mean structurally dependent on the rest of the metropolitan
# labor market? RAC + WAC + OD + zoning together, per the review:
#   - local jobs per resident worker            (WAC C000 / RAC C000)
#   - % resident workers working in their own jurisdiction      (OD)
#   - % local jobs filled by same-jurisdiction residents        (OD)
#   - % local LOW-WAGE jobs (SE01) filled by workers living
#     outside the jurisdiction                                  (OD group flows)
# summarized by exclusionary-zoning tercile for the paper (Table 1 extension),
# and per jurisdiction as groundwork for the separate "who houses whose
# workers?" jurisdictional paper.
#
# CONVENTIONS AND CAVEATS (state in any writeup):
#   - "Jurisdiction" is the tract's DOMINANT zoning jurisdiction (jurisd_main
#     from 61) -- an approximation for tracts split across jurisdictions.
#   - Flow measures use flows whose BOTH ends carry a jurisdiction (both
#     tracts in the zoning frame); the retained share is logged. Same-MSA
#     restriction is inherited from the OD frame itself (all CO flows read,
#     then filtered to frame tracts, which are Denver-MSA in-scope).
#   - Home-side shares are commuter-weighted; work-side shares job-weighted.
#   - SE01 = jobs earning <= $1,250/month (fixed nominal threshold).
#   - Descriptive only. No causal claims: jurisdictions' zoning, employment
#     bases, and commuting fields co-evolved.
#
# Output: output/models/p4_selfcontainment_tercile.csv       (paper numbers)
#         output/models/p4_selfcontainment_jurisdiction.csv  (future paper)
#         diagnostics/co_selfcontainment_coverage.csv
# ==============================================================================

source("60_co_setup.R")

## ---- frames ------------------------------------------------------------------
dat <- readRDS(file.path(DIR_CO_OUT, "co_analysis_home_panel.rds"))
xs <- dat |>
  filter(year == CO_ANCHOR_YEAR, in_scope, denver_msa,
         n_commuters >= P3_MIN_COMMUTERS) |>
  mutate(tercile = ntile(pct_reslow_of_res, 3)) |>
  select(tract_id, jurisd_main, tercile, pct_reslow_of_res)
stopifnot(!any(is.na(xs$jurisd_main)))
jmap <- xs |> select(tract_id, jurisd_main)

rac <- readRDS(file.path(DIR_CO_CLEAN,
                         sprintf("co_rac_tract_%d.rds", CO_ANCHOR_YEAR))) |>
  transmute(tract_id, workers_res = C000)
wac <- readRDS(file.path(DIR_CO_CLEAN,
                         sprintf("co_wac_tract_%d.rds", CO_ANCHOR_YEAR))) |>
  transmute(tract_id, jobs = C000, jobs_lowwage = CE01)

## ---- OD at the anchor year (reuse 51/62 caches; never re-download) -----------
od <- map(P3_OD_PARTS, function(part) {
  f_nat <- p3_od_cache(part, CO_ANCHOR_YEAR, "co")
  f_co  <- co_od_cache(part, CO_ANCHOR_YEAR)
  f <- if (file.exists(f_nat)) f_nat else f_co
  stopifnot(file.exists(f))
  readRDS(f)
}) |> bind_rows() |>
  group_by(w_tract, h_tract) |>
  summarise(S000 = sum(S000, na.rm = TRUE),
            SE01 = sum(SE01, na.rm = TRUE), .groups = "drop")

od2 <- od |>
  inner_join(jmap |> rename(h_jur = jurisd_main), by = c(h_tract = "tract_id")) |>
  inner_join(jmap |> rename(w_jur = jurisd_main), by = c(w_tract = "tract_id")) |>
  mutate(same_jur = h_jur == w_jur)
cov_row <- tibble(
  flows_total   = sum(od$S000),
  flows_in_frame = sum(od2$S000),
  pct_retained  = 100 * sum(od2$S000) / sum(od$S000))
write_codiag(cov_row, "co_selfcontainment_coverage")  # write_codiag appends .csv
message(sprintf("OD flows with both ends in the zoning frame: %.1f%%",
                cov_row$pct_retained))

## ---- tract-level pieces ------------------------------------------------------
home_side <- od2 |>
  group_by(h_tract) |>
  summarise(commuters = sum(S000),
            within_jur = sum(S000[same_jur]), .groups = "drop") |>
  rename(tract_id = h_tract)
work_side <- od2 |>
  group_by(w_tract) |>
  summarise(jobs_flows = sum(S000),
            jobs_local = sum(S000[same_jur]),
            low_flows  = sum(SE01),
            low_import = sum(SE01[!same_jur]), .groups = "drop") |>
  rename(tract_id = w_tract)

tr <- xs |>
  left_join(rac, by = "tract_id") |>
  left_join(wac, by = "tract_id") |>
  left_join(home_side, by = "tract_id") |>
  left_join(work_side, by = "tract_id") |>
  mutate(across(c(jobs, jobs_lowwage, jobs_flows, jobs_local,
                  low_flows, low_import), ~ replace_na(.x, 0)))

## ---- the paper numbers: by exclusionary tercile ------------------------------
by_terc <- tr |>
  group_by(tercile) |>
  summarise(
    n_tracts = n(),
    jobs_per_worker = mean(jobs / pmax(workers_res, 1), na.rm = TRUE),
    pct_work_own_jur = 100 * sum(within_jur, na.rm = TRUE) /
                             sum(commuters,  na.rm = TRUE),
    pct_jobs_filled_locally = 100 * sum(jobs_local) / pmax(sum(jobs_flows), 1),
    pct_lowwage_jobs_imported = 100 * sum(low_import) / pmax(sum(low_flows), 1),
    .groups = "drop")
write.csv(by_terc, file.path(DIR_CO_MOD, "p4_selfcontainment_tercile.csv"),
          row.names = FALSE)
print(as.data.frame(by_terc), digits = 3)

## ---- groundwork for the jurisdictional paper ---------------------------------
by_jur <- tr |>
  group_by(jurisd_main) |>
  summarise(
    n_tracts = n(),
    workers_housed = sum(workers_res, na.rm = TRUE),
    jobs_hosted = sum(jobs),
    lowwage_jobs = sum(jobs_lowwage),
    jobs_per_worker = jobs_hosted / pmax(workers_housed, 1),
    mean_pct_reslow = mean(pct_reslow_of_res, na.rm = TRUE),
    pct_residents_retained = 100 * sum(within_jur, na.rm = TRUE) /
                                   pmax(sum(commuters, na.rm = TRUE), 1),
    pct_jobs_filled_locally = 100 * sum(jobs_local) / pmax(sum(jobs_flows), 1),
    pct_lowwage_jobs_imported = 100 * sum(low_import) / pmax(sum(low_flows), 1),
    .groups = "drop") |>
  arrange(desc(jobs_hosted))
write.csv(by_jur, file.path(DIR_CO_MOD, "p4_selfcontainment_jurisdiction.csv"),
          row.names = FALSE)
cat("\nTop 12 jurisdictions by jobs hosted:\n")
print(as.data.frame(head(by_jur, 12)), digits = 3)
cat("\nMost exclusionary jurisdictions (mean % res-low >= 80, >= 3 tracts):\n")
print(as.data.frame(by_jur |>
        filter(mean_pct_reslow >= 80, n_tracts >= 3) |>
        arrange(desc(mean_pct_reslow)) |> head(12)), digits = 3)
message("77 complete.")
