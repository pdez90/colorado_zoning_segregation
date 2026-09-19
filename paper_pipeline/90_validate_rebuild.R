# ==============================================================================
# 90_validate_rebuild.R
# Does the rebuilt pipeline reproduce the committed results?
#
# Compares what 61-64 have just written on this machine against the reference
# copies committed to github.com/pdez90/colorado_zoning_segregation. Point it
# at a clone:
#   Sys.setenv(REF_DIR = "~/Downloads/colorado_zoning_segregation")
#   source("90_validate_rebuild.R")
#
# Reports, per comparison: how many rows matched, the largest absolute and
# relative difference, and the five worst rows. Writes the whole thing to
# diagnostics/90_rebuild_validation.csv.
#
# HOW TO READ IT. The OD job totals and the tract counts should match to the
# last digit -- they depend only on LODES and on geography, both of which are
# fully reconstructed. The segregation-dependent quantities (mean_wexp, the
# res-vs-wexp correlations, the model coefficients) should match to about
# three decimals; anything larger means the index is not being computed the
# same way and should be chased before anything else. Coefficients whose
# specification includes income_percapita may drift slightly more, because
# pre-2021 ACS income is crosswalked from 2010 tracts here (see 35). A drift
# that only shows up in income-adjusted models, and only in the third decimal,
# is the crosswalk; a drift everywhere is not.
# ==============================================================================

# 60 lives in the case-study folder and sources relatively, so borrow its wd
CO_DIR <- path.expand(Sys.getenv("CO_DIR", "~/Downloads/LODES/Colorado"))
if (!file.exists(file.path(CO_DIR, "60_co_setup.R")))
  stop("60_co_setup.R not found in ", CO_DIR,
       " -- set Sys.setenv(CO_DIR = \"<path>\").")
.owd <- setwd(CO_DIR)
source("60_co_setup.R")      # 60 itself hops to paper_pipeline and back
setwd(.owd)

REF <- path.expand(Sys.getenv("REF_DIR",
                              "~/Downloads/LODES/colorado_zoning_segregation"))
if (!dir.exists(REF))
  stop("Reference clone not found at ", REF,
       "\n  git clone https://github.com/pdez90/colorado_zoning_segregation\n",
       "  then set Sys.setenv(REF_DIR = \"<path>\").")

results <- list()

rd <- function(p) if (file.exists(p)) suppressWarnings(
  readr::read_csv(p, show_col_types = FALSE)) else NULL

compare <- function(label, rel_path, keys, values, tol) {
  ref <- rd(file.path(REF, rel_path))
  new <- rd(file.path(DIR_CO, rel_path))
  if (is.null(ref) || is.null(new)) {
    message(sprintf("  [skip] %-34s %s", label,
                    if (is.null(ref)) "no reference file" else
                      "not produced yet on this machine"))
    return(invisible(NULL))
  }
  keys   <- intersect(keys,   intersect(names(ref), names(new)))
  values <- intersect(values, intersect(names(ref), names(new)))
  if (!length(keys) || !length(values)) {
    message(sprintf("  [skip] %-34s columns do not line up", label))
    return(invisible(NULL))
  }
  j <- inner_join(ref |> select(all_of(c(keys, values))),
                  new |> select(all_of(c(keys, values))),
                  by = keys, suffix = c(".ref", ".new"))
  per_val <- map(values, function(v) {
    a <- j[[paste0(v, ".ref")]]; b <- j[[paste0(v, ".new")]]
    ok <- is.finite(a) & is.finite(b)
    tibble(comparison = label, variable = v,
           n_rows = sum(ok),
           n_ref_only = nrow(ref) - nrow(j),
           max_abs_diff = if (any(ok)) max(abs(a[ok] - b[ok])) else NA_real_,
           max_rel_diff = if (any(ok))
             max(abs(a[ok] - b[ok]) / pmax(abs(a[ok]), 1e-9)) else NA_real_,
           tol = tol,
           verdict = dplyr::case_when(
             !any(ok) ~ "NO OVERLAP",
             max(abs(a[ok] - b[ok])) <= tol ~ "match",
             TRUE ~ "DIFFERS"))
  }) |> bind_rows()
  results[[label]] <<- per_val
  for (i in seq_len(nrow(per_val)))
    message(sprintf("  %-7s %-34s %-24s max|d| = %.6g  (n = %d)",
                    per_val$verdict[i], label, per_val$variable[i],
                    per_val$max_abs_diff[i], per_val$n_rows[i]))
  bad <- per_val |> filter(verdict == "DIFFERS")
  if (nrow(bad)) for (v in bad$variable) {
    a <- j[[paste0(v, ".ref")]]; b <- j[[paste0(v, ".new")]]
    worst <- j[order(-abs(a - b))[1:min(5, nrow(j))], keys, drop = FALSE]
    worst$ref <- a[order(-abs(a - b))[1:min(5, nrow(j))]]
    worst$new <- b[order(-abs(a - b))[1:min(5, nrow(j))]]
    message("    worst rows for ", v, ":")
    print(as.data.frame(worst), row.names = FALSE)
  }
  invisible(per_val)
}

message("\n== Stage 1: LODES + geography (should be exact) ==")
compare("OD jobs by year", "diagnostics/62_cross_msa_jobs_dropped_by_year.csv",
        keys = "year",
        values = c("jobs_both_metro", "jobs_same_msa", "pct_dropped_cross_msa"),
        tol = 1e-6)

message("\n== Stage 2: the segregation index (the load-bearing rebuild) ==")
compare("wexp by CBSA-year", "diagnostics/62_home_panel_by_cbsa_year.csv",
        keys = c("cbsa", "year"),
        values = c("n_tracts", "med_commuters", "mean_wexp"),
        tol = 5e-4)
compare("res-vs-wexp correlation", "diagnostics/62_res_vs_wexp_correlation_by_year.csv",
        keys = "year", values = "r_resseg_wexp", tol = 5e-4)

message("\n== Stage 3: the analysis panel ==")
compare("panel completeness", "diagnostics/63_home_panel_completeness_by_year.csv",
        keys = "year",
        values = c("n_tracts", "pctna_wexp_whiteblack_wac_half",
                   "pctna_d_whiteblack_rac_half", "pctna_pct_res_low",
                   "pctna_pct_reslow_of_res", "pctna_income_percapita"),
        tol = 0.05)
compare("zoning coverage by county", "diagnostics/61_coverage_by_county.csv",
        keys = intersect(names(rd(file.path(REF, "diagnostics/61_coverage_by_county.csv"))),
                         c("county_fips", "county", "COUNTYFP")),
        values = c("n_tracts", "n_covered", "mean_cover"),
        tol = 1e-6)

message("\n== Stage 4: published estimates ==")
compare("main model coefficients", "output/models/co_model_coefficients.csv",
        keys = c("model_id", "term"),
        values = c("estimate", "std.error", "p.value", "n_obs"),
        tol = 5e-3)
## ---- coefficients, split by what the term is ---------------------------------
# A single max|d| over all 700+ coefficients hides the thing you actually want
# to know: whether the ZONING and SEGREGATION terms reproduce. Income is a
# reconstruction (see 35) and enters as a quadratic, which is tail-sensitive,
# so it is expected to drift on its own coefficient. Read the three rows.
local({
  fa <- file.path(REF,    "output/models/co_model_coefficients.csv")
  fb <- file.path(DIR_CO, "output/models/co_model_coefficients.csv")
  if (!file.exists(fa) || !file.exists(fb)) return(invisible(NULL))
  a <- rd(fa) |> filter(!is.na(estimate))
  b <- rd(fb) |> filter(!is.na(estimate))
  j <- inner_join(a |> select(model_id, term, estimate, p.value),
                  b |> select(model_id, term, estimate, p.value),
                  by = c("model_id", "term"), suffix = c(".ref", ".new")) |>
    mutate(d = abs(estimate.ref - estimate.new),
           kind = case_when(
             grepl("income_percapita", term) ~ "income (reconstructed)",
             grepl(paste0("pct_res_low|pct_reslow|pct_adu|zoning_entropy|",
                          "pct_job_zone|pct_industrial|pct_commercial|",
                          "d_whiteblack|d_hisp|d_college|res_seg|wexp"),
                   term) ~ "zoning / segregation",
             TRUE ~ "other controls"))
  message("\n-- coefficient differences by term type --")
  print(as.data.frame(
    j |> group_by(kind) |>
      summarise(n = n(),
                median_abs_diff = round(median(d), 4),
                max_abs_diff    = round(max(d), 4),
                within_0.005    = sprintf("%d/%d", sum(d <= 0.005), n()),
                .groups = "drop")), row.names = FALSE)

  flips <- j |>
    filter(kind == "zoning / segregation",
           (p.value.ref < .05) != (p.value.new < .05))
  message("\n-- significance flips (p < .05) on zoning / segregation terms: ",
          nrow(flips), " --")
  if (nrow(flips))
    print(as.data.frame(flips |>
      transmute(model_id, term = substr(term, 1, 44),
                ref = round(estimate.ref, 4), p_ref = round(p.value.ref, 4),
                new = round(estimate.new, 4), p_new = round(p.value.new, 4))),
      row.names = FALSE)
})

compare("group model coefficients", "output/models/co_group_model_coefficients.csv",
        keys = c("model_id", "term"),
        values = c("estimate", "std.error", "p.value"), tol = 5e-3)
compare("2023 descriptives", "output/models/co_descriptives_2023.csv",
        keys = "reslow_tercile",
        values = c("n", "pct_res_low", "pct_adu_res", "zoning_entropy",
                   "pct_job_zone", "d_whiteblack_rac_half",
                   "wexp_whiteblack_wac_half", "mean_dist_km", "eff_n_dest",
                   "jobs_housing_ratio", "pct_black_rac",
                   "income_percapita_k"),
        tol = 5e-3)

## ---- summary -----------------------------------------------------------------
all_res <- bind_rows(results)
if (nrow(all_res)) {
  write_codiag(all_res, "90_rebuild_validation")
  message("\n== SUMMARY ==")
  print(as.data.frame(all_res |> count(verdict)), row.names = FALSE)
  if (any(all_res$verdict == "DIFFERS")) {
    message("\nStill differing:")
    print(as.data.frame(all_res |> filter(verdict == "DIFFERS") |>
                          select(comparison, variable, max_abs_diff)),
          row.names = FALSE)
  } else message("\nEverything that could be compared matches.")
} else message("Nothing to compare yet -- run 61-64 first.")

message("90_validate_rebuild.R complete.")
