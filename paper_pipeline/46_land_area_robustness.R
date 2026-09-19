# ==============================================================================
# 46_land_area_robustness.R
# Tract-size robustness of the preferred (position-adjusted) specification.
#
# Local environments in the segregation index are tract-bounded and tract area
# grows toward the metropolitan periphery, so the exclusionary interaction could
# reflect geometry. This script refits 68's B_plus_regional_position rung
#   wexp ~ X * z_pct_reslow_of_res + COVS + POS | county_fips
# (P0) adding log land area (P1) and log land area x residential segregation
# (P2), and reports the correlation of tract area with CBD distance.
#
# Output: diagnostics/46_land_area_robustness.csv
#   cd ~/Downloads/LODES/paper_pipeline
#   Rscript -e "source('46_land_area_robustness.R')"
# ==============================================================================

source("50_p3_setup.R")
suppressPackageStartupMessages(library(fixest))

CO_DIR <- path.expand("~/Downloads/LODES/Colorado")
dat <- readRDS(file.path(CO_DIR, "output/co_analysis_home_panel.rds"))

## ---- 1. find the land-area column --------------------------------------------
AREA_CANDS <- c("aland_km2", "ALAND_km2", "aland", "ALAND", "area_km2", "land_km2")
acol <- AREA_CANDS[AREA_CANDS %in% names(dat)][1]

if (is.na(acol)) {
  message("land area not in the panel -- searching clean/ for the 61 output ...")
  cands <- list.files(file.path(CO_DIR, "clean"), pattern = "[.]rds$",
                      full.names = TRUE)
  got <- NULL
  for (f in cands) {
    d <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(d) || !is.data.frame(d)) next
    hit <- AREA_CANDS[AREA_CANDS %in% names(d)][1]
    if (!is.na(hit) && "tract_id" %in% names(d)) {
      message("  found ", hit, " in ", basename(f))
      got <- d |> select(tract_id, dplyr::all_of(hit)) |> distinct(tract_id, .keep_all = TRUE)
      acol <- hit; break
    }
  }
  if (is.null(got)) stop("No land-area column found in the panel or in clean/. ",
                         "Re-run 61_co_zoning_tract.R first.")
  dat <- dat |> left_join(got, by = "tract_id")
}
message("land-area column: ", acol)

## ---- 2. the frame, exactly as 68 builds it ------------------------------------
zscore <- function(x) { x[!is.finite(x)] <- NA_real_; as.numeric(scale(x)) }
acc <- readRDS(file.path(CO_DIR, "clean/co_accessibility_2023.rds"))
xs <- dat |>
  filter(year == 2023L, in_scope, denver_msa, n_commuters >= 20L) |>
  left_join(acc, by = "tract_id") |>
  mutate(land_km2 = suppressWarnings(as.numeric(.data[[acol]])),
         log_area = ifelse(land_km2 > 0, log(land_km2), NA_real_)) |>
  mutate(across(c(wexp_whiteblack_wac_half, d_whiteblack_rac_half,
                  pct_reslow_of_res, pct_black_rac, pct_lowincome_rac,
                  log_worker_density_rac, income_percapita_k,
                  income_percapita_k_sq, dist_cbd_km, dist_empctr_km,
                  log_area),
                zscore, .names = "z_{.col}"))
message("frame: ", nrow(xs), " tracts | land area present for ",
        sum(is.finite(xs$land_km2)))

## ---- 3. the correlation the manuscript quotes ---------------------------------
ok <- is.finite(xs$land_km2) & is.finite(xs$dist_cbd_km)
r_raw <- cor(xs$land_km2[ok], xs$dist_cbd_km[ok])
r_log <- cor(xs$log_area[ok], xs$dist_cbd_km[ok], use = "complete.obs")
message("\ncorrelation of tract land area with CBD distance:")
message(sprintf("  raw km2   r = %.3f", r_raw))
message(sprintf("  log km2   r = %.3f", r_log))

## ---- 4. the three specifications ----------------------------------------------
COVS <- paste(c("z_pct_black_rac", "z_pct_lowincome_rac",
                "z_log_worker_density_rac", "z_income_percapita_k",
                "z_income_percapita_k_sq"), collapse = " + ")
X   <- "z_d_whiteblack_rac_half"
POS <- sprintf(paste("z_dist_cbd_km + z_dist_empctr_km +",
                     "%s:z_dist_cbd_km + %s:z_dist_empctr_km"), X, X)
KEY <- paste0(X, ":z_pct_reslow_of_res")

fit <- function(extra, id) tryCatch({
  f <- feols(as.formula(sprintf(
    "z_wexp_whiteblack_wac_half ~ %s * z_pct_reslow_of_res + %s + %s%s | county_fips",
    X, COVS, POS, extra)), data = xs, cluster = ~jurisd_main)
  ct <- as.data.frame(summary(f)$coeftable); ct$term <- rownames(ct)
  names(ct)[1:4] <- c("estimate", "std.error", "statistic", "p.value")
  ct |> filter(term %in% c(KEY, "z_log_area", paste0(X, ":z_log_area"))) |>
    mutate(model_id = id, n_obs = f$nobs)
}, error = function(e) { message("!! ", id, ": ", conditionMessage(e)); NULL })

res <- bind_rows(
  fit("",                                   "P0_preferred"),
  fit(" + z_log_area",                      "P1_plus_log_area"),
  fit(sprintf(" + z_log_area + %s:z_log_area", X), "P2_area_interacted"))

message("\n", strrep("=", 74))
message("TRACT-SIZE ROBUSTNESS -- current values")
message(strrep("=", 74), "\n")
print(as.data.frame(res |> transmute(model_id, term = substr(term, 1, 42),
  estimate = round(estimate, 4), se = round(std.error, 4),
  p = round(p.value, 4), n = n_obs)), row.names = FALSE)

key <- res |> filter(term == KEY)
if (nrow(key) >= 3) {
  se0 <- key$std.error[key$model_id == "P0_preferred"]
  se2 <- key$std.error[key$model_id == "P2_area_interacted"]
  message(sprintf(
    "\nstandard error on the key interaction: %.4f preferred -> %.4f with the",
    se0, se2))
  message(sprintf("area interaction, a %.0f%% increase.", 100 * (se2 / se0 - 1)))
}

write_diag(res, "46_land_area_robustness")

message("46_land_area_robustness.R complete.")
