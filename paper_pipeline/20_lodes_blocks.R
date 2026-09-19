# ==============================================================================
# 20_lodes_blocks.R
# Block-level LODES RAC and WAC for Colorado, 2011-2023, aggregated to the
# margins the segregation index needs. Also writes the tract-level OD caches
# the Colorado scripts expect at raw/od_tract_{part}_{yr}_co.rds.
#
# Produces:
#   clean/blocks_rac_{yr}.rds   block_id, tract_id, <margin columns>
#   clean/blocks_wac_{yr}.rds   same, workplace side
#   raw/od_tract_{part}_{yr}_co.rds   h_tract, w_tract, S000, SE01-03, SI01-03
#
# Downloads are cached; re-running skips anything already on disk.
# ==============================================================================

source("50_p3_setup.R")
suppressPackageStartupMessages(library(lehdr))

STATE <- "co"

keep_cols <- function(df, want) {
  have <- intersect(want, names(df))
  df |> select(all_of(c(intersect(c("block_id", "tract_id"), names(df)), have)))
}

## ---- RAC / WAC at block level ------------------------------------------------
for (yr in P3_YEARS) {
  for (ds in c("rac", "wac")) {
    f <- file.path(DIR_CLEAN, sprintf("blocks_%s_%s.rds", ds, yr))
    if (file.exists(f)) { message("  ", basename(f), " exists"); next }
    message("Downloading ", toupper(ds), " blocks ", yr, " ...")
    d <- lehdr::grab_lodes(state = STATE, year = yr, version = "LODES8",
                           lodes_type = ds, job_type = "JT01",
                           segment = "S000", agg_geo = "block")
    idcol <- if (ds == "rac") "h_geocode" else "w_geocode"
    if (!idcol %in% names(d)) idcol <- grep("geocode", names(d), value = TRUE)[1]
    d <- d |>
      mutate(block_id = as.character(.data[[idcol]]),
             tract_id = substr(block_id, 1, 11)) |>
      keep_cols(if (ds == "rac") P2_RAC_COLS else P2_WAC_COLS)
    saveRDS(d, f)
    message("  wrote ", basename(f), " (", nrow(d), " blocks)")
  }
}

## ---- tract-level OD ----------------------------------------------------------
for (yr in P3_YEARS) {
  for (part in P3_OD_PARTS) {
    f <- p3_od_cache(part, yr, STATE)
    if (file.exists(f)) { message("  ", basename(f), " exists"); next }
    message("Downloading OD ", part, " ", yr, " ...")
    d <- tryCatch(
      lehdr::grab_lodes(state = STATE, year = yr, version = "LODES8",
                        lodes_type = "od", job_type = "JT01",
                        segment = "S000", state_part = part,
                        agg_geo = "tract"),
      error = function(e) { message("  skipped (", conditionMessage(e), ")"); NULL })
    if (is.null(d)) next
    nm <- names(d)
    hc <- grep("^h_tract|^h_geocode", nm, value = TRUE)[1]
    wc <- grep("^w_tract|^w_geocode", nm, value = TRUE)[1]
    d <- d |>
      mutate(h_tract = substr(as.character(.data[[hc]]), 1, 11),
             w_tract = substr(as.character(.data[[wc]]), 1, 11)) |>
      group_by(h_tract, w_tract) |>
      summarise(across(any_of(P3_OD_COLS), ~ sum(.x, na.rm = TRUE)),
                .groups = "drop")
    saveRDS(d, f)
    message("  wrote ", basename(f), " (", nrow(d), " tract pairs)")
  }
}

message("20_lodes_blocks.R complete.")
