# ==============================================================================
# 45_manuscript_number_audit.R
#
# Dump EVERY numeric value the pipeline currently produces into one tidy file,
# so each of the manuscript's numeric claims can be checked against it
# mechanically rather than by eye.
#
# Reads every .csv in:
#   Colorado/output/models/      the committed result tables
#   Colorado/diagnostics/        the per-script diagnostics
#   LODES/diagnostics/           the upstream diagnostics
#
# Writes:
#   Colorado/output/models/ALL_CURRENT_VALUES.csv
#     file | row_key | column | value
#   where row_key is built from that file's non-numeric (label) columns, so a
#   value can be traced back to the row it came from.
#
#   cd ~/Downloads/LODES/paper_pipeline
#   Rscript -e "source('45_manuscript_number_audit.R')"
# ==============================================================================

source("50_p3_setup.R")

CO_DIR <- path.expand(Sys.getenv("CO_DIR", "~/Downloads/LODES/Colorado"))
DIRS <- c(file.path(CO_DIR, "output/models"),
          file.path(CO_DIR, "diagnostics"),
          file.path(path.expand(Sys.getenv("LODES_ROOT", "~/Downloads/LODES")), "diagnostics"))

files <- unlist(lapply(DIRS[dir.exists(DIRS)], function(d)
  list.files(d, pattern = "[.]csv$", full.names = TRUE)))
files <- files[!grepl("ALL_CURRENT_VALUES", files)]
message("scanning ", length(files), " csv files")

one <- function(f) {
  d <- tryCatch(suppressWarnings(
         readr::read_csv(f, show_col_types = FALSE, progress = FALSE)),
       error = function(e) NULL)
  if (is.null(d) || !nrow(d) || !ncol(d)) return(NULL)

  is_num <- vapply(d, is.numeric, logical(1))
  # label columns: everything non-numeric, plus any numeric column that looks
  # like an identifier (year, tercile, decile, bin, group) rather than a result
  idlike <- grepl("^(year|tercile|decile|bin|group|dec|res_decile|quadrant)$",
                  tolower(names(d)))
  lab <- names(d)[!is_num | idlike]
  val <- names(d)[is_num & !idlike]
  if (!length(val)) return(NULL)

  key <- if (length(lab))
    apply(d[lab], 1, function(r) paste(trimws(as.character(r)), collapse = " | "))
  else as.character(seq_len(nrow(d)))

  out <- lapply(val, function(v)
    tibble(file = basename(f), row_key = key, column = v,
           value = suppressWarnings(as.numeric(d[[v]]))))
  bind_rows(out) |> filter(is.finite(value))
}

all <- bind_rows(lapply(files, one))
message("collected ", format(nrow(all), big.mark = ","), " numeric values")

# a rounded copy makes text-matching against the manuscript practical
all <- all |>
  mutate(r2 = round(value, 2), r3 = round(value, 3), r4 = round(value, 4))

out_f <- file.path(CO_DIR, "output/models/ALL_CURRENT_VALUES.csv")
write.csv(all, out_f, row.names = FALSE)
message("wrote ", out_f)

## ---- a short readout of the quantities the manuscript quotes most -------------
show <- function(lbl, f, pat_key = NULL, pat_col = NULL) {
  s <- all |> filter(file == f)
  if (!is.null(pat_key)) s <- s |> filter(grepl(pat_key, row_key))
  if (!is.null(pat_col)) s <- s |> filter(grepl(pat_col, column))
  if (!nrow(s)) { message("  [", lbl, "] not found in ", f); return(invisible()) }
  message("\n-- ", lbl)
  print(as.data.frame(s |> transmute(row_key = substr(row_key, 1, 58),
                                     column, value = round(value, 4))),
        row.names = FALSE)
}

message("\n", strrep("=", 78))
message("THE QUANTITIES THE MAIN TEXT QUOTES")
message(strrep("=", 78))

show("accessibility ladder (zoning attenuation)", "p4_accessibility_ladder.csv",
     "pct_reslow_of_res")
show("network-accessibility ladder (Table S9)", "p4_network_access_ladder.csv",
     ":z_")
show("spatial inference (Moran, Conley)", "p4_spatial_inference.csv")
show("2023 descriptives by tercile", "co_descriptives_2023.csv")
show("opportunity vs realized (92/8)", "p4_opportunity_vs_realized.csv",
     "exclusionary_tercile")
show("measure sensitivity (decay, truncation)", "p4_measure_sensitivity.csv",
     "preferred")
show("alternative dimensions", "p4_alt_dimensions.csv")
show("tercile descriptives", "p4_tercile_descriptives.csv")

message("\n", strrep("=", 78))
message("Stage ALL_CURRENT_VALUES.csv to check every manuscript number against")
message("the pipeline's current output, rather than only the ones we expect to")
message("have moved.")
message(strrep("=", 78))
