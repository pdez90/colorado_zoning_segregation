# ==============================================================================
# 91_validate_all_outputs.R
# Compare EVERY committed result CSV against the one this machine produced.
#
# 90_validate_rebuild.R checks a hand-picked set of quantities in detail. This
# sweeps the rest: for each file in the reference's output/models, it lines the
# two copies up on whatever non-numeric columns they share and reports the
# largest difference across every numeric column, plus anything not yet
# produced.
#
#   bash run_all.sh validate        # runs 90 then this
#   # or, directly:
#   Sys.setenv(REF_DIR = "~/Downloads/LODES/colorado_zoning_segregation")
#   source("91_validate_all_outputs.R")
#
# Writes diagnostics/91_all_outputs_validation.csv.
#
# Deliberately base R only -- no dplyr, readr or purrr. The comparison logic is
# the thing being trusted here, so it is written so it can be exercised
# anywhere, including outside this pipeline.
#
# HOW TO READ IT
#   match         every numeric column agrees within tolerance on every row
#   DIFFERS       names the worst column and the gap -- the thing to chase
#   PARTIAL       fewer than 90% of reference rows found a partner
#   SHAPE         could not be lined up: different row counts, no usable key.
#                 Usually an output-format change, not a numbers change
#   NOT PRODUCED  the script that writes it has not been run here
# ==============================================================================

CO_DIR <- path.expand(Sys.getenv("CO_DIR", "~/Downloads/LODES/Colorado"))
if (!exists("DIR_CO")) {
  if (!file.exists(file.path(CO_DIR, "60_co_setup.R")))
    stop("60_co_setup.R not found in ", CO_DIR)
  .owd <- setwd(CO_DIR); source("60_co_setup.R"); setwd(.owd)
}

REF <- path.expand(Sys.getenv("REF_DIR",
                              "~/Downloads/LODES/colorado_zoning_segregation"))
if (!dir.exists(REF)) stop("Reference clone not found at ", REF)

TOL <- 5e-3

## ---- the comparison, base R --------------------------------------------------
read_csv_safe <- function(p) {
  if (!file.exists(p)) return(NULL)
  tryCatch(utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE),
           error = function(e) NULL)
}

compare_tables <- function(a, b, tol = TOL) {
  shared <- intersect(names(a), names(b))
  if (!length(shared))
    return(list(n_joined = NA_integer_, worst_col = NA_character_,
                max_abs_diff = NA_real_, verdict = "SHAPE"))

  is_num <- function(d, cols) cols[vapply(d[cols], is.numeric, logical(1))]
  nums <- intersect(is_num(a, shared), is_num(b, shared))
  keys <- setdiff(shared, nums)

  if (!length(nums))
    return(list(n_joined = NA_integer_, worst_col = NA_character_,
                max_abs_diff = NA_real_, verdict = "no numeric columns"))

  j <- NULL
  if (length(keys)) {
    ka <- a[keys]; kb <- b[keys]
    # a non-unique key turns the merge into a cross product and the report
    # into nonsense; fall back to position instead
    if (anyDuplicated(ka) == 0L && anyDuplicated(kb) == 0L) {
      j <- merge(a[c(keys, nums)], b[c(keys, nums)], by = keys,
                 suffixes = c(".ref", ".new"))
      if (nrow(j) == 0L) j <- NULL
    }
  }
  if (is.null(j) && nrow(a) == nrow(b)) {
    j <- cbind(stats::setNames(a[nums], paste0(nums, ".ref")),
               stats::setNames(b[nums], paste0(nums, ".new")))
  }
  if (is.null(j))
    return(list(n_joined = NA_integer_, worst_col = NA_character_,
                max_abs_diff = NA_real_, verdict = "SHAPE"))

  worst <- 0; wcol <- NA_character_
  for (v in nums) {
    x <- suppressWarnings(as.numeric(j[[paste0(v, ".ref")]]))
    y <- suppressWarnings(as.numeric(j[[paste0(v, ".new")]]))
    ok <- is.finite(x) & is.finite(y)
    if (!any(ok)) next
    d <- max(abs(x[ok] - y[ok]))
    if (is.finite(d) && d > worst) { worst <- d; wcol <- v }
  }
  verdict <- if (nrow(j) < nrow(a) * 0.9) "PARTIAL" else
             if (worst <= tol) "match" else "DIFFERS"
  list(n_joined = nrow(j), worst_col = wcol, max_abs_diff = worst,
       verdict = verdict)
}

compare_file <- function(fname) {
  a <- read_csv_safe(file.path(REF,    "output/models", fname))
  b <- read_csv_safe(file.path(DIR_CO, "output/models", fname))
  if (is.null(a)) return(NULL)
  if (is.null(b))
    return(data.frame(file = fname, n_ref = nrow(a), n_new = NA_integer_,
                      n_joined = NA_integer_, worst_col = NA_character_,
                      max_abs_diff = NA_real_, verdict = "NOT PRODUCED",
                      stringsAsFactors = FALSE))
  r <- compare_tables(a, b)
  data.frame(file = fname, n_ref = nrow(a), n_new = nrow(b),
             n_joined = r$n_joined, worst_col = r$worst_col,
             max_abs_diff = r$max_abs_diff, verdict = r$verdict,
             stringsAsFactors = FALSE)
}

## ---- run ---------------------------------------------------------------------
files <- sort(basename(list.files(file.path(REF, "output/models"),
                                  pattern = "\\.csv$")))
message("comparing ", length(files), " committed result files against ",
        file.path(DIR_CO, "output/models"), "\n")

res <- do.call(rbind, lapply(files, compare_file))
ord <- c("DIFFERS", "PARTIAL", "SHAPE", "NOT PRODUCED", "no numeric columns",
         "match")
res <- res[order(match(res$verdict, ord), -replace(res$max_abs_diff,
                                                   is.na(res$max_abs_diff), 0)), ]
shown <- res; shown$max_abs_diff <- signif(shown$max_abs_diff, 4)
print(shown, row.names = FALSE)

message("\n== SUMMARY ==")
print(as.data.frame(table(verdict = res$verdict)), row.names = FALSE)

np <- res[res$verdict == "NOT PRODUCED", ]
if (nrow(np))
  message("\n", nrow(np), " file(s) not produced yet. Run scripts 65-87:\n",
          "  bash run_all.sh rest")

bad <- res[res$verdict %in% c("DIFFERS", "PARTIAL", "SHAPE"), ]
if (nrow(bad)) {
  message("\nneeds attention:")
  print(bad[, c("file", "worst_col", "max_abs_diff", "verdict")],
        row.names = FALSE)
} else if (!nrow(np)) {
  message("\nEvery committed result file reproduces within ", TOL, ".")
}

write_codiag(res, "91_all_outputs_validation")
message("91_validate_all_outputs.R complete.")
