# ==============================================================================
# patch_60_paths.R
# 60_co_setup.R hard-codes two input paths from the old machine. This finds
# where they actually live now and rewrites those two lines in place, keeping
# a .bak copy. Run it once, from anywhere:
#
#   source("~/Downloads/LODES/paper_pipeline/patch_60_paths.R")
#
# It changes nothing else in the file, and it refuses to write a path that
# does not exist -- so if it reports a miss, find the file and add its folder
# to the candidate list below rather than editing 60 by hand.
# ==============================================================================

f60 <- path.expand("~/Downloads/LODES/Colorado/60_co_setup.R")
if (!file.exists(f60)) stop("60_co_setup.R not found at ", f60)

## ---- 1. the zoning shapefile -------------------------------------------------
ZON_NAME <- "ALL_DenverMSA_10.3.23.shp"
ZON_CANDIDATES <- path.expand(c(
  "~/Wellbeing/data/Zoning",      # where it is on this machine
  "~/Wellbeing/Zoning",           # where 60 was written to look
  "~/Downloads",
  "~/Downloads/Zoning",
  "~/Downloads/LODES/Colorado/zoning"))

zon_hit <- NA_character_
for (d in ZON_CANDIDATES) {
  p <- file.path(d, ZON_NAME)
  if (file.exists(p)) { zon_hit <- p; break }
}
if (is.na(zon_hit)) {                      # last resort: search the home dir
  hits <- list.files(path.expand("~"), pattern = paste0("^", ZON_NAME, "$"),
                     recursive = TRUE, full.names = TRUE)
  if (length(hits)) zon_hit <- hits[1]
}

## ---- 2. the TIGER tract folder -----------------------------------------------
TIG_CANDIDATES <- path.expand(c(
  "~/Downloads/LODES/TIGER2024_TRACT_UNZIPPED/tl_2024_08_tract",   # 10 writes here
  "~/Downloads/tl_2024_08_tract",
  "~/Downloads/LODES/Colorado/tl_2024_08_tract"))
tig_hit <- TIG_CANDIDATES[file.exists(file.path(TIG_CANDIDATES,
                                                "tl_2024_08_tract.shp"))][1]

## ---- 3. rewrite --------------------------------------------------------------
src <- readLines(f60, warn = FALSE)
orig <- src
esc <- function(p) gsub("\\\\", "/", p)

if (!is.na(zon_hit)) {
  i <- grep("^CO_ZONING_SHP\\s*<-", src)
  if (length(i) == 1) {
    src[i] <- sprintf('CO_ZONING_SHP <- path.expand("%s")',
                      sub(path.expand("~"), "~", esc(zon_hit), fixed = TRUE))
    message("  CO_ZONING_SHP -> ", zon_hit)
  } else message("  !! could not find a single CO_ZONING_SHP line")
} else message("  !! zoning shapefile not found -- 61 and 70 will stop. ",
               "Locate ", ZON_NAME, " and add its folder to ZON_CANDIDATES.")

if (!is.na(tig_hit)) {
  i <- grep("^CO_TIGER_DIR\\s*<-", src)
  if (length(i) == 1) {
    # the assignment spans two lines in the original; drop the continuation
    j <- if (i < length(src) && grepl("^\\s*\"", src[i + 1])) i + 1 else i
    src <- c(src[seq_len(i - 1)],
             sprintf('CO_TIGER_DIR  <- path.expand("%s")',
                     sub(path.expand("~"), "~", esc(tig_hit), fixed = TRUE)),
             if (j < length(src)) src[(j + 1):length(src)])
    message("  CO_TIGER_DIR  -> ", tig_hit)
  } else message("  !! could not find a single CO_TIGER_DIR line")
} else message("  note: no TIGER folder yet -- run 10_geography.R, which ",
               "writes one, then re-run this patch. (61 falls back to ",
               "tigris' 2023 vintage otherwise, which is NOT what the ",
               "published zoning measures used.)")

if (identical(src, orig)) {
  message("Nothing to change.")
} else {
  file.copy(f60, paste0(f60, ".bak"), overwrite = TRUE)
  writeLines(src, f60)
  message("Patched ", f60, " (backup at ", basename(f60), ".bak)")
}
