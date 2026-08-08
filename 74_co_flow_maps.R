# ==============================================================================
# 74_co_flow_maps.R      [PAPER 4, step 14]
# Commute-flow visualizations: where the most- and least-segregated
# neighborhoods' residents actually go to work.
#
#   Fig F1  Paired commute-shed maps: origins = top vs bottom DECILE of
#           residential D (in-scope Denver frame, 2023); curved flow lines
#           from home to workplace tract, colored by the DESTINATION's
#           workplace-location segregation quintile, width/alpha ~ flow.
#           Verified headline annotation: flow-weighted destination D =
#           0.018 (least-segregated origins) vs 0.034 (most) -- nearly 2x.
#
#   Fig F2  The paper's thesis as a picture: two tracts with near-identical
#           residential D whose destination fans differ completely.
#           Verified exemplar pair (n_commuters >= 500, res D within .004,
#           max wexp gap): 08005006864 (Arapahoe; res D .166, wexp .021)
#           vs 08059012058 (JeffCo foothills; res D .168, wexp .063 -- 3x).
#           The pair is re-derived from the data below, so it updates
#           automatically if the panel changes.
#
# Needs: 62's OD caches, 63's panel, Paper 2 seg panel, centroids, 66's WAC
#        cache (employment centers), tract geometry (61).
# Output: output/figures/p4_figF1_commute_sheds.png
#         output/figures/p4_figF2_exemplar_pair.png
# ==============================================================================

source("60_co_setup.R")
library(patchwork)

## ---- pieces ------------------------------------------------------------------
fr <- readRDS(file.path(DIR_CO_OUT, "co_analysis_home_panel.rds")) |>
  filter(year == CO_ANCHOR_YEAR, in_scope, denver_msa,
         n_commuters >= P3_MIN_COMMUTERS)
cent <- readRDS(file.path(DIR_CLEAN, "tract_centroids_km.rds")) |>
  transmute(tract_id = as.character(GEOID), CBSA_Code, X = X_km * 1000,
            Y = Y_km * 1000)
seg23 <- readRDS(file.path(DIR_P2_OUT, "p2_tract_segregation_panel.rds")) |>
  filter(year == CO_ANCHOR_YEAR) |>
  transmute(tract_id, d_wac = d_whiteblack_wac_half)
geom <- readRDS(file.path(DIR_CO_CLEAN, "co_tract_geom.rds"))

od <- map(P3_OD_PARTS, function(part) {
  for (f in c(p3_od_cache(part, CO_ANCHOR_YEAR, "co"),
              co_od_cache(part, CO_ANCHOR_YEAR)))
    if (file.exists(f)) return(readRDS(f))
  NULL
}) |> compact() |> bind_rows() |>
  group_by(w_tract, h_tract) |>
  summarise(S000 = sum(S000, na.rm = TRUE), .groups = "drop") |>
  inner_join(cent, by = c("h_tract" = "tract_id")) |>
  rename(cbsa_h = CBSA_Code, xh = X, yh = Y) |>
  inner_join(cent, by = c("w_tract" = "tract_id")) |>
  rename(cbsa_w = CBSA_Code, xw = X, yw = Y) |>
  filter(cbsa_h == CO_CBSA_MAIN, cbsa_w == CO_CBSA_MAIN) |>
  left_join(seg23, by = c("w_tract" = "tract_id"))

wac <- readRDS(file.path(DIR_CO_CLEAN,
                         sprintf("co_wac_tract_%d.rds", CO_ANCHOR_YEAR))) |>
  transmute(tract_id = as.character(tract_id), jobs = C000) |>
  inner_join(cent, by = "tract_id") |>
  filter(CBSA_Code %in% CO_CBSA_KEEP)
ctr <- wac |> filter(jobs >= quantile(jobs, 0.98))

# extent: urbanised core (same convention as 70)
cxy <- suppressWarnings(st_coordinates(st_centroid(
  geom |> inner_join(fr |> select(tract_id), by = "tract_id"))))
xlim <- quantile(cxy[, 1], c(.01, .93)) + c(-7000, 7000)
ylim <- quantile(cxy[, 2], c(.02, .98)) + c(-7000, 7000)
BLUES <- c("#c6dbef", "#9ecae1", "#6baed6", "#3182bd", "#08519c")

## ---- basemap anchors ---------------------------------------------------------
counties_sf <- geom |>
  mutate(county = substr(tract_id, 1, 5)) |>
  group_by(county) |> summarise(.groups = "drop")
jlab <- fr |>
  inner_join(tibble(tract_id = geom$tract_id,
                    cx = st_coordinates(suppressWarnings(st_centroid(geom)))[, 1],
                    cy = st_coordinates(suppressWarnings(st_centroid(geom)))[, 2]),
             by = "tract_id") |>
  mutate(name = c(DENVER = "Denver", AURORA = "Aurora", Lakewood = "Lakewood",
                  Arvada = "Arvada", THORNTON = "Thornton",
                  WESTMINSTER = "Westminster", CENTENNIAL = "Centennial",
                  LITTLETON = "Littleton", Parker = "Parker",
                  Broomfield = "Broomfield",
                  ENGLEWOOD = "Englewood")[jurisd_main]) |>
  filter(!is.na(name)) |>
  group_by(name) |>
  summarise(n = n(), lx = mean(cx), ly = mean(cy), .groups = "drop") |>
  filter(n >= 8)
cbd <- st_sfc(st_point(c(-105.0000, 39.7531)), crs = 4326) |>
  st_transform(CRS_METERS) |> st_coordinates()
roads <- tryCatch(
  tigris::primary_secondary_roads("CO", year = 2023) |>
    st_transform(CRS_METERS) |> filter(RTTYP %in% c("I", "U")),
  error = function(e) { message("roads skipped: ", conditionMessage(e)); NULL })

base_layers <- function() {
  out <- list(
    geom_sf(data = geom, fill = "white", color = "#eceae4", linewidth = .18),
    geom_sf(data = counties_sf, fill = NA, color = "#b9b8b0", linewidth = .55))
  if (!is.null(roads))
    out <- c(out, list(geom_sf(data = roads, color = "#d8d0c2",
                               linewidth = .45)))
  c(out, list(
    geom_point(data = ctr, aes(X, Y), shape = 24, size = 2.4, fill = "white",
               color = "#52514e", stroke = .8),
    annotate("point", x = cbd[1], y = cbd[2], shape = 8, size = 3.2,
             color = "#0b0b0b"),
    annotate("text", x = cbd[1] + 2500, y = cbd[2] + 2500,
             label = "downtown Denver", fontface = "bold", size = 2.6,
             hjust = 0),
    geom_text(data = jlab, aes(lx, ly, label = name), size = 2.7,
              color = "#52514e")))
}

theme_flow <- theme_void(base_size = 10) +
  theme(plot.title = element_text(size = 11, hjust = 0),
        plot.subtitle = element_text(size = 8, color = "#52514e"),
        legend.position = "bottom")

shed_panel <- function(origins, title) {
  f_all <- od |> filter(h_tract %in% origins, !is.na(d_wac))
  fw <- with(f_all, sum(d_wac * S000) / sum(S000))  # statistic on ALL flows,
  # matching the paper's wexp construct (same-tract commutes included);
  # arcs below drop same-tract flows only because geom_curve cannot draw them
  f <- f_all |> filter(S000 >= 5, h_tract != w_tract) |>
    mutate(q = ntile(d_wac, 5), w = S000 / max(S000)) |>
    arrange(S000)
  ggplot() +
    base_layers() +
    geom_sf(data = geom |> filter(tract_id %in% origins),
            fill = "#efc9b8", alpha = .40, color = "#c8825f",
            linewidth = .3) +
    geom_curve(data = f, aes(x = xh, y = yh, xend = xw, yend = yw,
                             color = factor(q), linewidth = w,
                             alpha = pmin(.75, .12 + .5 * sqrt(w))),
               curvature = 0.18) +
    scale_color_manual(values = BLUES,
                       labels = c("least segregated", "2nd", "3rd", "4th",
                                  "most segregated"),
                       name = "destination workplace-location D (quintile)") +
    scale_linewidth(range = c(.12, 1.1), guide = "none") +
    scale_alpha_identity() +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    labs(title = title,
         subtitle = sprintf(
           "origins shaded; arcs: flows >= 5 workers between distinct tracts (statistic includes all flows) | flow-weighted destination D = %.3f",
           fw)) +
    theme_flow
}

fr_s <- fr |> arrange(d_whiteblack_rac_half)
n10  <- ceiling(nrow(fr_s) / 10)
pF1 <- (shed_panel(head(fr_s, n10)$tract_id,
                   "A. Least-segregated decile of neighborhoods") |
        shed_panel(tail(fr_s, n10)$tract_id,
                   "B. Most-segregated decile of neighborhoods")) +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave(file.path(DIR_CO_FIG, "p4_figF1_commute_sheds.png"), pF1,
       width = 13.4, height = 7, dpi = 350, bg = "white")

## ---- F2: exemplar pair (re-derived, not hard-coded) --------------------------
cand <- fr |>
  filter(n_commuters >= 500,
         d_whiteblack_rac_half >=
           quantile(d_whiteblack_rac_half, .6, na.rm = TRUE)) |>
  select(tract_id, res = d_whiteblack_rac_half,
         wexp = wexp_whiteblack_wac_half)
pairs <- tidyr::crossing(a = seq_len(nrow(cand)), b = seq_len(nrow(cand))) |>
  filter(a < b) |>
  mutate(dres = abs(cand$res[a] - cand$res[b]),
         dwex = abs(cand$wexp[a] - cand$wexp[b])) |>
  filter(dres < 0.004) |> arrange(desc(dwex)) |> slice(1)
tA <- cand$tract_id[pairs$a]; tB <- cand$tract_id[pairs$b]
if (cand$wexp[pairs$a] > cand$wexp[pairs$b]) { tmp <- tA; tA <- tB; tB <- tmp }
sA <- cand |> filter(tract_id == tA); sB <- cand |> filter(tract_id == tB)
message(sprintf("Exemplar pair: A %s (res %.3f, wexp %.4f) vs B %s (res %.3f, wexp %.4f)",
                tA, sA$res, sA$wexp, tB, sB$res, sB$wexp))

fan <- function(t, col) {
  f <- od |> filter(h_tract == t, S000 >= 3, h_tract != w_tract) |>
    mutate(w = S000 / max(S000)) |> arrange(S000)
  geom_curve(data = f, aes(x = xh, y = yh, xend = xw, yend = yw,
                           linewidth = w,
                           alpha = pmin(.8, .18 + .55 * sqrt(w))),
             color = col, curvature = 0.18)
}
pF2 <- ggplot() +
  base_layers() +
  fan(tA, "#2a78d6") + fan(tB, "#eb6834") +
  geom_sf(data = geom |> filter(tract_id == tA), fill = "#2a78d6",
          color = "#0b0b0b", alpha = .85, linewidth = .5) +
  geom_sf(data = geom |> filter(tract_id == tB), fill = "#eb6834",
          color = "#0b0b0b", alpha = .85, linewidth = .5) +
  scale_linewidth(range = c(.15, 1.5), guide = "none") +
  scale_alpha_identity() +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  labs(title = "Same residential segregation, different labor markets",
       subtitle = sprintf(paste(
         "Tract A (blue): residential D %.3f -> exposure %.3f | Tract B",
         "(orange): residential D %.3f -> exposure %.3f | triangles =",
         "major employment centers"),
         sA$res, sA$wexp, sB$res, sB$wexp)) +
  theme_flow
ggsave(file.path(DIR_CO_FIG, "p4_figF2_exemplar_pair.png"), pF2,
       width = 9.4, height = 9.8, dpi = 350, bg = "white")

message("74 complete. Flow maps in ", DIR_CO_FIG)
