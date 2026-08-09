# ==============================================================================
# 80_co_decentralization_network.R  [PAPER 4, step 20 -- SI extensions 7 + 8]
#
#  Ext. 7 -- EMPLOYMENT DECENTRALIZATION 2011-2023. Where did Denver's jobs
#    move over the LODES record, and does the spatial structure underlying the
#    position result look stable, strengthening, or weakening? Per year:
#    share of MSA jobs within 5 km of the CBD, share in the 2023-defined major
#    employment centers (top 2% of tracts by 2023 jobs -- held FIXED so the
#    series measures movement of jobs, not movement of the definition),
#    job-weighted mean distance of employment from the CBD, and the
#    flow-weighted mean distance of commuting DESTINATIONS from the CBD.
#    Descriptive: there is no historical zoning; 2020-21 carry the WFH flag.
#
#  Ext. 8 -- LABOR-MARKET CATCHMENTS (light network structure, no igraph
#    dependency). Treat tracts as nodes and flows as edges; assign each origin
#    its PRIMARY destination jurisdiction (largest flow share). The map of
#    primary destinations reveals labor-market subregions that administrative
#    boundaries miss; the table overlays zoning (mean exclusionary share of
#    each catchment). Full community detection stays in the separate
#    network paper; this is the descriptive core.
#
# Output: output/models/p4_decentralization_series.csv
#         output/models/p4_network_catchments.csv
#         output/figures/p4_fig_decentralization.png
#         output/figures/p4_fig_catchments.png
# ==============================================================================

source("60_co_setup.R")

## ---- shared ------------------------------------------------------------------
dat <- readRDS(file.path(DIR_CO_OUT, "co_analysis_home_panel.rds"))
xs <- dat |>
  filter(year == CO_ANCHOR_YEAR, in_scope, denver_msa,
         n_commuters >= P3_MIN_COMMUTERS) |>
  select(tract_id, jurisd_main, pct_reslow_of_res) |>
  mutate(tercile = ntile(pct_reslow_of_res, 3))
acc <- readRDS(file.path(DIR_CO_CLEAN, "co_accessibility_2023.rds")) |>
  select(tract_id, dist_cbd_km)
cent <- readRDS(file.path(DIR_CLEAN, "tract_centroids_km.rds")) |>
  transmute(tract_id = as.character(GEOID), CBSA_Code)
msa_tracts <- cent |> filter(CBSA_Code == CO_CBSA_MAIN) |> pull(tract_id)

grab_wac <- function(yr) {
  cache <- file.path(DIR_CO_CLEAN, sprintf("co_wac_tract_%s.rds", yr))
  if (file.exists(cache)) return(readRDS(cache))
  message("Downloading WAC tract co ", yr)
  df <- grab_lodes(state = "co", year = yr, version = "LODES8",
                   lodes_type = "wac", job_type = "JT01", agg_geo = "tract") |>
    mutate(tract_id = as.character(w_tract)) |>
    select(tract_id, any_of(c("C000", sprintf("CE%02d", 1:3),
                              sprintf("CD%02d", 1:4), sprintf("CNS%02d", 1:20))))
  saveRDS(df, cache)
  df
}
read_od_yr <- function(yr) {
  map(P3_OD_PARTS, function(part) {
    for (f in c(p3_od_cache(part, yr, "co"), co_od_cache(part, yr)))
      if (file.exists(f)) return(readRDS(f))
    NULL
  }) |> compact() |> bind_rows() |>
    group_by(w_tract, h_tract) |>
    summarise(S000 = sum(S000, na.rm = TRUE), .groups = "drop")
}

## ---- ext. 7: the decentralization series -------------------------------------
wac23 <- grab_wac(CO_ANCHOR_YEAR)
ctr23 <- wac23 |>
  filter(tract_id %in% msa_tracts) |>
  filter(C000 >= quantile(C000, 0.98, na.rm = TRUE)) |>
  pull(tract_id)
message(length(ctr23), " employment-center tracts (2023 definition, held fixed)")

series <- map(CO_YEARS, function(yr) {
  w <- grab_wac(yr) |>
    filter(tract_id %in% msa_tracts) |>
    left_join(acc, by = "tract_id") |>
    filter(!is.na(dist_cbd_km))
  o <- read_od_yr(yr)
  dest_cbd <- if (!is.null(o) && nrow(o)) {
    o |> filter(h_tract %in% msa_tracts, w_tract %in% msa_tracts) |>
      left_join(acc, by = c(w_tract = "tract_id")) |>
      summarise(v = weighted.mean(dist_cbd_km, S000, na.rm = TRUE)) |> pull(v)
  } else NA_real_
  tibble(year = yr,
         jobs_msa = sum(w$C000),
         pct_jobs_within_5km_cbd = 100 * sum(w$C000[w$dist_cbd_km <= 5]) /
                                         sum(w$C000),
         pct_jobs_in_2023_centers = 100 * sum(w$C000[w$tract_id %in% ctr23]) /
                                          sum(w$C000),
         jobwt_mean_cbd_km = weighted.mean(w$dist_cbd_km, w$C000),
         flowwt_dest_cbd_km = dest_cbd,
         pandemic_flag = yr %in% c(2020L, 2021L))
}) |> bind_rows()
write.csv(series, file.path(DIR_CO_MOD, "p4_decentralization_series.csv"),
          row.names = FALSE)
print(as.data.frame(series), digits = 4)

pD <- series |>
  select(year, pandemic_flag,
         `% of jobs within 5 km of CBD` = pct_jobs_within_5km_cbd,
         `% of jobs in 2023-defined centers` = pct_jobs_in_2023_centers,
         `job-weighted mean km from CBD` = jobwt_mean_cbd_km,
         `flow-weighted destination km from CBD` = flowwt_dest_cbd_km) |>
  pivot_longer(-c(year, pandemic_flag)) |>
  ggplot(aes(year, value)) +
  geom_line(color = "#08519c", linewidth = .7) +
  geom_point(aes(shape = pandemic_flag), color = "#08519c", size = 1.8) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1),
                     labels = c("normal year", "pandemic (WFH affects LODES)"),
                     name = NULL) +
  facet_wrap(~name, scales = "free_y") +
  scale_x_continuous(breaks = seq(2011, 2023, 4)) +
  labs(title = "Employment centralization in the Denver MSA, 2011-2023",
       subtitle = paste("Employment-center set fixed at its 2023 definition;",
                        "descriptive -- no historical zoning exists."),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 9.5) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 8))
ggsave(file.path(DIR_CO_FIG, "p4_fig_decentralization.png"), pD,
       width = 8.6, height = 5.6, dpi = 350, bg = "white")

## ---- ext. 8: labor-market catchments -----------------------------------------
jmap <- xs |> select(tract_id, jurisd_main)
od23 <- read_od_yr(CO_ANCHOR_YEAR) |>
  inner_join(xs |> select(tract_id, tercile, pct_reslow_of_res),
             by = c(h_tract = "tract_id")) |>
  inner_join(jmap |> rename(w_jur = jurisd_main), by = c(w_tract = "tract_id"))
catch <- od23 |>
  group_by(h_tract, tercile, pct_reslow_of_res) |>
  summarise(flows = sum(S000),
            primary_jur = names(which.max(tapply(S000, w_jur, sum))),
            pct_primary = 100 * max(tapply(S000, w_jur, sum)) / sum(S000),
            primary_nondenver = {
              s <- tapply(S000[w_jur != "DENVER"], w_jur[w_jur != "DENVER"], sum)
              if (length(s)) names(which.max(s)) else NA_character_ },
            pct_nondenver_primary = {
              s <- tapply(S000[w_jur != "DENVER"], w_jur[w_jur != "DENVER"], sum)
              if (length(s)) 100 * max(s) / sum(S000) else NA_real_ },
            .groups = "drop")
message(sprintf(
  "Monocentricity at the flow level: Denver is the primary destination for %d of %d origin tracts (%.0f%%).",
  sum(catch$primary_jur == "DENVER"), nrow(catch),
  100 * mean(catch$primary_jur == "DENVER")))
catch_tab <- catch |>
  group_by(primary_jur) |>
  summarise(n_origin_tracts = n(),
            commuters = sum(flows),
            mean_pct_primary = weighted.mean(pct_primary, flows),
            mean_pct_reslow_origins = weighted.mean(pct_reslow_of_res, flows),
            .groups = "drop") |>
  arrange(desc(commuters))
write.csv(catch_tab, file.path(DIR_CO_MOD, "p4_network_catchments.csv"),
          row.names = FALSE)
print(as.data.frame(head(catch_tab, 12)), digits = 3)

# Denver is the primary destination almost everywhere (monocentricity), so the
# informative map is the SECOND-ORDER structure: the primary NON-Denver
# destination, which reveals the suburban labor-market subregions.
nd_tab <- catch |>
  filter(!is.na(primary_nondenver)) |>
  count(primary_nondenver, wt = flows, sort = TRUE)
geom <- readRDS(file.path(DIR_CO_CLEAN, "co_tract_geom.rds"))
top8 <- head(nd_tab$primary_nondenver, 8)
map_df <- geom |>
  inner_join(catch |>
               filter(!is.na(primary_nondenver)) |>
               mutate(shown = ifelse(primary_nondenver %in% top8,
                                     primary_nondenver, "other")),
             by = c("tract_id" = "h_tract"))
counties_sf <- geom |>
  mutate(county = substr(tract_id, 1, 5)) |>
  group_by(county) |> summarise(.groups = "drop")
cxy <- suppressWarnings(st_coordinates(st_centroid(map_df)))
xlim <- quantile(cxy[, 1], c(.01, .99)) + c(-5000, 5000)
ylim <- quantile(cxy[, 2], c(.01, .99)) + c(-5000, 5000)
pal <- c("#08519c", "#e6550d", "#31a354", "#756bb1", "#c51b8a",
         "#6baed6", "#fdae6b", "#74c476", "grey80")
names(pal) <- c(top8, "other")
pC <- ggplot(map_df) +
  geom_sf(aes(fill = shown), color = "white", linewidth = .05) +
  geom_sf(data = counties_sf, fill = NA, color = "grey35", linewidth = .45) +
  scale_fill_manual(values = pal,
                    name = "primary non-Denver destination") +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  labs(title = "Second-order labor-market catchments",
       subtitle = paste(
         "Denver is the primary destination for nearly every tract",
         "(flow-level monocentricity);\neach origin tract is colored by its",
         "largest NON-Denver destination jurisdiction, 2023.")) +
  theme_void(base_size = 9.5) +
  theme(legend.position = "right",
        legend.key.size = unit(0.35, "cm"),
        legend.text = element_text(size = 7))
ggsave(file.path(DIR_CO_FIG, "p4_fig_catchments.png"), pC,
       width = 8.8, height = 7.2, dpi = 350, bg = "white")
message("80 complete.")
