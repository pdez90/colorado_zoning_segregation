# ==============================================================================
# 65_co_figures.R      [COLORADO ZONING case study, step 5]
# Figures (Denver MSA, 2023 anchor unless noted):
#   Fig 1  three-panel tract map: exclusionary zoning share | residential D |
#          OD-weighted workplace-segregation exposure
#   Fig 2  binned scatter: res seg vs wexp, colored by exclusionary-zoning
#          tercile (the zoning analogue of Paper 3's Fig 2)
#   Fig 3  coefficient forest: ZB interaction terms (headline moderation)
#   Fig 4  jurisdiction scatter: % Residential_Low vs mean residential D,
#          point size = resident workers (zoning is WRITTEN here)
#   Fig 5  trajectories: 2011->2023 change in res D by exclusionary tercile
#
# Needs: 61 (geometry + jurisdiction summary), 63 (panels), 64 (coefficients).
# Output: output/figures/co_fig*.png
# ==============================================================================

source("60_co_setup.R")
library(patchwork)

dat  <- readRDS(file.path(DIR_CO_OUT, "co_analysis_home_panel.rds"))
geom <- readRDS(file.path(DIR_CO_CLEAN, "co_tract_geom.rds"))
jur  <- readRDS(file.path(DIR_CO_CLEAN, "co_zoning_jurisd.rds"))
coefs <- read.csv(file.path(DIR_CO_MOD, "co_model_coefficients.csv"))

theme_co <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

xs <- dat |>
  filter(year == CO_ANCHOR_YEAR, in_scope, denver_msa,
         n_commuters >= P3_MIN_COMMUTERS) |>
  mutate(reslow_tercile = ntile(pct_res_low, 3))

## ---- Fig 1: three-panel map --------------------------------------------------
gmap <- geom |> inner_join(
  xs |> select(tract_id, pct_res_low, d_whiteblack_rac_half,
               wexp_whiteblack_wac_half),
  by = "tract_id")
panel <- function(var, ttl, lab) {
  ggplot(gmap) +
    geom_sf(aes(fill = .data[[var]]), color = NA) +
    scale_fill_viridis_c(name = lab) +
    labs(title = ttl) + theme_void(base_size = 10) +
    theme(legend.position = "bottom", legend.key.height = unit(3, "mm"))
}
p1 <- panel("pct_res_low", "Exclusionary zoning", "% zoned Res_Low") +
  panel("d_whiteblack_rac_half", "Residential segregation",
        "D (RAC, beta=0.5)") +
  panel("wexp_whiteblack_wac_half", "Workplace-seg exposure",
        "OD-weighted D (WAC)")
ggsave(file.path(DIR_CO_FIG, "co_fig1_maps.png"), p1,
       width = 13, height = 5.5, dpi = 300, bg = "white")

## ---- Fig 2: binned scatter by exclusionary tercile ---------------------------
p2 <- xs |>
  mutate(bin = ntile(d_whiteblack_rac_half, 20)) |>
  group_by(reslow_tercile, bin) |>
  summarise(res = mean(d_whiteblack_rac_half),
            wexp = mean(wexp_whiteblack_wac_half), .groups = "drop") |>
  ggplot(aes(res, wexp, color = factor(reslow_tercile))) +
  geom_point(size = 1.2) +
  geom_smooth(method = "lm", se = FALSE, linewidth = .6) +
  scale_color_manual(values = c("#1b6ca8", "#888888", "#c8442c"),
                     labels = c("least", "middle", "most"),
                     name = "Exclusionary-zoning\ntercile (% Res_Low)") +
  labs(title = "Residential segregation vs. workplace-segregation exposure (2023)",
       subtitle = "Denver MSA tracts, binned means by exclusionary-zoning tercile",
       x = "residential D (RAC, beta=0.5)",
       y = "OD-weighted workplace D exposure") + theme_co
ggsave(file.path(DIR_CO_FIG, "co_fig2_binned_scatter_zoning.png"), p2,
       width = 8, height = 5.5, dpi = 300, bg = "white")

## ---- Fig 3: ZB interaction forest --------------------------------------------
fb <- coefs |>
  filter(grepl("^ZB_mod_whiteblack_", model_id), grepl(":", term))
if (nrow(fb) > 0) {
  p3 <- fb |>
    mutate(mod = sub("^ZB_mod_whiteblack_", "", model_id)) |>
    ggplot(aes(estimate, reorder(mod, estimate))) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = .3) +
    geom_pointrange(aes(xmin = conf.low, xmax = conf.high)) +
    labs(title = "Zoning moderation of the res -> workplace-exposure link",
         subtitle = "Coefficient on res_seg x zoning measure; SD units; jurisdiction-clustered CI",
         x = "interaction coefficient", y = NULL) + theme_co
  ggsave(file.path(DIR_CO_FIG, "co_fig3_interaction_forest.png"), p3,
         width = 7.5, height = 4.5, dpi = 300, bg = "white")
}

## ---- Fig 4: jurisdiction scatter ---------------------------------------------
jd <- xs |>
  group_by(jurisd_main) |>
  summarise(res_d = weighted.mean(d_whiteblack_rac_half, n_commuters,
                                  na.rm = TRUE),
            workers = sum(n_commuters), .groups = "drop") |>
  inner_join(jur, by = c("jurisd_main" = "jurisd")) |>
  filter(workers >= 1000)
p4 <- ggplot(jd, aes(pct_res_low, res_d)) +
  geom_point(aes(size = workers), alpha = .6, color = "#1b6ca8") +
  ggrepel::geom_text_repel(aes(label = jurisd_main), size = 2.6,
                           max.overlaps = 12) +
  geom_smooth(method = "lm", se = TRUE, linewidth = .6, color = "#c8442c") +
  scale_size_area(max_size = 10, guide = "none") +
  labs(title = "Jurisdiction zoning regimes and residential segregation (2023)",
       subtitle = "Denver-MSA jurisdictions; point size = resident commuters",
       x = "% of zoned area Residential_Low",
       y = "worker-weighted mean tract residential D") + theme_co
ggsave(file.path(DIR_CO_FIG, "co_fig4_jurisdictions.png"), p4,
       width = 8.5, height = 6, dpi = 300, bg = "white")

## ---- Fig 5: 2011 -> 2023 trajectories by exclusionary tercile ----------------
tercile_2023 <- xs |> select(tract_id, reslow_tercile)
tr <- dat |>
  filter(denver_msa, in_scope, n_commuters >= P3_MIN_COMMUTERS) |>
  inner_join(tercile_2023, by = "tract_id") |>
  group_by(year, reslow_tercile) |>
  summarise(res  = weighted.mean(d_whiteblack_rac_half, n_commuters,
                                 na.rm = TRUE),
            wexp = weighted.mean(wexp_whiteblack_wac_half, n_commuters,
                                 na.rm = TRUE), .groups = "drop") |>
  pivot_longer(c(res, wexp), names_to = "sphere")
p5 <- ggplot(tr, aes(year, value, color = factor(reslow_tercile),
                     linetype = sphere)) +
  geom_line(linewidth = .7) + geom_point(size = .9) +
  scale_color_manual(values = c("#1b6ca8", "#888888", "#c8442c"),
                     labels = c("least", "middle", "most"),
                     name = "Exclusionary\ntercile") +
  scale_linetype_manual(values = c(res = "solid", wexp = "22"),
                        labels = c(res = "residential D",
                                   wexp = "workplace exposure"),
                        name = NULL) +
  labs(title = "Segregation trajectories by 2023 zoning regime, 2011-2023",
       subtitle = "Commuter-weighted means; pandemic years interpret with care (WFH)",
       x = NULL, y = "White-Black D (beta=0.5)") + theme_co
ggsave(file.path(DIR_CO_FIG, "co_fig5_trends_by_tercile.png"), p5,
       width = 8.5, height = 5.5, dpi = 300, bg = "white")

message("65 complete. Figures in ", DIR_CO_FIG)
