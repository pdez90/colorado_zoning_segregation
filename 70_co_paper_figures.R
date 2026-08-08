# ==============================================================================
# 70_co_paper_figures.R      [PAPER 4, step 10]
# Publication figures for the JAPA manuscript.
#
#   Fig 1  A. generalized zoning districts (the raw layer, dissolved by class)
#          B. tract-level exclusionary measure
#   Fig 2  three spheres side by side: exclusionary zoning | residential
#          segregation | OD-weighted workplace-segregation exposure
#   Fig 3  THE ARGUMENT: moderation coefficient by rung (naive / +metropolitan
#          position / +local job surface) for each zoning measure
#   Fig 4  A. binned residential-vs-workplace coupling by exclusionary tercile
#          B. native-unit gradients (commute km, jobs-housing ratio)
#   SI S1  zoning coverage screen; SI S2 ADU jurisdictional patchwork;
#   SI S3  accessibility surface (the confound, mapped)
#
# Palette: single-hue sequential (light -> dark) for all magnitude maps so they
# survive grayscale printing; a colorblind-checked categorical set for the
# land-use classes; low-density residential deliberately given a recessive
# neutral so the job-permitting land reads against it.
#
# Needs: 61 (geometry+zoning), 63 (panel), 68 (accessibility + ladder), 69.
# Output: output/figures/p4_fig1..fig4.png, p4_figS1..S3.png
# ==============================================================================

source("60_co_setup.R")
library(patchwork)

geom  <- readRDS(file.path(DIR_CO_CLEAN, "co_tract_geom.rds"))
zoning<- readRDS(file.path(DIR_CO_CLEAN, "co_tract_zoning.rds"))
dat   <- readRDS(file.path(DIR_CO_OUT, "co_analysis_home_panel.rds"))
acc   <- readRDS(file.path(DIR_CO_CLEAN, "co_accessibility_2023.rds"))
ladder<- read.csv(file.path(DIR_CO_MOD, "p4_accessibility_ladder.csv"))

xs <- dat |>
  filter(year == CO_ANCHOR_YEAR, in_scope, denver_msa,
         n_commuters >= P3_MIN_COMMUTERS) |>
  left_join(acc, by = "tract_id")
g <- geom |> inner_join(xs, by = "tract_id")

# zoom to the urbanised core: the rural rectangles otherwise dominate the frame
ctr <- suppressWarnings(st_coordinates(st_centroid(g)))
xlim <- quantile(ctr[, 1], c(.01, .93)) + c(-9000, 9000)
ylim <- quantile(ctr[, 2], c(.02, .98)) + c(-9000, 9000)

SEQ <- "Blues"
theme_map <- theme_void(base_size = 9) +
  theme(legend.position = "bottom", legend.key.height = unit(3, "mm"),
        legend.key.width = unit(7, "mm"),
        plot.title = element_text(size = 10, hjust = 0))
theme_p4 <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10.5))

qmap <- function(df, var, title, legend_title) {
  df$.q <- cut(df[[var]], breaks = quantile(df[[var]], seq(0, 1, .2),
                                            na.rm = TRUE),
               include.lowest = TRUE, dig.lab = 2)
  ggplot(df) +
    geom_sf(aes(fill = .q), color = "white", linewidth = .05) +
    scale_fill_brewer(palette = SEQ, name = legend_title, na.value = "grey85") +
    coord_sf(xlim = xlim, ylim = ylim) +
    labs(title = title) + theme_map
}

## ---- Fig 1: the zoning data --------------------------------------------------
zon <- st_read(CO_ZONING_SHP, quiet = TRUE) |> st_zm(drop = TRUE) |>
  st_transform(CRS_METERS) |> st_make_valid()
zon <- zon[!st_is_empty(zon), ] |> mutate(GenZone2 = trimws(GenZone2))
cls_map <- c(Residential_Low = "Residential: low density",
  Residential_Med = "Residential: medium/high",
  Residential_MedHigh = "Residential: medium/high",
  Residential_High = "Residential: medium/high",
  MobileHome = "Mobile home",
  MixedUseRes_Low = "Mixed use w/ residential",
  MixedUseRes_Med = "Mixed use w/ residential",
  MixedUseRes_MedHigh = "Mixed use w/ residential",
  MixedUseRes_High = "Mixed use w/ residential",
  MixedUse_Conditional = "Mixed use w/ residential",
  Commercial = "Commercial", Industrial = "Industrial",
  OpenSpace = "Open space / civic", Civic = "Open space / civic",
  Uncertain = "Uncertain")
zon$cls <- factor(cls_map[zon$GenZone2],
  levels = c("Residential: low density", "Residential: medium/high",
             "Mixed use w/ residential", "Commercial", "Industrial",
             "Open space / civic", "Mobile home", "Uncertain"))
zd <- zon |> group_by(cls) |> summarise(.groups = "drop")
LU_COLS <- c("#e8dcc8", "#2a78d6", "#1baf7a", "#eb6834", "#4a3aa7",
             "#008300", "#e87ba4", "#b9b8b0")

p1a <- ggplot(zd) +
  geom_sf(aes(fill = cls), color = NA) +
  scale_fill_manual(values = LU_COLS, name = NULL, drop = FALSE) +
  coord_sf(xlim = xlim, ylim = ylim) +
  guides(fill = guide_legend(ncol = 2, override.aes = list(color = "grey40"))) +
  labs(title = sprintf("A. Generalized zoning districts (Oct %d)",
                       CO_ANCHOR_YEAR)) + theme_map
p1b <- qmap(g, "pct_reslow_of_res", "B. Exclusionary zoning, by census tract",
            "% of residential land zoned low-density (quintiles)")
ggsave(file.path(DIR_CO_FIG, "p4_fig1_zoning.png"), p1a + p1b,
       width = 11.5, height = 6.6, dpi = 400, bg = "white")

## ---- Fig 2: the three spheres ------------------------------------------------
p2 <- qmap(g, "pct_reslow_of_res", "A. Exclusionary zoning", "quintiles") +
  qmap(g, "d_whiteblack_rac_half", "B. Residential segregation", "quintiles") +
  qmap(g, "wexp_whiteblack_wac_half",
       "C. Workplace-segregation exposure", "quintiles")
ggsave(file.path(DIR_CO_FIG, "p4_fig2_three_spheres.png"), p2,
       width = 14.5, height = 5.6, dpi = 400, bg = "white")

## ---- Fig 3: THE ARGUMENT -----------------------------------------------------
LAB <- c(pct_reslow_of_res = "Exclusionary\n(% res. land low-density)",
         pct_res_low = "Exclusionary\n(% all zoned land)",
         zoning_entropy = "Use mix (entropy)",
         pct_job_zone = "Job-permitting land",
         pct_adu_res = "ADU-permitting\n(see measurement caveat)")
SPECLAB <- c(A_total = "Naive (no position control)",
             B_plus_regional_position = "+ metropolitan position (preferred)",
             C_plus_local_jobsurface = "+ local job surface (over-control)")
f3 <- ladder |>
  filter(grepl("^z_d_whiteblack_rac_half:z_", term),
         grepl("__(A_total|B_plus_regional_position|C_plus_local_jobsurface)$",
               model_id)) |>
  mutate(zoning = sub("__.*", "", model_id), spec = sub(".*__", "", model_id)) |>
  # keep ONLY each model's own zoning interaction (rung C also carries the
  # job-surface and distance interactions in the same table)
  filter(zoning %in% names(LAB),
         term == paste0("z_d_whiteblack_rac_half:z_", zoning)) |>
  mutate(zlab = factor(LAB[zoning], levels = rev(LAB)),
         slab = factor(SPECLAB[spec], levels = SPECLAB))
p3 <- ggplot(f3, aes(estimate, zlab, color = slab)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = .3, color = "grey50") +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high), size = .4,
                  position = position_dodge(width = .55)) +
  scale_color_manual(values = c("#2a78d6", "#eb6834", "#b9b8b0"), name = NULL) +
  labs(title = "Estimated zoning moderation shrinks by 60-72% once metropolitan position is included",
       subtitle = "Coefficient on residential segregation x zoning measure; SD units; jurisdiction-clustered 95% CI",
       x = "interaction coefficient", y = NULL) +
  theme_p4 + theme(legend.position = "bottom")
ggsave(file.path(DIR_CO_FIG, "p4_fig3_ladder.png"), p3,
       width = 9, height = 5.4, dpi = 400, bg = "white")

## ---- Fig 4: coupling + native-unit gradients ---------------------------------
xt <- xs |> mutate(terc = factor(ntile(pct_reslow_of_res, 3), 1:3,
                                 c("least", "middle", "most")))
p4a <- xt |>
  mutate(bin = ntile(d_whiteblack_rac_half, 20)) |>
  group_by(terc, bin) |>
  summarise(res = mean(d_whiteblack_rac_half),
            wexp = mean(wexp_whiteblack_wac_half), .groups = "drop") |>
  ggplot(aes(res, wexp, color = terc)) +
  geom_point(size = 1.3) + geom_smooth(method = "lm", se = FALSE, linewidth = .6) +
  scale_color_manual(values = c("#c6dbef", "#6baed6", "#08519c"),
                     name = "Exclusionary tercile") +
  labs(title = "A. Coupling of the two spheres (unadjusted)",
       x = "residential dissimilarity", y = "workplace-segregation exposure") +
  theme_p4 + theme(legend.position = "bottom")
p4b <- xt |>
  select(terc, `Mean commute (km)` = mean_dist_km,
         `Jobs per resident worker` = jobs_housing_ratio,
         `Distance to CBD (km)` = dist_cbd_km) |>
  pivot_longer(-terc) |>
  group_by(terc, name) |>
  summarise(m = mean(value, na.rm = TRUE), .groups = "drop") |>
  ggplot(aes(m, terc, fill = terc)) +
  geom_col(width = .62) +
  geom_text(aes(label = sprintf("%.1f", m)), hjust = -0.15, size = 3,
            color = "#0b0b0b") +
  facet_wrap(~name, scales = "free_x") +
  scale_fill_manual(values = c("#c6dbef", "#6baed6", "#08519c"), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, .18))) +
  labs(title = "B. What the exclusionary tercile looks like on the ground",
       x = NULL, y = NULL) + theme_p4
ggsave(file.path(DIR_CO_FIG, "p4_fig4_coupling_gradients.png"),
       p4a / p4b + plot_layout(heights = c(1.25, 1)),
       width = 9, height = 8.4, dpi = 400, bg = "white")

## ---- SI figures --------------------------------------------------------------
cov_all <- geom |> inner_join(zoning, by = "tract_id")
pS1 <- ggplot(cov_all) +
  geom_sf(aes(fill = pmin(cover_zoned, 1)), color = "white", linewidth = .05) +
  scale_fill_distiller(palette = SEQ, direction = 1,
                       name = "zoned share of tract area") +
  coord_sf(xlim = xlim, ylim = ylim) +
  labs(title = sprintf("Figure S1. Zoning coverage; analysis keeps tracts >= %.0f%%",
                       100 * CO_MIN_ZONED_COVER)) + theme_map
ggsave(file.path(DIR_CO_FIG, "p4_figS1_coverage.png"), pS1,
       width = 7, height = 6.4, dpi = 400, bg = "white")

pS2 <- qmap(g, "pct_adu_res",
            "Figure S2. ADU permissiveness is a jurisdictional patchwork",
            "% of residential land allowing ADUs (quintiles)")
ggsave(file.path(DIR_CO_FIG, "p4_figS2_adu.png"), pS2,
       width = 7, height = 6.4, dpi = 400, bg = "white")

pS3 <- qmap(g, "log_jobs_grav",
            "Figure S3. The confound: gravity-weighted job accessibility",
            "log jobs accessible (quintiles)")
ggsave(file.path(DIR_CO_FIG, "p4_figS3_accessibility.png"), pS3,
       width = 7, height = 6.4, dpi = 400, bg = "white")

message("70 complete. Figures in ", DIR_CO_FIG)
