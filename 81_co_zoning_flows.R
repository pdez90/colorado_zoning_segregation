# ==============================================================================
# 81_co_zoning_flows.R      [PAPER 4, step 21 -- SI extension: zoning-to-zoning]
# How does local zoning create metropolitan interdependence between where
# workers are housed and where jobs are allowed? Every OD flow has a zoning
# regime at its residential origin and another at its workplace destination;
# this script treats commuting as the connector between land-use regimes.
#
#  PART A -- the zoning-to-zoning flow matrix. Origins classified by
#    exclusionary tercile (pct_reslow_of_res), destinations by tercile of
#    employment-permitting land (pct_job_zone). For each 3x3 cell: share of
#    commuters, ratio of observed flow to the flow expected under
#    independence of origin and destination classes (O/E; >1 = the two
#    regimes are commuting complements), and mean flow distance.
#
#  PART B -- tract-pair gravity model (PPML / fepois on all ordered pairs of
#    frame tracts, zeros included). Flow_ij ~ log origin workers + log
#    destination jobs + log distance + same-jurisdiction indicator
#    + origin zoning (exclusionary share, use-mix entropy)
#    + destination zoning (employment-permitting share)
#    + the interdependence terms:
#        excl_origin x jobzone_destination   (complementary-regime pairing)
#        entropy_origin x log distance       (do mixed-use origins generate
#                                             more distance-contained flows?)
#    Errors clustered on origin jurisdiction. Descriptive gravity structure,
#    not causal identification -- zoning at both ends co-evolved with the
#    flow field.
#
# Output: output/models/p4_zoning_flow_matrix.csv
#         output/models/p4_gravity_models.csv
#         output/figures/p4_fig_zoning_flows.png
# ==============================================================================

source("60_co_setup.R")
library(fixest)

## ---- frames ------------------------------------------------------------------
dat <- readRDS(file.path(DIR_CO_OUT, "co_analysis_home_panel.rds"))
xs <- dat |>
  filter(year == CO_ANCHOR_YEAR, in_scope, denver_msa,
         n_commuters >= P3_MIN_COMMUTERS) |>
  transmute(tract_id, jurisd_main,
            pct_reslow_of_res, pct_job_zone, zoning_entropy,
            terc_excl = ntile(pct_reslow_of_res, 3))
xs$terc_jobz <- ntile(xs$pct_job_zone, 3)
cent <- readRDS(file.path(DIR_CLEAN, "tract_centroids_km.rds")) |>
  transmute(tract_id = as.character(GEOID), X_km, Y_km)
rac <- readRDS(file.path(DIR_CO_CLEAN,
                         sprintf("co_rac_tract_%d.rds", CO_ANCHOR_YEAR))) |>
  transmute(tract_id, workers = C000)
wac <- readRDS(file.path(DIR_CO_CLEAN,
                         sprintf("co_wac_tract_%d.rds", CO_ANCHOR_YEAR))) |>
  transmute(tract_id, jobs = C000)
od <- map(P3_OD_PARTS, function(part) {
  for (f in c(p3_od_cache(part, CO_ANCHOR_YEAR, "co"),
              co_od_cache(part, CO_ANCHOR_YEAR)))
    if (file.exists(f)) return(readRDS(f))
  NULL
}) |> compact() |> bind_rows() |>
  group_by(w_tract, h_tract) |>
  summarise(S000 = sum(S000, na.rm = TRUE), .groups = "drop")

## ---- PART A: the 3x3 regime matrix -------------------------------------------
odm <- od |>
  inner_join(xs |> select(tract_id, terc_excl), by = c(h_tract = "tract_id")) |>
  inner_join(xs |> select(tract_id, terc_jobz), by = c(w_tract = "tract_id")) |>
  left_join(cent |> rename(xh = X_km, yh = Y_km), by = c(h_tract = "tract_id")) |>
  left_join(cent |> rename(xw = X_km, yw = Y_km), by = c(w_tract = "tract_id")) |>
  mutate(dist_km = sqrt((xh - xw)^2 + (yh - yw)^2))
tot <- sum(odm$S000)
mat <- odm |>
  group_by(terc_excl, terc_jobz) |>
  summarise(flow = sum(S000),
            mean_dist_km = weighted.mean(dist_km, S000), .groups = "drop") |>
  group_by(terc_excl) |> mutate(orig_share = sum(flow) / tot) |> ungroup() |>
  group_by(terc_jobz) |> mutate(dest_share = sum(flow) / tot) |> ungroup() |>
  mutate(share = flow / tot,
         expected = orig_share * dest_share,
         obs_over_exp = share / expected) |>
  select(terc_excl, terc_jobz, flow, share, obs_over_exp, mean_dist_km)
write.csv(mat, file.path(DIR_CO_MOD, "p4_zoning_flow_matrix.csv"),
          row.names = FALSE)
print(as.data.frame(mat), digits = 3)

pM <- ggplot(mat, aes(factor(terc_jobz), factor(terc_excl))) +
  geom_tile(aes(fill = obs_over_exp), color = "white", linewidth = .6) +
  geom_text(aes(label = sprintf("O/E %.2f\n%.1f km", obs_over_exp,
                                mean_dist_km)), size = 2.9) +
  scale_fill_gradient2(low = "#c6dbef", mid = "white", high = "#08519c",
                       midpoint = 1, name = "observed / expected flow") +
  scale_x_discrete(labels = c("least", "middle", "most")) +
  scale_y_discrete(labels = c("least", "middle", "most")) +
  labs(title = "Which land-use regimes does commuting connect?",
       subtitle = paste("Commuter flows by exclusionary tercile of the origin",
                        "and employment-zoning tercile of the destination;",
                        "\nO/E > 1: the regimes are commuting complements.",
                        "Cell text: O/E ratio and mean flow distance."),
       x = "destination: employment-permitting land (tercile)",
       y = "origin: exclusionary residential zoning (tercile)") +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(), legend.position = "bottom")
ggsave(file.path(DIR_CO_FIG, "p4_fig_zoning_flows.png"), pM,
       width = 6.8, height = 6.2, dpi = 350, bg = "white")

## ---- PART B: tract-pair gravity (PPML) ---------------------------------------
zsc <- function(x) as.numeric(scale(x))
fr <- xs |>
  left_join(rac, by = "tract_id") |>
  left_join(wac, by = "tract_id") |>
  left_join(cent, by = "tract_id") |>
  filter(!is.na(X_km))
pairs <- expand_grid(h_tract = fr$tract_id, w_tract = fr$tract_id) |>
  filter(h_tract != w_tract) |>
  left_join(fr |> select(h_tract = tract_id, o_workers = workers,
                         o_excl = pct_reslow_of_res, o_ent = zoning_entropy,
                         o_jur = jurisd_main, xh = X_km, yh = Y_km),
            by = "h_tract") |>
  left_join(fr |> select(w_tract = tract_id, d_jobs = jobs,
                         d_jobz = pct_job_zone, d_jur = jurisd_main,
                         xw = X_km, yw = Y_km),
            by = "w_tract") |>
  left_join(od, by = c("w_tract", "h_tract")) |>
  mutate(S000 = replace_na(S000, 0),
         dist_km = sqrt((xh - xw)^2 + (yh - yw)^2),
         log_dist = log(dist_km),
         same_jur = as.integer(o_jur == d_jur),
         z_o_excl = zsc(o_excl), z_o_ent = zsc(o_ent),
         z_d_jobz = zsc(d_jobz)) |>
  filter(is.finite(log_dist), o_workers > 0, d_jobs > 0)
message(sprintf("Gravity frame: %s ordered tract pairs, %.1f%% with flow > 0",
                format(nrow(pairs), big.mark = ","),
                100 * mean(pairs$S000 > 0)))

tidy1 <- function(fit, id, keep) {
  ct <- as.data.frame(summary(fit)$coeftable)
  ct$term <- rownames(ct); rownames(ct) <- NULL
  names(ct)[1:4] <- c("estimate", "std.error", "statistic", "p.value")
  ct |> filter(term %in% keep) |> mutate(model_id = id, n_obs = fit$nobs)
}
K1 <- c("log(o_workers)", "log(d_jobs)", "log_dist", "same_jur")
K2 <- c(K1, "z_o_excl", "z_o_ent", "z_d_jobz")
K3 <- c(K2, "z_o_excl:z_d_jobz", "z_o_ent:log_dist")
g1 <- fepois(S000 ~ log(o_workers) + log(d_jobs) + log_dist + same_jur,
             data = pairs, cluster = ~o_jur)
g2 <- fepois(S000 ~ log(o_workers) + log(d_jobs) + log_dist + same_jur +
               z_o_excl + z_o_ent + z_d_jobz,
             data = pairs, cluster = ~o_jur)
g3 <- fepois(S000 ~ log(o_workers) + log(d_jobs) + log_dist + same_jur +
               z_o_excl + z_o_ent + z_d_jobz +
               z_o_excl:z_d_jobz + z_o_ent:log_dist,
             data = pairs, cluster = ~o_jur)
grav <- bind_rows(tidy1(g1, "G1_baseline", K1),
                  tidy1(g2, "G2_zoning", K2),
                  tidy1(g3, "G3_interdependence", K3))
write.csv(grav, file.path(DIR_CO_MOD, "p4_gravity_models.csv"),
          row.names = FALSE)
print(as.data.frame(grav), digits = 3)
message("81 complete.")
