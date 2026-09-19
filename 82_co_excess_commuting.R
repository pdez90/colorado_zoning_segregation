# ==============================================================================
# 82_co_excess_commuting.R      [PAPER 4, step 22 -- SI extension: excess commuting]
# Do exclusionary areas commute farther than the metropolitan jobs-housing
# geography requires? Classic excess-commuting decomposition (Hamilton 1982;
# White 1988): solve the transportation problem that assigns resident workers
# (RAC) to jobs (WAC) minimizing aggregate commute distance, and compare the
# minimum-required commute with actual OD commuting.
#
#   excess % = (T_actual - T_minimum) / T_actual
#
# Two readings the paper's 13.5 -> 19.2 km gradient cannot distinguish on its
# own, and this can:
#   (a) exclusionary neighborhoods are simply FAR from jobs -> their minimum
#       commute is also long; excess similar to everywhere else;
#   (b) their residents travel farther than job geography requires -> high
#       excess concentrated in exclusionary areas.
#
# Implementation notes (state in SI):
#   - Distances are centroid Euclidean km, IDENTICAL metric for actual and
#     minimum (within-tract = 0 for both), so the comparison is internally
#     consistent even though both understate network distance.
#   - Jobs exceed same-frame resident workers, so job masses are scaled
#     proportionally to total workers (standard balancing; stated).
#   - Solved for all workers and separately for low- (CE01/SE01) and
#     high-earnings (CE03/SE03) workers -- matching workers to jobs of their
#     own earnings band.
#   - Per-origin decomposition: the optimal assignment yields a minimum mean
#     commute for each origin tract; tercile summaries compare actual vs
#     minimum vs excess. NOTE the minimum is a system optimum -- per-origin
#     values depend on the whole configuration -- so tercile rows are
#     descriptive of the optimum's geography, not tract-level entitlements.
#   - Requires the 'transport' package (fast network-flow OT solver):
#       install.packages("transport")
#
# Output: output/models/p4_excess_commuting.csv          (overall + by group)
#         output/models/p4_excess_commuting_tercile.csv
#         output/figures/p4_fig_excess_commuting.png
# ==============================================================================

source("60_co_setup.R")
if (!requireNamespace("transport", quietly = TRUE))
  stop("Package 'transport' is required: install.packages(\"transport\")")

## ---- frames ------------------------------------------------------------------
dat <- readRDS(file.path(DIR_CO_OUT, "co_analysis_home_panel.rds"))
xs <- dat |>
  filter(year == CO_ANCHOR_YEAR, in_scope, denver_msa,
         n_commuters >= P3_MIN_COMMUTERS) |>
  transmute(tract_id, tercile = ntile(pct_reslow_of_res, 3))
cent <- readRDS(file.path(DIR_CLEAN, "tract_centroids_km.rds")) |>
  transmute(tract_id = as.character(GEOID), X_km, Y_km)
rac <- readRDS(file.path(DIR_CO_CLEAN,
                         sprintf("co_rac_tract_%d.rds", CO_ANCHOR_YEAR)))
wac <- readRDS(file.path(DIR_CO_CLEAN,
                         sprintf("co_wac_tract_%d.rds", CO_ANCHOR_YEAR)))
od <- map(P3_OD_PARTS, function(part) {
  for (f in c(p3_od_cache(part, CO_ANCHOR_YEAR, "co"),
              co_od_cache(part, CO_ANCHOR_YEAR)))
    if (file.exists(f)) return(readRDS(f))
  NULL
}) |> compact() |> bind_rows() |>
  group_by(w_tract, h_tract) |>
  summarise(across(any_of(c("S000", "SE01", "SE03")), ~ sum(.x, na.rm = TRUE)),
            .groups = "drop")

fr <- xs |> inner_join(cent, by = "tract_id")
D <- as.matrix(dist(fr[, c("X_km", "Y_km")]))   # symmetric, 0 diagonal
idx <- setNames(seq_len(nrow(fr)), fr$tract_id)

## ---- solver ------------------------------------------------------------------
# Supply and demand are the ORIGIN and DESTINATION MARGINS of the observed flow
# matrix between frame tracts, for the same earnings band (White 1988). Actual
# and minimum commuting are therefore computed for exactly the same workers and
# jobs, totals balance by construction, and no rescaling is needed.
solve_group <- function(od_col, label) {
  a <- od |>
    filter(h_tract %in% fr$tract_id, w_tract %in% fr$tract_id,
           .data[[od_col]] > 0) |>
    mutate(flow = .data[[od_col]],
           dkm = D[cbind(idx[h_tract], idx[w_tract])])
  s <- numeric(nrow(fr)); d <- numeric(nrow(fr))
  so <- tapply(a$flow, a$h_tract, sum); s[idx[names(so)]] <- so
  de <- tapply(a$flow, a$w_tract, sum); d[idx[names(de)]] <- de
  stopifnot(isTRUE(all.equal(sum(s), sum(d))))
  plan <- transport::transport(s, d, costm = D, method = "networkflow")
  # plan: from (origin idx), to (dest idx), mass
  t_min_total <- sum(plan$mass * D[cbind(plan$from, plan$to)])
  min_by_o <- tapply(plan$mass * D[cbind(plan$from, plan$to)], plan$from, sum) /
              tapply(plan$mass, plan$from, sum)
  o_min <- tibble(tract_id = fr$tract_id[as.integer(names(min_by_o))],
                  min_km = as.numeric(min_by_o))
  t_act <- weighted.mean(a$dkm, a$flow)
  o_act <- a |> group_by(tract_id = h_tract) |>
    summarise(act_km = weighted.mean(dkm, flow),
              n = sum(flow), .groups = "drop")
  t_min <- t_min_total / sum(s)
  list(summary = tibble(group = label,
                        workers = sum(s),
                        actual_km = t_act, minimum_km = t_min,
                        excess_km = t_act - t_min,
                        excess_pct = 100 * (t_act - t_min) / t_act),
       origins = o_act |> left_join(o_min, by = "tract_id") |>
         mutate(group = label))
}

res_all  <- solve_group("S000", "all")
res_low  <- solve_group("SE01", "low_earnings")
res_high <- solve_group("SE03", "high_earnings")
summ <- bind_rows(res_all$summary, res_low$summary, res_high$summary)
write.csv(summ, file.path(DIR_CO_MOD, "p4_excess_commuting.csv"),
          row.names = FALSE)
print(as.data.frame(summ), digits = 4)

## ---- tercile decomposition (all workers) -------------------------------------
terc <- res_all$origins |>
  inner_join(xs, by = "tract_id") |>
  group_by(tercile) |>
  summarise(n_tracts = n(),
            actual_km  = weighted.mean(act_km, n, na.rm = TRUE),
            minimum_km = weighted.mean(min_km, n, na.rm = TRUE),
            excess_km  = actual_km - minimum_km,
            excess_pct = 100 * excess_km / actual_km, .groups = "drop")
write.csv(terc, file.path(DIR_CO_MOD, "p4_excess_commuting_tercile.csv"),
          row.names = FALSE)
print(as.data.frame(terc), digits = 3)

pE <- terc |>
  select(tercile, actual = actual_km, `minimum required` = minimum_km) |>
  pivot_longer(-tercile) |>
  ggplot(aes(factor(tercile, 1:3, c("least", "middle", "most")), value,
             fill = name)) +
  geom_col(position = position_dodge(width = .7), width = .62) +
  geom_text(aes(label = sprintf("%.1f", value)),
            position = position_dodge(width = .7), vjust = -0.4, size = 2.9) +
  scale_fill_manual(values = c(actual = "#08519c",
                               `minimum required` = "#c6dbef"), name = NULL) +
  labs(title = "Actual vs minimum-required commuting, by exclusionary tercile",
       subtitle = paste("Minimum from the optimal assignment of resident",
                        "workers to jobs (transportation problem);",
                        "\nthe gap is excess commuting. Centroid km,",
                        "identical metric for both."),
       x = "exclusionary-zoning tercile of home tract",
       y = "mean commute (km)") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")
ggsave(file.path(DIR_CO_FIG, "p4_fig_excess_commuting.png"), pE,
       width = 6.8, height = 5.2, dpi = 350, bg = "white")
message("82 complete.")
