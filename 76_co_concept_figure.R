# ==============================================================================
# 76_co_concept_figure.R      [PAPER 4, step 16]
# Figure 1 for the reorganized manuscript: the paper's argument as a picture.
#   Panel A -- the conceptual chain the paper studies:
#       residential neighborhood
#         -> accessible labor market  ("the menu" of reachable jobs;
#            shaped by zoning, metropolitan position, transportation network)
#         -> actual commuting destinations
#         -> workplace-location segregation
#     with the co-evolution caution printed under the diagram: these are
#     descriptive relationships, not a causal path model.
#   Panel B -- the menu result, simplified from Figure S9's tract scatter:
#     mean segregation of the ACCESSIBLE labor market vs mean REALIZED
#     workplace exposure by decile of residential segregation, from
#     output/models/p4_opportunity_vs_realized.csv (written by 75). The two
#     lines rise together and the shaded sorting band between them stays thin:
#     realized exposure tracks the composition of what is reachable (92% of
#     variance), not sorting within it (8%).
#
# Requires: 75 has run (the decile CSV exists).
# Output:   output/figures/p4_fig1_concept.png   (main-text Figure 1)
# ==============================================================================

source("60_co_setup.R")
library(patchwork)

opp_file <- file.path(DIR_CO_MOD, "p4_opportunity_vs_realized.csv")
stopifnot(file.exists(opp_file))

## ---- Panel A: the conceptual chain -------------------------------------------
box <- function(p, x, y, w, h, label, size = 3.1)
  p + annotate("rect", xmin = x - w/2, xmax = x + w/2, ymin = y - h/2,
               ymax = y + h/2, fill = "#eaf1f8", colour = "#08519c",
               linewidth = .45) +
      annotate("text", x = x, y = y, label = label, size = size, lineheight = .95)
arrow_v <- function(p, x, y0, y1)
  p + annotate("segment", x = x, xend = x, y = y0, yend = y1,
               linewidth = .5, colour = "#08519c",
               arrow = arrow(length = unit(5, "pt"), type = "closed"))

pA <- ggplot() + xlim(0, 10) + ylim(-0.6, 10) + theme_void()
pA <- box(pA, 5.6, 9.1, 6.4, 1.15, "Residential neighborhood")
pA <- box(pA, 5.6, 6.0, 6.4, 1.5,
          "Accessible labor market\n(“the menu” of reachable jobs)")
pA <- box(pA, 5.6, 3.2, 6.4, 1.15, "Actual commuting destinations")
pA <- box(pA, 5.6, 0.9, 6.4, 1.15, "Workplace-location segregation")
pA <- arrow_v(pA, 5.6, 8.50, 6.85)
pA <- arrow_v(pA, 5.6, 5.22, 3.85)
pA <- arrow_v(pA, 5.6, 2.60, 1.55)
# what shapes the menu: three side inputs converging on the second box
pA <- pA +
  annotate("text", x = 1.05, y = c(8.6, 7.6, 6.5), hjust = 0.5, size = 2.7,
           fontface = "italic", colour = "grey25",
           label = c("Zoning", "Metropolitan\nposition", "Transportation\nnetwork")) +
  annotate("segment",
           x = c(1.75, 1.85, 1.85), y = c(8.35, 7.45, 6.55),
           xend = 2.32, yend = 6.55, colour = "grey55", linewidth = .4,
           arrow = arrow(length = unit(3.5, "pt"), type = "closed")) +
  annotate("text", x = 5.2, y = -0.45, size = 2.5, fontface = "italic",
           colour = "grey40", lineheight = .95,
           label = paste("Descriptive relationships: zoning, infrastructure,",
                         "employment location,\nand residential sorting co-evolved."))

## ---- Panel B: the menu result, by decile -------------------------------------
opp <- read.csv(opp_file) |>
  filter(split == "res_seg_decile", !is.na(accessible)) |>
  mutate(group = as.integer(group))
pB <- ggplot(opp, aes(group)) +
  geom_ribbon(aes(ymin = accessible, ymax = realized), fill = "#c6dbef",
              alpha = .7) +
  geom_line(aes(y = accessible, colour = "menu"), linetype = 2,
            linewidth = .6) +
  geom_point(aes(y = accessible, colour = "menu"), shape = 15, size = 1.9) +
  geom_line(aes(y = realized, colour = "realized"), linewidth = .7) +
  geom_point(aes(y = realized, colour = "realized"), size = 2) +
  scale_colour_manual(
    values = c(realized = "#08519c", menu = "#6baed6"),
    labels = c(realized = "realized workplace exposure",
               menu = "segregation of the accessible labor market"),
    breaks = c("realized", "menu"), name = NULL) +
  scale_x_continuous(breaks = 1:10) +
  labs(x = "decile of residential segregation",
       y = "workplace-location segregation",
       subtitle = paste("Shaded band: sorting within the menu",
                        "(8% of variance in realized exposure)")) +
  theme_minimal(base_size = 9.5) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        legend.direction = "vertical", legend.margin = margin(t = -4),
        plot.subtitle = element_text(size = 7.5, colour = "grey30"))

fig <- pA + pB + plot_layout(widths = c(1, 1.15)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold"))
ggsave(file.path(DIR_CO_FIG, "p4_fig1_concept.png"), fig,
       width = 10.2, height = 4.7, dpi = 350, bg = "white")
message("76 complete: p4_fig1_concept.png written.")
