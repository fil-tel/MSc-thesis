library(ggplotify)
library(dplyr)
library(ggplot2)
library(cowplot)

l1p_plots <- readRDS("plots_1p.rds")
l2p_plots <- readRDS("plots_2p.rds")

for (name in names(l1p_plots)) {
  p1 <- l1p_plots[[name]]
  p2 <- l2p_plots[[name]]
  grid <- plot_grid(p1, p2, labels = "AUTO", ncol = 1)
  ggsave(filename = paste0("figs/" , name, ".png"), plot = grid, width = 10, height = 13)
}
