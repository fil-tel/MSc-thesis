library(GenomicRanges)
library(tidyverse)
library(plyranges)
source("utils.R")

if(!dir.exists("figs")) dir.create("figs")

# define tracts

tracts_gr <- GRanges(seqnames = "chr1", ranges = IRanges(start = c(100, 100, 130, 90, 90, 50), end = c(150, 150, 150, 170, 140, 80)))
mcols(tracts_gr)$name <- paste0("IND_", 1:length(tracts_gr))
mcols(tracts_gr)$col <- as.factor(c(1,1,2,3,4,5))

# tick_positions <- sort(unique(c(start(tracts_gr), end(tracts_gr))))

tracts_gr %>% as_tibble %>% 
  ggplot(aes(x = start, xend = end, y = name, yend = name, colour = col)) +
  geom_segment(linewidth = 3) + scale_y_discrete(limits=rev) +
  labs(x = "position [bp]", y = "samples") +
  coord_cartesian(xlim = c(0, 200)) +
  # scale_x_continuous(breaks = tick_positions) +
  # theme(panel.grid = element_blank()) + 
  guides(colour = "none")+
  ggtitle("Example of introgressed tracts on chromosomes") -> p1

p1

# function to get the long df to plot
get_long_df <- function(bin_mat){
  bin_long <- bin_mat %>% as.data.frame() %>% rownames_to_column(var = "name")%>% 
    pivot_longer(cols=-name,names_to = "bin", values_to = "value") %>% mutate(
      fill_group = ifelse(value == 1, bin, "0"),
      name = factor(name, levels = rev(unique(name)))  # To match tracts plot order
    )
  bin_long$bin <- factor(bin_long$bin, levels = colnames(bin_mat))
  bin_long$name <- factor(bin_long$name, levels = rev(unique(tracts_gr$name)))
  bin_long
}


# uniqueness approach -----------------------------------------------------


bin_mat_unique <- t(get_bin_mat_unique(tracts_gr))
# long format
bin_unique_long <- get_long_df(bin_mat = bin_mat_unique)

# function to get trhe plots for unique aproach 
plot_unique <- function(bin_mat_long){
  bin_levels <- unique(bin_mat_long$bin)
  bin_colors <- setNames(RColorBrewer::brewer.pal(length(bin_levels), "Set1"), bin_levels)
  bin_colors_0 <- c("0" = "white", bin_colors) 
  p2 <- ggplot(bin_mat_long, aes(x = bin, y = name, fill = fill_group)) +
    geom_tile(color = "grey80") +
    geom_text(aes(label = value), color = "black", size = 4) +
    scale_fill_manual(values = bin_colors_0) +
    labs(x = "Segment coordinates", y = NULL) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      panel.grid = element_blank()
    ) + ggtitle("Uniqueness approach matrix")
  names(bin_colors) <- unique(tracts_gr$col)
  tracts_gr %>% as_tibble %>% 
    ggplot(aes(x = start, xend = end, y = name, yend = name, colour = col)) +
    geom_segment(linewidth = 3) + scale_y_discrete(limits=rev) +
    labs(x = "Position [bp]", y = "Samples") +
    coord_cartesian(xlim = c(0, 200)) +
    # scale_x_continuous(breaks = tick_positions) +
    # theme(panel.grid = element_blank()) +
    scale_colour_manual(values = bin_colors)+
    guides(colour = "none")+
    ggtitle("Example of introgressed tracts on chromosomes") -> p1
  
  list(p1, p2)  
}

p_unique <- plot_unique(bin_mat_long = bin_unique_long)
p_uni_grid <- cowplot::plot_grid(plotlist = p_unique, labels = "AUTO")
ggsave(plot = p_uni_grid, filename = "figs/uni_grid.png", width = 12, height = 5, units = "in")

# windows approach --------------------------------------------------------


bin_mat_windows <- t(get_bin_mat_windows(tracts_gr, window_size = 25, step_size = 25, len_chr = 200))
bin_windows_long <- get_long_df(bin_mat = bin_mat_windows)

# function to get trhe plots for unique aproach 
plot_windows <- function(bin_mat_long){
  bin_levels <- unique(bin_mat_long$bin)
  bin_colors <- setNames(RColorBrewer::brewer.pal(length(bin_levels), "Set1"), bin_levels)
  bin_colors_0 <- c("0" = "white", bin_colors)
  alpha_v <- c("0" = 1, setNames(rep(0.2, length(bin_levels)), bin_levels))
  p2 <- ggplot(bin_mat_long) +
    geom_tile(color = "grey80", aes(x = bin, y = name, fill = fill_group, alpha = fill_group)) +
    scale_alpha_manual(values = alpha_v)+
    geom_text(aes(x = bin, y = name, label = value), color = "black", size = 4) +
    scale_fill_manual(values = bin_colors_0) +
    labs(x = "Genomic bin", y = NULL) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      panel.grid = element_blank()
    ) + ggtitle("Windows approach matrix")
  regions <- tibble(
    xmin = seq(26, 152, 25),
    xmax = seq(51, 177, 25),
    Bins = factor(unique(bin_mat_long$bin), levels = unique(bin_mat_long$bin))
  )
  tracts_gr %>% as_tibble %>% 
    ggplot() +
    geom_segment(aes(x = start, xend = end, y = name, yend = name),linewidth = 3) + scale_y_discrete(limits=rev) +
    labs(x = "Position [bp]", y = "Samples") +
    coord_cartesian(xlim = c(0, 200)) +
    # scale_x_continuous(breaks = tick_positions) +
    theme(panel.grid = element_blank()) +
    ggtitle("Example of introgressed tracts on chromosomes") + 
    geom_vline(xintercept = seq(1, 201, by = 25), color="red", lty=2, alpha = 0.5)+
    geom_rect(data = regions,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = Bins),
              inherit.aes = FALSE, alpha = 0.2) +
    scale_fill_manual(values = bin_colors)+
    guides(fill = "none")-> p1
  
  list(p1, p2)  
}


p_windows <- plot_windows(bin_windows_long)
p_wind_grid <- cowplot::plot_grid(plotlist = p_windows, labels = "AUTO")
ggsave(plot = p_wind_grid, filename = "figs/wind_grid.png", width = 12, height = 5, units = "in")

# subset approach ---------------------------------------------------------

bin_mat_subset <- t(get_bin_mat_subset(tracts_gr))
bin_subset_long <- get_long_df(bin_mat = bin_mat_subset)

# function to get trhe plots for unique aproach 
plot_subset <- function(bin_mat_long){
  tick_positions <- sort(unique(c(start(tracts_gr), end(tracts_gr))))
  bin_levels <- unique(bin_mat_long$bin)
  bin_colors <- setNames(RColorBrewer::brewer.pal(length(bin_levels), "Set1"), bin_levels)
  bin_colors_0 <- c("0" = "white", bin_colors)
  alpha_v <- c("0" = 1, setNames(rep(0.2, length(bin_levels)), bin_levels)) 
  p2 <- ggplot(bin_mat_long) +
    geom_tile(color = "grey80", aes(x = bin, y = name, fill = fill_group, alpha = fill_group)) +
    scale_alpha_manual(values = alpha_v)+
    geom_text(aes(x = bin, y = name, label = value), color = "black", size = 4) +
    scale_fill_manual(values = bin_colors_0) +
    labs(x = "Genomic bin", y = NULL) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      panel.grid = element_blank()
    )+ggtitle("Sub-tracts approach matrix")
  regions <- tibble(
    xmin = c(50, 90, 100, 130, 140, 150),
    xmax = c(80, 100, 130, 140, 150, 170),
    Bins = factor(unique(bin_mat_long$bin), levels = unique(bin_mat_long$bin))
  )
  tracts_gr %>% as_tibble %>% 
    ggplot() +
    geom_segment(aes(x = start, xend = end, y = name, yend = name),linewidth = 3) + scale_y_discrete(limits=rev) +
    labs(x = "Position [bp]", y = "Samples") +
    coord_cartesian(xlim = c(0, 200)) +
    # scale_x_continuous(breaks = tick_positions) +
    theme(panel.grid = element_blank()) +
    ggtitle("Example of introgressed tracts on chromosomes") + 
    geom_vline(xintercept = tick_positions, color="red", lty=2, alpha = 0.5)+
    geom_rect(data = regions,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = Bins),
              inherit.aes = FALSE, alpha = 0.2) +
    scale_fill_manual(values = bin_colors) +
    guides(fill = "none")-> p1
  
  list(p1, p2)  
}

p_subset <- plot_subset(bin_subset_long)
p_sub_grid <- cowplot::plot_grid(plotlist = p_subset, labels = "AUTO")
ggsave(plot = p_sub_grid, filename = "figs/sub_grid.png", width = 12, height = 5, units = "in")


# recombination sites -----------------------------------------------------

bin_mat_sites <- t(get_bin_mat_sites(tracts_gr))
bin_sites_long <- get_long_df(bin_mat = bin_mat_sites)

plot_sites <- function(bin_mat_long){
  tick_positions <- sort(unique(c(start(tracts_gr), end(tracts_gr))))
  bin_levels <- unique(bin_mat_long$bin)
  bin_colors <- setNames(RColorBrewer::brewer.pal(length(bin_levels), "Set1"), bin_levels)
  bin_colors_0 <- c("0" = "white", bin_colors)
  alpha_v <- c("0" = 1, setNames(rep(0.2, length(bin_levels)), bin_levels)) 
  p2 <- ggplot(bin_mat_long) +
    geom_tile(color = "grey80", aes(x = bin, y = name, fill = fill_group, alpha = fill_group)) +
    scale_alpha_manual(values = alpha_v)+
    geom_text(aes(x = bin, y = name, label = value), color = "black", size = 4) +
    scale_fill_manual(values = bin_colors_0) +
    labs(x = "Recombination breakpoints coordinates", y = NULL) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      panel.grid = element_blank()
    )+ggtitle("Recombination breakpoints approach matrix")
  
  tracts_gr %>% as_tibble %>% 
    ggplot() +
    geom_segment(aes(x = start, xend = end, y = name, yend = name),linewidth = 3) + scale_y_discrete(limits=rev) +
    labs(x = "Position [bp]", y = "Samples") +
    coord_cartesian(xlim = c(0, 200)) +
    # scale_x_continuous(breaks = tick_positions) +
    theme(panel.grid = element_blank()) +
    ggtitle("Example of introgressed tracts on chromosomes") + 
    geom_vline(xintercept = tick_positions, color=bin_colors, lty=2, alpha = 0.7) +
    guides(fill = "none")-> p1
  
  list(p1, p2)  
}

p_sites <- plot_sites(bin_mat_long = bin_sites_long)
p_site_grid <- cowplot::plot_grid(plotlist = p_sites, labels = "AUTO")
ggsave(plot = p_site_grid, filename = "figs/site_grid.png", width = 12, height = 5, units = "in")


# ALL together

grid_all <- cowplot::plot_grid(plotlist = c(p_unique, p_windows, p_subset, p_sites), ncol = 2, labels = "AUTO")
ggsave(plot = grid_all, filename = "figs/grid_all.png", width = 12, height = 20, units = "in")


# table for latex

# df <- df %>% rename("chrom"=seqnames) %>% select(-c(strand, width))
# print(xtable(df, caption = "Test", label = "method:tab_tracts"), include.rownames=FALSE)
