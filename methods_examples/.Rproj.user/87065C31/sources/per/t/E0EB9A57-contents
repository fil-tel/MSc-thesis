library(GenomicRanges)
library(tidyverse)
library(plyranges)
source("utils.R")

if(!dir.exists("figs")) dir.create("figs")

# define tracts

tracts_gr <- GRanges(seqnames = "chr1", ranges = IRanges(start = c(100, 100, 130, 90, 90, 50), end = c(150, 150, 150, 170, 140, 80)))
mcols(tracts_gr)$name <- paste0("IND_", 1:length(tracts_gr))
mcols(tracts_gr)$col <- as.factor(c(1,1,2,3,4,5))


# unique
bin_mat_unique <- get_bin_mat_unique(tracts_gr)
tfs_uni <- compute_tfs(bin_mat_unique)
p_unique <- plot_tfs(tfs_uni, "Uniqueness approach")

# windows
bin_mat_windows <- get_bin_mat_windows(tracts_gr, window_size = 25, step_size = 25, len_chr = 200)
tfs_wind <- compute_tfs(bin_mat_windows)
p_wind <- plot_tfs(tfs_wind, "Windows approach")

# subset
bin_mat_subset <- get_bin_mat_subset(tracts_gr)
tfs_sub <- compute_tfs(bin_mat_subset)
p_sub <- plot_tfs(tfs_sub, "Sub-tracts approach")

# rec sites
bin_mat_sites <- get_bin_mat_sites(tracts_gr)
tfs_sit <- compute_tfs(bin_mat_sites)
p_sit <- plot_tfs(tfs_sit, "Recombination breakpoints approach")


cowplot::plot_grid(p_unique, p_wind, p_sub, p_sit, labels = "AUTO")
ggsave("figs/tfs.png", width = 8, height = 8, units = "in")


