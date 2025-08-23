# here I compute the 2d tfs
# between old East samples and old Western samples
# to show how there are common recombination sites between them
# suggesting at least some common origin of tracts

suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
  library(here)
  library(readr)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(GenomicRanges)
  library(plyranges)
  library(cowplot)
  library(latticeExtra)
  library(sf)
  library(stars)
})

source(here::here("utils.R"))

if (!dir.exists("figs"))
  dir.create("figs")

## Empirical data

# metadata
# read all cause I also want the Asian
metadata <- read_metadata_all()
# tracts
tracts <- read_tracts_all(set = "Ancient")

# Add age and coverage to the data frame
tracts_df <- dplyr::select(metadata, sampleId, ageAverage, coverage) %>% inner_join(tracts, by = c("sampleId" = "ID"))

tracts_df$age_group <- cut(tracts_df$ageAverage, breaks = c(Inf, 30e3, 12e3, 10e3, 5e3, 2e3, 0))

group_levels <- levels(tracts_df$age_group)

tracts_df <- tracts_df %>%
  mutate(
    age_group = as.character(age_group),
    age_group = ifelse(is.na(age_group), "present-day", age_group),
    age_group = factor(age_group, levels = c("present-day", group_levels))
  )
# change name sample id to make work with function afterwards
tracts_df <- dplyr::rename(tracts_df, name = sampleId)

tracts_gr <- tracts_df %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE)

seqlengths(tracts_gr) <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)[names(seqlengths(tracts_gr))]
genome(tracts_gr) <- "hg19"

# sample that are in the binary matrix
metadata <- metadata %>% filter(sampleId %in% unique(tracts_df$name))

remove(tracts_df, tracts)


old_gr <- tracts_gr %>% filter(ageAverage > 30000)

# eur_ids <- c("Kostenki", "SII", "SIII", "SIV", "SI")
# extract 3 samples from each individuals
eur_ids <- c("Kostenki", "SII", "SIII")
asi_ids <- c("UstIshim", "Yana", "Yana2")

old_gr <- old_gr %>% filter(name %in% c(eur_ids, asi_ids))
# bin_mat <- get_bin_mat_sites_emp(tracts_gr = old_gr)

bin_mat <- get_bin_mat_sites_emp(tracts_gr = old_gr)

# d_mat <- 1/dist(t(bin_mat))

# qgraph::qgraph(d_mat, layout='spring', vsize=3)

# asian are cols (x-axis)
# european are rows (y-axis)
tfs <- get_2d_tfs(bin_mat = bin_mat,
                  pop1_ids = eur_ids,
                  pop2_ids = asi_ids)
#
# pheatmap::pheatmap(tfs, cluster_rows = FALSE, cluster_cols = FALSE)
#
# tfs[is.na(tfs)] <- 0
# cols<-function(n) {
#   colorRampPalette(c("#FFC0CB", "#CC0000"))(20)
# }
#
# p <- lattice::cloud(tfs, panel.3d.cloud = panel.3dbars, col="white",                      # white borders for bars
#                xbase = 1, ybase = 1, zlim = c(0, max(tfs)),                              # No space around the bars
#                scales = list(arrows = FALSE, just = "left"), xlab = "European", ylab = "Asian", zlab = "Count recombination sites",
#
#                col.facet = level.colors(tfs, at = do.breaks(range(tfs), 20),
#                                         col.regions = cols,                                   # color ramp for filling the bars
#                                         colors = TRUE),
#                colorkey = list(col = cols, at = do.breaks(range(tfs), 20)),
#                screen = list(z = 60, x = -120, y = -180))
#


tfs_df_sites_old <- as_tibble(tfs) %>%
  mutate(row = row_number() - 1) %>%
  pivot_longer(-row, names_to = "col", values_to = "Count") %>%
  mutate(col = as.integer(stringr::str_remove(col, "V")), mode="sites", type="anc")

p_sites_old <- ggplot(tfs_df_sites_old, aes(x = col, y = row, fill = Count)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(x = "ASI", y = "EUR") + theme(aspect.ratio = 1) + ggtitle("2D-TFS: Europeans vs. Asians", "Recombination sites approach")
p_sites_old

ggsave(plot = p_sites_old, filename = "figs/p_sites_emp.png")
# subset ------------------------------------------------------------------


bin_mat_subset <- get_bin_mat_subset_emp(tracts_gr = old_gr)
tfs_subset <- get_2d_tfs(bin_mat = bin_mat_subset,
                         pop1_ids = eur_ids,
                         pop2_ids = asi_ids)
# tfs_subset

tfs_subset_df_old <- as_tibble(tfs_subset) %>%
  mutate(row = row_number() - 1) %>%
  pivot_longer(-row, names_to = "col", values_to = "Count") %>%
  mutate(col = as.integer(stringr::str_remove(col, "V")), mode="sub", type="anc")

p_sub_old <- ggplot(tfs_subset_df_old, aes(x = col, y = row, fill = Count)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(x = "ASI", y = "EUR") + theme(aspect.ratio = 1) + ggtitle("2D-TFS: Ancient Europeans vs. Asians", "Sub-tracts approach")
p_sub_old

ggsave(plot = p_sub_old, filename = "figs/p_sub_emp.png")


# location of samples -----------------------------------------------------

set.seed(34)
metadata_old <- metadata %>% filter(sampleId %in% c(eur_ids, asi_ids)) %>% st_as_sf(coords = c("longitude", "latitude")) %>%  st_set_crs(4326) %>% st_jitter(amount = 2)

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
sf::st_agr(world) <- "constant"
# bbox <- st_as_sfc(st_bbox(c(xmin = -25, xmax = 120, ymin = 150, ymax = 150), crs = st_crs(world)))
# eurasia <- world %>% filter(continent %in% c("Europe", "Asia"))
bbox <- st_as_sfc(st_bbox(c(
  xmin = 0,
  xmax = 140,
  ymin = 10,
  ymax = 90
), crs = st_crs(world)))
eurasia <- world %>% filter(continent %in% c("Asia", "Europe"))
eurasia_cut <- st_crop(st_make_valid(eurasia), bbox)

p <- ggplot() +
  geom_sf(data = eurasia_cut) +
  geom_sf(data = metadata_old, mapping = aes(color = region, alpha=ageAverage)) +
  coord_sf(crs = 3035) + ggtitle("European and Asian samples location")+guides(color=guide_legend(title="Region"), alpha=guide_legend(title="Sample age (kya)"))+
  scale_alpha_continuous(range = c(0.2, 1))

p

ggsave(plot = p, filename = "figs/sample_map.png")

# Modern samples

# metadata
# read all cause I also want the Asian
metadata <- read_metadata_all()
metadata <- metadata %>% filter(endsWith(region, "Asia") | endsWith(region, "Europe"))
# tracts
tracts <- read_tracts_all(set = "Modern")

# Add age and coverage to the data frame
tracts_df <- dplyr::select(metadata, sampleId, ageAverage, coverage) %>% inner_join(tracts, by = c("sampleId" = "ID"))


# change name sample id to make work with function afterwards
tracts_df <- dplyr::rename(tracts_df, name = sampleId)

tracts_gr <- tracts_df %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE)

seqlengths(tracts_gr) <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)[names(seqlengths(tracts_gr))]
genome(tracts_gr) <- "hg19"

# sample that are in the binary matrix
metadata <- metadata %>% filter(sampleId %in% unique(tracts_df$name))

remove(tracts_df, tracts)

set.seed(123)
# eur_ids <- c("Kostenki", "SII", "SIII", "SIV", "SI")
# extract 3 samples from each individuals
eur_ids <- metadata %>% filter(endsWith(region, "Europe")) %>% sample_n(size = 100) %>% .[["sampleId"]] 
asi_ids <- metadata %>% filter(endsWith(region,  "Asia")) %>% sample_n(size = 100) %>% .[["sampleId"]]

tracts_gr <- tracts_gr %>% filter(name %in% c(eur_ids, asi_ids))
# bin_mat <- get_bin_mat_sites_emp(tracts_gr = old_gr)

if(!file.exists("bin_mat_sites_mod.rds")){
  bin_mat <- get_bin_mat_sites_emp(tracts_gr = tracts_gr)
  saveRDS(bin_mat, file = "bin_mat_sites_mod.rds")
}else{
  bin_mat <- readRDS("bin_mat_sites_mod.rds")
}

# asian are cols (x-axis)
# european are rows (y-axis)
tfs <- get_2d_tfs(bin_mat = bin_mat,
                  pop1_ids = eur_ids,
                  pop2_ids = asi_ids)


tfs_df_sites_mod <- as_tibble(tfs) %>%
  mutate(row = row_number() - 1) %>%
  pivot_longer(-row, names_to = "col", values_to = "Count") %>%
  mutate(col = as.integer(stringr::str_remove(col, "V")), mode="sites", type="mod")

p_sites_mod <- ggplot(tfs_df_sites_mod, aes(x = col, y = row, fill = Count)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(x = "ASI", y = "EUR") + theme(aspect.ratio = 1) + ggtitle("2D-TFS: Present-day Europeans vs. Asians", "Recombination sites approach")
p_sites_mod

ggsave(plot = p_sites_mod, filename = "figs/p_sites_emp_modern.png")
# subset ------------------------------------------------------------------

if(!file.exists("bin_mat_sub_mod.rds")){
  bin_mat_subset <- get_bin_mat_subset_emp(tracts_gr = tracts_gr)
  saveRDS(bin_mat_subset, file = "bin_mat_sub_mod.rds")
}else{
  bin_mat_subset <- readRDS("bin_mat_sub_mod.rds")
}

tfs_subset <- get_2d_tfs(bin_mat = bin_mat_subset,
                         pop1_ids = eur_ids,
                         pop2_ids = asi_ids)
# tfs_subset

tfs_subset_df_mod <- as_tibble(tfs_subset) %>%
  mutate(row = row_number() - 1) %>%
  pivot_longer(-row, names_to = "col", values_to = "Count") %>%
  mutate(col = as.integer(stringr::str_remove(col, "V")), mode="sub", type="mod")

p_sub_mod <- ggplot(tfs_subset_df_mod, aes(x = col, y = row, fill = Count)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(x = "ASI", y = "EUR") + theme(aspect.ratio = 1) + ggtitle("2D-TFS: Europeans vs. Asians", "Sub-tracts approach")
p_sub_mod

ggsave(plot = p_sub_mod, filename = "figs/p_sub_emp_modern.png")

# Grid

p_grid <- cowplot::plot_grid(p_sub_mod, p_sites_mod, labels = "AUTO")
ggsave(
  plot = p_grid,
  filename = "figs/p_grid_emp.png",
  width = 10,
  height = 5,
  units = "in"
)



# plot modern and ancient together

# cretae df

full_df <- rbind(tfs_df_sites_mod, tfs_df_sites_old, tfs_subset_df_mod, tfs_subset_df_old)
# 
# full_df %>% ggplot() + geom_tile(aes(x = col, y = row, fill = Count)) + facet_wrap(type~mode, scales = "free")+
#   scale_fill_viridis_c() +
#   labs(x = "ASI", y = "EUR") + theme(aspect.ratio = 1) + ggtitle("2D-TFS: Europeans vs. Asians", "Sub-tracts approach")


plots <- full_df %>%
  group_split(type, mode) %>%
  map(~{
    group_label <- ifelse(unique(.x$type)=="anc", "Ancient samples", "Present-day samples")
    approach <- ifelse(unique(.x$mode)=="sites", "Recombination sites approach", "Sub-tracts approach")
    ggplot(.x) +
      geom_tile(aes(x = col, y = row, fill = Count)) +
      scale_fill_viridis_c() +
      # scale_fill_gradient(low = "lightblue", high = "blue")+
      labs(x = "ASI", y = "EUR", fill = "Count") +
      ggtitle(paste("2D-TFS:", group_label), subtitle = approach) +
      theme(aspect.ratio = 1)
  })

p_def <- cowplot::plot_grid(plotlist = plots, ncol = 2, labels = "AUTO")

ggsave(plot = p_def, filename = "figs/plot_mod_anc_euras.png")

