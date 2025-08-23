# here i apply my tract genealogy graph
# on empirical spatial sample
# this is mainly an example and exploration rather than a proper analysis

suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
  library(here)
  library(readr)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(GenomicRanges)
  library(plyranges)
  library(cowplot)
  library(sf)
  library(sfnetworks)
  library(igraph)
  library(ggraph)
})

source(here::here("utils.R"))

if(!dir.exists("figs")) dir.create("figs")

# get metadata
metadata <- read_metadata_all() %>% 
  filter(region %in% c("SouthernEurope", "WesternEurope", "WesternAsia", "NorthernEurope", "CentralEasternEurope",
                       "WesternAsia", "NorthAsia", "SouthAsia", "CentralAsia", "SouthEastAsia"))

tracts_df <- read_tracts_all(set="Ancient")
# Add age and coverage to the data frame
tracts_df <- dplyr::select(metadata, sampleId, ageAverage, coverage) %>% inner_join(tracts_df, by = c("sampleId" = "ID"))

tracts_df$age_group <- cut(
  tracts_df$ageAverage,
  breaks = c(Inf, 30e3, 12e3, 10e3, 5e3, 2e3, 0)
)

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

tracts_gr_selected <- tracts_gr %>% filter(age_group=="(1e+04,1.2e+04]")

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
sf::st_agr(world) <- "constant"
# bbox <- st_as_sfc(st_bbox(c(xmin = -25, xmax = 120, ymin = 150, ymax = 150), crs = st_crs(world)))
# eurasia <- world %>% filter(continent %in% c("Europe", "Asia"))
bbox <- st_as_sfc(st_bbox(c(xmin = -25, xmax = 65, ymin = 25, ymax = 70), crs = st_crs(world)))
western_eurasia <- st_crop(st_make_valid(world), bbox)

metadata %>% filter(sampleId %in% unique(tracts_gr_selected$name), ) %>% 
  filter(!is.na(latitude) & !is.na(longitude)) %>%
  st_as_sf(coords = c("longitude", "latitude")) %>%
  st_set_crs(4326) %>% st_jitter(amount = 1.5) %>% 
  ggplot() +
  geom_sf(data = western_eurasia) +
  geom_sf(aes(color = sampleId, size=ageAverage)) +
  coord_sf(crs = 3035) + guides(color="none")

bin_mat_sites <- get_bin_mat_sites_emp(tracts_gr_selected)
# bin_mat_subset <- get_bin_mat_subset_emp(tracts_gr = tracts_gr_selected)

# get a location that all of them have from the subset map

# bin_mat_subset[rowSums(bin_mat_subset)==ncol(bin_mat_subset)-1,]
# select location

location <- GRanges(seqnames = "chr4", ranges = IRanges(start = 35491388 , end = 35524710))
# select tracts intersecting that location to compute genealogy
tracts_gen <- tracts_gr_selected[queryHits(findOverlaps(tracts_gr_selected, location))]

get_adj_mat_emp <- function(tract_df){
  # extract the tracts intersecting with that genomic location (pos)
  my_tracts <- tract_df %>% as.data.frame()
  names <- my_tracts$name
  my_tracts <- my_tracts %>% select(start, end)
  # convert to matrix
  tract_mat <- my_tracts %>% as.matrix()
  # get vector of starts and ends
  starts_ends <- tract_mat %>% c() %>% unique() %>% sort()
  tract_mat_new <- apply(tract_mat, MARGIN = 1, function(row) starts_ends %in% row %>% as.numeric) %>% t
  colnames(tract_mat_new) <- starts_ends
  # matrix multiplication
  adj_mat <- tract_mat_new%*%t(tract_mat_new)
  # since I do not want the link with itself (2) I set all the values = to 2 to 0
  adj_mat[lower.tri(adj_mat, diag = TRUE)] <- 0
  diag(adj_mat) <- 0
  colnames(adj_mat) <- rownames(adj_mat) <- names
  adj_mat
}

# t_start <- 35491388
adj_mat <- get_adj_mat_emp(tract_df = tracts_gen)
gg <- igraph::graph_from_adjacency_matrix(adj_mat, mode = "undirected")
plot(gg)


# sfnetworks --------------------------------------------------------------

# convert cooridnates to sf points

metadata_2_plot <- metadata %>% filter(sampleId %in% unique(tracts_gr_selected$name))
set.seed(123)
metadata_sf <- metadata_2_plot %>% st_as_sf(coords = c("longitude", "latitude")) %>%  st_set_crs(4326) %>% st_jitter(amount = 2)

# get edges from graph
edge_df <- igraph::as_edgelist(gg) %>% as.data.frame()
colnames(edge_df) <- c("from", "to")

net <- sfnetwork(nodes = metadata_sf, edges = edge_df, directed = FALSE)
# 
# ggraph(net, layout = "manual", node.position = st_coordinates(metadata_sf)) +
#   geom_edge_link(color = "steelblue", alpha = 0.6) +
#   geom_node_point(aes(color = region), size = 3) +
#   theme_minimal()


world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
sf::st_agr(world) <- "constant"
# bbox <- st_as_sfc(st_bbox(c(xmin = -25, xmax = 120, ymin = 150, ymax = 150), crs = st_crs(world)))
# eurasia <- world %>% filter(continent %in% c("Europe", "Asia"))
bbox <- st_as_sfc(st_bbox(c(xmin = -10, xmax = 60, ymin = 25, ymax = 70), crs = st_crs(world)))
western_eurasia <- st_crop(st_make_valid(world), bbox)

# activate nodes
net = net %>%
  activate("nodes") 
net <- tidygraph::convert(
  net, 
  to_spatial_explicit, 
  .clean = TRUE
)

# metadata_sf %>%

edge_df <- as_data_frame(net, what = "edges")

edge_df <- edge_df %>%
  mutate(pair = paste0(pmin(from, to), "_", pmax(from, to))) %>%
  add_count(pair, name = "pair_count") %>%  
  mutate(is_identical = ifelse(pair_count > 1, "Yes", "No"))

net <- net %>%
  activate("edges") %>%
  mutate(
    pair = edge_df$pair,
    pair_count = edge_df$pair_count,
    is_identical = edge_df$is_identical
  )


p <- ggplot() +
  geom_sf(data = western_eurasia) +
  geom_sf(data = st_as_sf(net, "edges"), mapping = aes(linetype= is_identical)) + 
  geom_sf(data = metadata_sf, mapping = aes(color=sampleId, size = ageAverage))+guides(color=guide_legend(title="Sample ID"), linetype=guide_legend(title="Identical tract"))+
  coord_sf(crs = 3035) +   scale_size_continuous(
    name = "Sample age",  
    range = c(1, 3)
  )

p

ggsave(plot = p, filename = "figs/spatial_graph.png", height = 10, width = 13)
