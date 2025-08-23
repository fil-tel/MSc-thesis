###############################
# Here I explore a way of being able to link tracts
# in order to obtain a tract geneaology

library(ggplot2)
library(plotly)
library(dplyr)
library(ggtree)
library(cowplot)
library(ape)
library(GenomicRanges)
library(igraph)

# devtools::install_github("bodkan/slendr", ref = "msprime-simplification-fix")
library(slendr)
init_env()

anc <- population("ancestor", time = 700e3, N = 10000, remove = 640e3)
afr <- population("AFR", parent = anc, time = 650e3, N = 10000)
nea <- population("NEA", parent = anc, time = 650e3, N = 10, remove = 39e3) # note the extremely low Ne to force coalescence
ooa <- population("OOA", parent = afr, time = 80e3, N = 5000, remove = 39e3)
eur <- population("EUR", parent = ooa, time = 40e3, N = 5000)

gf <- gene_flow(from = nea, to = ooa, rate = 0.03, start = 55000, end = 50000)

model <- compile_model(
  populations = list(anc, afr, nea, ooa, eur),
  gene_flow = gf,
  generation_time = 30,
  serialize = FALSE
)

samples <- rbind(
  schedule_sampling(model, times = 50e3, list(nea, 1)),
  schedule_sampling(model, times = c(50, 40, 30, 20, 10, 5, 0) * 1e3, list(ooa, 5), list(eur, 5)),
  schedule_sampling(model, times = 0, list(afr, 1))
)


plot_model(model, sizes = FALSE, order = c("AFR", "ancestor", "EUR", "OOA", "NEA"), proportions = TRUE)

ts <- msprime(model, sequence_length = 5e6, recombination_rate = 1e-8, samples = samples, random_seed = 123)

# get "squashed" tracts (contiguous segments merged)
squashed_tracts <- ts_tracts(ts, census = 55000, squashed = TRUE)

# Let's subset the tree sequence to only those individuals who carry introgressed
# tracts (when we plot trees below, they will be less cluttered because we'll
# be only showing individuals who are interesting, i.e. those with tracts).
subset <- c("NEA_1", "AFR_1", unique(squashed_tracts$name))
ts <- ts_simplify(ts, simplify_to = subset, filter_nodes = FALSE)

# Function which extracts a tree from a tree sequence at a given position in the genome
# and visualizes it
plot_tree <- function(ts, pos, method = c("slendr", "tskit")) {
  method <- match.arg(method)
  
  # plot the tree using the native R functionality (ape and ggtree)
  if (method == "slendr") {
    tree <- ts_phylo(ts, pos, mode = "position")
    
    p <- ggtree(tree) +
      geom_tiplab(aes(color = !grepl("NEA", label))) +
      geom_nodelab(geom = "label", size = 3) +
      theme_tree2() +
      guides(color = "none")
    revts(p) + # ggtree hack to plot x-scale correctly (and reversed)
      geom_vline(xintercept = -55000, linetype = "dashed", color = "red") + # upper limit on introgression TMRCA
      expand_limits(x = 100000) # extend x-axis to fit tree tip labels
  } else { # plot it using slendr's interface to tskit drawing functionality
    tree <- ts_tree(ts, i = pos, mode = "position")
    ts_draw(tree, labels = TRUE, width = 2000, height = 800)
  }
}

# color tracts by time -- this will give the same plot as p0
squashed_tracts %>%
  mutate(chrom = paste(name, " (node", node_id, ")"),
         chrom = factor(chrom, levels = unique(chrom[order(time)])),
         time = factor(time, levels = sort(unique(time), decreasing = TRUE))) %>%
  ggplot(aes(x = left, xend = right, y = chrom, yend = chrom, color = time)) +
  geom_segment(linewidth = 3) +
  labs(x = "position [bp]", y = "haplotype") +
  coord_cartesian(xlim = c(0, ts$sequence_length)) +
  theme(panel.grid = element_blank()) +
  ggtitle("Neanderthal tracts on sampled chromosomes") +
  facet_grid(time ~ ., scales = "free_y") -> p1

p1

ggsave("figs/tracts_all.png", plot = p1, width = 300, height = 250, units = "mm")

# example locus from the figure above to demonstrate the tree-vs-tracts idea

t_start <- 1053331
p1 <- p1+geom_vline(xintercept = t_start, color="red", lty=2)

pt <- plot_tree(ts, t_start, method = "slendr")
plot_grid(p1, pt)


# extract tracts intersecting that location

my_tracts <- squashed_tracts %>% filter(t_start-left>=0, t_start-right<0)

# sanity check that it makes sense
my_tracts %>%
  mutate(chrom = paste(name, " (node", node_id, ")"),
         chrom = factor(chrom, levels = unique(chrom[order(time)])),
         time = factor(time, levels = sort(unique(time), decreasing = TRUE))) %>%
  ggplot(aes(x = left, xend = right, y = chrom, yend = chrom, color = time)) +
  geom_segment(linewidth = 3) +
  labs(x = "position [bp]", y = "haplotype") +
  coord_cartesian(xlim = c(0, ts$sequence_length)) +
  theme(panel.grid = element_blank()) +
  ggtitle("Neanderthal tracts on sampled chromosomes (present-day at the bottom)") +
  facet_grid(time ~ ., scales = "free_y") + geom_vline(xintercept = t_start, color="red", lty=2, alpha=0.3)-> p0

p0
pt <- plot_tree(ts, t_start, method = "slendr")
p_grid <- plot_grid(p0, pt)
ggsave("figs/tracts_and_tree.png", plot = p_grid, width = 400, height = 300, units = "mm")


# function  to get the adjacency matrix from tracts
# to build tract graph
get_adj_mat <- function(pos, tract_df){
  # extract the tracts intersecting with that genomic location (pos)
  my_tracts <- tract_df %>% filter(pos-left>=0, pos-right<0) %>% as.data.frame()
  nodes <- my_tracts$node_id
  my_tracts <- my_tracts %>% select(left, right)
  # convert to matrix
  tract_mat <- my_tracts %>% as.matrix()
  # get vector of starts and ends
  starts_ends <- tract_mat %>% c() %>% unique() %>% sort()
  tract_mat_new <- apply(tract_mat, MARGIN = 1, function(row) starts_ends %in% row %>% as.numeric) %>% t
  colnames(tract_mat_new) <- starts_ends
  # matrix multiplication
  adj_mat <- tract_mat_new%*%t(tract_mat_new)
  # since I do not want the link with itself (2) I set all the values = to 2 to 0
  # adj_mat[lower.tri(adj_mat, diag = TRUE)] <- 0
  diag(adj_mat) <- 0
  colnames(adj_mat) <- rownames(adj_mat) <- nodes
  adj_mat
}

adj_mat <- get_adj_mat(pos = t_start, tract_df = my_tracts)
graph <- graph_from_adjacency_matrix(adj_mat, mode="undirected")


# function to extract a subtree for a specified node
get_subtree <- function(tree, node_id){
  # convert node on the tree from slendr to node in the object
  edge_tibble <- tree %>% as_tibble()
  # edge_tibble
  node_new <- edge_tibble %>% filter(label==as.character(node_id)) %>% .[["parent"]]
  node_new
  subtrees_l <- subtrees(tree)
  id <- sapply(subtrees_l, function(sub_tree) sub_tree$name==node_new) %>% which()
  subtrees_l[[id]]
}

# extract big tree
tree <- ts_phylo(ts, t_start, mode = "position")
t <- get_subtree(tree, 1134)

palette_tips <- c(
  "1 (NEA_1)" = "red",
  "0 (NEA_1)" = "red",
  "61 (EUR_20)" = "orange",
  "49 (EUR_14)" = "orange",
  "71 (EUR_24)" = "orange",
  "70 (EUR_24)" = "green",
  "48 (EUR_14)" = "green",
  "37 (EUR_8)" = "green",
  "27 (EUR_3)" = "green",
  "12 (OOA_6)" ="green"
)
p <- ggtree(t) +
  geom_tiplab(aes(color = label)) +
  geom_nodelab(geom = "label", size = 3) +
  theme_tree2() +
  guides(color = "none")+scale_color_manual(values = palette_tips)
p <- revts(p) + # ggtree hack to plot x-scale correctly (and reversed)
  geom_vline(xintercept = -55000,
             linetype = "dashed",
             color = "red") + # upper limit on introgression TMRCA
  expand_limits(x = 10000) # extend x-axis to fit tree tip labels
p

ggsave("figs/tree_genealogy.png", plot = p, width = 190, height = 143, units = "mm")


# save graph
V(graph)$color <- c("gray", "green",  "green",  "green",  "green", "orange",  "orange", "gray", "green",  "orange" )
png("figs/my_graph.png", height = 1000, width = 1000, units = "px")
plot(graph, edge.width=3, edge.arrow.size = 10, vertex.label.cex = 2)
dev.off()
