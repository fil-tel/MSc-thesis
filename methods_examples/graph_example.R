library(GenomicRanges)
library(tidyverse)
library(plyranges)
library(igraph)
library(cowplot)
source("utils.R")

if(!dir.exists("figs")) dir.create("figs")

# define tracts

tracts_gr <- GRanges(seqnames = "chr1", ranges = IRanges(start = c(100, 100, 130, 90, 90, 50), end = c(150, 150, 150, 170, 140, 80)))
mcols(tracts_gr)$name <- paste0("IND_", 1:length(tracts_gr))
mcols(tracts_gr)$tract_id <- paste0("T_", 1:length(tracts_gr))

p_ind <- setNames(RColorBrewer::brewer.pal(length(tracts_gr), "Set1"), tracts_gr$tract_id)


# build graph

get_adj_mat_ex <- function(tract_gr){
  # extract the tracts intersecting with that genomic location (pos)
  my_tracts <- tract_gr %>% as.data.frame()
  nodes <- my_tracts$tract_id
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
  colnames(adj_mat) <- rownames(adj_mat) <- nodes
  adj_mat
}

adj_mat <- get_adj_mat_ex(tracts_gr)
graph <- graph_from_adjacency_matrix(adj_mat, mode = "undirected")


V(graph)$color <- p_ind

png("figs/my_graph.png", height = 1000, width = 1000, units = "px")
plot(graph, edge.width = 4, vertex.label.cex = 2, vertex.size=30)
dev.off()

graph_1_cwpl <- ggdraw() + draw_image("figs/my_graph.png")

links_1 <- data.frame(ystart = c("T_1", "T_1", "T_2"),
                    yend = c("T_2", "T_3", "T_3"),
                    xstart=c( 150, 150, 150),
                    xend = c( 150, 150, 150))

links_2 <- data.frame(ystart = c("T_1", "T_4"),
                      yend = c("T_2", "T_5"),
                      xstart=c(100, 90),
                      xend = c(100, 90))
  
tracts_gr %>% as_tibble %>% 
  ggplot() +
  geom_segment(aes(x = start, xend = end, y = tract_id, yend = tract_id, colour = tract_id),linewidth = 3) + scale_y_discrete(limits=rev) +
  labs(x = "Position [bp]", y = "Samples") +
  coord_cartesian(xlim = c(0, 200)) +
  # scale_x_continuous(breaks = tick_positions) +
  # theme(panel.grid = element_blank()) + 
  guides(colour = "none")+
  ggtitle("Example of introgressed tracts on chromosomes")+
  scale_colour_manual(values = p_ind) + 
  geom_curve(data=links_1,
    aes(x = xstart, y = ystart, xend = xend, yend = yend), curvature = -0.8) +
  geom_curve(data=links_2,
             aes(x = xstart, y = ystart, xend = xend, yend = yend), curvature = 0.8) -> p1

ggsave(plot = p1, filename = "figs/edges_on_plot.png")

grid_plot <- plot_grid(p1, graph_1_cwpl, labels = "AUTO")

ggsave(plot = grid_plot, filename = "figs/graph_grid.png", width = 11, height = 5, units = "in")

# latex table
# df <- tracts_gr %>% as.data.frame() %>% rename("chrom"=seqnames) %>% select(-c(strand, width))
# print(xtable(df, caption = "Test", label = "method:tab_tracts"), include.rownames=FALSE)

bmatrix = function(x, digits=NULL, ...) {
  library(xtable)
  default_args = list(include.colnames=TRUE, only.contents=TRUE,
                      include.rownames=TRUE, hline.after=NULL, comment=FALSE,
                      print.results=FALSE)
  passed_args = list(...)
  calling_args = c(list(x=xtable(x, digits=digits)),
                   c(passed_args,
                     default_args[setdiff(names(default_args), names(passed_args))]))
  cat("\\begin{bmatrix}\n",
      do.call(print.xtable, calling_args),
      "\\end{bmatrix}\n")
}

bmatrix(tract_mat_new, digits = 0)

