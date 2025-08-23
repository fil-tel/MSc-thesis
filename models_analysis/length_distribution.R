library(readr)
library(ggplot2)
library(cowplot)
library(parallel)
library(dplyr)
source("utils.R")

if(!dir.exists("figs_length")) dir.create("figs_length")
if(!dir.exists("length_df")) dir.create("length_df")

suffices <- c(paste0(1:3, c("_A")), paste0(1:3, c("_B")))


get_length_df <- function(length_vec){
  h_ind <- length_vec %>% hist(breaks = seq(0, 7e6, 20000), plot = FALSE)
  length_df <- data.frame(x = as.integer(h_ind$mids), y = h_ind$density, line = dexp(h_ind$mids, rate = 1 / mean(length_vec)))
  length_df
}

get_length_2_plot <- function(suffix) {
  dir_m <- sprintf("model%s/", suffix)
  ncores <- detectCores()-1
  tracts_files <- list.files(path = dir_m,
                             pattern =  "*tracts.tsv.gz",
                             full.names = TRUE)
  names_list <- sapply(strsplit(tracts_files, "_"), function(x)
    paste0("Ne_", paste0(x[4:5], collapse = "_"))) %>% gsub(pattern = ".tsv", replacement =
                                                              "")
  tracts_list <- setNames(mclapply(tracts_files, read_tsv, mc.cores = ncores), names_list)
  tracts_list_gr <- lapply(tracts_list, convert_2_gr)
  if(unlist(strsplit(suffix, split = "_"))[2]=="B") {
    tracts_list_gr <- mclapply(tracts_list_gr, function(tracts_gr) {
      tracts_grl <- lapply(unique(tracts_gr$name), function(id) {
        # extract tract for individual
        ind_gr <- tracts_gr %>% filter(name == id)
        # create a df for the metacolumns, they will be lost otherwise
        # and remove some columns because not needed or wrong
        df <- mcols(ind_gr) %>% as.data.frame() %>% select(-one_of("node_id", "length", "haplotype"))
        # reduce tracts to simuate what IBDmix would do
        # since we do not have the haplotype in output
        res_gr <- reduce(ind_gr)
        # add metadata colummns
        mcols(res_gr) <- df[1:length(res_gr), ]
        res_gr$length <- width(res_gr)
        res_gr
      })
      do.call(c, tracts_grl)
    }, mc.cores = ncores)
  }
  # tracts_list_gr
  lengths <- lapply(tracts_list_gr, function(x) x$length)
  indices <- unique(gsub(pattern = "rep(\\d+)", replacement = "", names(tracts_list_gr)))
  df <- do.call(rbind, lapply(indices, function(ind) {
    res <- get_length_df(unlist(lengths[grepl(pattern = ind, names(tracts_list_gr))]))
    res$ne <- as.numeric(unlist(strsplit(ind, split = "_"))[2])
    res
  }))
  sprintf("model_%s is done\n", suffix)
  df$model <- sprintf("model_%s", suffix)
  df
}


if (!file.exists("length_df/lengths_df_no_filter.tsv")) {
  final_df_A <- do.call(rbind, lapply(suffices[1:3], get_length_2_plot))
  final_df_B <- do.call(rbind, lapply(suffices[4:6], get_length_2_plot))
  final_df <- rbind(final_df_A, final_df_B)
  write_tsv(final_df, file = "length_df/lengths_df_no_filter.tsv")
}

final_df <- read_tsv("length_df/lengths_df_no_filter.tsv")

p_models <- ggplot() + geom_line(
  data = final_df,
  mapping = aes(x, y, color = factor(ne)),
  linetype = 2
) + facet_wrap( ~ model, ncol = 2)+
  xlab("Tract length bin [bp]") + ylab("Density") + theme(aspect.ratio = 1) +
  labs(color = "Ne OOA", linetype="")+coord_cartesian(xlim = c(0, 2e6))

p_models

ggsave(plot = p_models, filename = "figs_length/length_dist_models.png", width = 7)

p_length_cut <- ggplot() + geom_line(
  data = final_df,
  mapping = aes(x, y, color = factor(ne)),
  linetype = 2
) +  facet_wrap( ~ model, ncol = 2)+
  xlab("Tract length bin [bp]") + ylab("Density") + theme(aspect.ratio = 1) +
  labs(color = "Ne OOA", linetype="")+coord_cartesian(xlim = c(0, 500e3))

p_length_cut+theme(legend.position = "none")

ggsave(plot = p_length_cut, filename = "figs_length/length_dist_cut.png")

p_length_cut <- p_length_cut

p_models <- p_models+theme(legend.position = "none")

cowplot::plot_grid(p_models, p_length_cut, labels = "AUTO")

ggsave(filename = "figs_length/length_grid.png", width = 16, height = 10)

