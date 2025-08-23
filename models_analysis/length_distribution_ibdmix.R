library(BSgenome.Hsapiens.UCSC.hg19)
library(readr)
library(ggplot2)
library(cowplot)
library(parallel)
library(dplyr)
source("utils.R")

if(!dir.exists("figs_length")) dir.create("figs_length")
if(!dir.exists("length_df")) dir.create("length_df")

suffices <- c(paste0(1:3, c("_A")), paste0(1:3, c("_B")))
# size to filter out tracts
size_cutoff <- 50e3


get_length_df <- function(length_vec){
  h_ind <- length_vec %>% hist(breaks = seq(49999, 7e6, 20000), plot = FALSE)
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
  tracts_list_gr <- mclapply(tracts_list_gr, function(tracts_gr) {
      tracts_gr$name <- gsub(pattern = "_hap_\\d+", replacement = "", x = tracts_gr$name)
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
        filter(res_gr, length>=size_cutoff)
      })
      do.call(c, tracts_grl)
    }, mc.cores = ncores)
  
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


if (!file.exists("length_df/lengths_df_ibdmix_like.tsv")) {
  final_df_A <- do.call(rbind, lapply(suffices[1:3], get_length_2_plot))
  final_df_B <- do.call(rbind, lapply(suffices[4:6], get_length_2_plot))
  final_df <- rbind(final_df_A, final_df_B)
  write_tsv(final_df, file = "length_df/lengths_df_ibdmix_like.tsv")
}

final_df <- read_tsv("length_df/lengths_df_ibdmix_like.tsv")


p_models <- ggplot() + geom_line(
  data = filter(final_df, !x %in% c(10000, 30000)),
  mapping = aes(x, y, color = factor(ne)),
  linetype = 2
) + facet_wrap( ~ model, ncol = 2)+
  xlab("Tract length bin [bp]") + ylab("Density") + theme(aspect.ratio = 1) +
  labs(title = "Distribution of Tract Lengths", subtitle = "Models", color = "Ne OOA", linetype="")

p_models
# empirical data ----------------------------------------------------------


# metadata
metadata <- read_metadata()
# tracts
# by default is only western Eurasian
tracts <- read_tracts(set = "Modern")

# Add age and coverage to the data frame
tracts_df <- dplyr::select(metadata, sampleId, ageAverage, coverage) %>% inner_join(tracts, by = c("sampleId" = "ID"))

# change name sample id to make work with function afterwards
tracts_df <- dplyr::rename(tracts_df, name = sampleId)


# sample 20 ppl
# sample_ids <- sample(unique(tracts_df$name), 20)
# tracts_df <- tracts_df %>% filter(name %in% sample_ids)

tracts_gr <- tracts_df %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE)

seqlengths(tracts_gr) <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)[names(seqlengths(tracts_gr))]
genome(tracts_gr) <- "hg19"

metadata <- metadata %>% filter(sampleId %in% unique(tracts_df$name))

# tracts_gr_filt <- tracts_gr %>% filter(name %in% sample(unique(tracts_gr$name), 50))
emp_df <- get_length_df(tracts_gr$length)
# emp_df <- get_length_df(tracts_gr_filt$length)


p_length_full <- ggplot() + geom_line(
  data = filter(final_df, !x %in% c(10000, 30000)),
  mapping = aes(x, y, color = factor(ne)),
  linetype = 2
) + geom_line(data = filter(emp_df, !x %in% c(10000, 30000)), aes(x = x, y = y,  linetype = "Empirical data")) +
  facet_wrap( ~ model, ncol = 2)+scale_linetype_manual(values = c("Empirical data" = "dashed"))+
  xlab("Tract length bin [bp]") + ylab("Density") + theme(aspect.ratio = 1) +
  labs(color = "Ne OOA", linetype="")+coord_cartesian(xlim = c(50e3, 5e6))

p_length_full

ggsave(plot = p_length_full, filename = "figs_length/length_dist_emp_full.png")


p_length_cut <- ggplot() + geom_line(
  data = filter(final_df, !x %in% c(10000, 30000)),
  mapping = aes(x, y, color = factor(ne)),
  linetype = 2
) + geom_line(data = filter(emp_df, !x %in% c(10000, 30000)), aes(x = x, y = y,  linetype = "Empirical data")) +
  facet_wrap( ~ model, ncol = 2)+scale_linetype_manual(values = c("Empirical data" = "dashed"))+
  xlab("Tract length bin [bp]") + ylab("Density") + theme(aspect.ratio = 1) +
  labs(color = "Ne OOA", linetype="")+coord_cartesian(xlim = c(50e3, 5e5))

p_length_cut+theme(legend.position = "none")

ggsave(plot = p_length_cut, filename = "figs_length/length_dist_emp_cut.png")

p_length_cut <- p_length_cut

p_length_full <- p_length_full+theme(legend.position = "none")

cowplot::plot_grid(p_length_full, p_length_cut, labels = "AUTO")

ggsave(filename = "figs_length/length_grid_emp.png", width = 16, height = 10)

