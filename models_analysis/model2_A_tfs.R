library(dplyr)
library(readr)
library(parallel)
library(ggplot2)
library(slendr)
init_env()
source("utils.R")


# model 1 A ---------------------------------------------------------------
if(!dir.exists("tfs_df")) dir.create("tfs_df")
if(!dir.exists("rds_files")) dir.create("rds_files")

suffix <- "2_A"
dir_m <- sprintf("model%s/", suffix)
ncores <- detectCores()-1

tracts_files <- list.files(path = dir_m, pattern =  "*tracts.tsv.gz", full.names = TRUE)
names_list <- sapply(strsplit(tracts_files, "_"), function(x) paste0("Ne_", paste0(x[4:5], collapse = "_"))) %>% gsub(pattern=".tsv", replacement="")
tracts_list <- setNames(lapply(tracts_files, read_tsv), names_list)
tracts_list_gr <- lapply(tracts_list, convert_2_gr)

# subset
if (!file.exists(sprintf("rds_files/bin_mat_sub_%s.rds", suffix))) {
  bin_mat_list_sub <- mclapply(tracts_list_gr, get_bin_mat_subset, mc.cores = ncores)
  saveRDS(bin_mat_list_sub, sprintf("rds_files/bin_mat_sub_%s.rds", suffix))
} else {
  bin_mat_list_sub <- readRDS(sprintf("rds_files/bin_mat_sub_%s.rds", suffix))
}

# unique
if (!file.exists(sprintf("rds_files/bin_mat_uni_%s.rds", suffix))) {
  bin_mat_list_uni <- mclapply(tracts_list_gr, get_bin_mat_unique, mc.cores = ncores)
  saveRDS(bin_mat_list_uni, sprintf("rds_files/bin_mat_uni_%s.rds", suffix))
} else {
  bin_mat_list_uni <- readRDS(sprintf("rds_files/bin_mat_uni_%s.rds", suffix))
}

# window
if (!file.exists(sprintf("rds_files/bin_mat_wind_%s.rds", suffix))) {
  bin_mat_list_wind <- mclapply(tracts_list_gr, function(x) get_bin_mat_windows(tracts_gr = x, len_chr = 100e6), mc.cores = ncores)
  saveRDS(bin_mat_list_wind, sprintf("rds_files/bin_mat_wind_%s.rds", suffix))
} else {
  bin_mat_list_wind <- readRDS(sprintf("rds_files/bin_mat_wind_%s.rds", suffix))
}

# sites
if (!file.exists(sprintf("rds_files/bin_mat_sites_%s.rds", suffix))) {
  bin_mat_list_sites <- mclapply(tracts_list_gr, function(x) get_bin_mat_sites(tracts_gr = x), mc.cores = ncores)
  saveRDS(bin_mat_list_sites, sprintf("rds_files/bin_mat_sites_%s.rds", suffix))
} else {
  bin_mat_list_sites <- readRDS(sprintf("rds_files/bin_mat_sites_%s.rds", suffix))
}

# to pass a list of binary matrix
# return a df with all the rthe TFS for each replicate and Ne 
# for the model
get_df_2_plot <- function(bin_mat_list){
  bin_mat_list <- Filter(Negate(is.null), bin_mat_list)
  tfs_df <- as.data.frame(sapply(bin_mat_list, compute_tfs))
  tfs_df$bin <- 1:nrow(tfs_df)
  
  # split it 
  ne_groups <- unique(sub("(_rep\\d+)$", "", names(tfs_df)))
  mat_list <- lapply(ne_groups, function(ne) {
    tfs_df[ , grepl(paste0("^", ne, "_rep"), names(tfs_df)), drop = FALSE]
  })
  names(mat_list) <- ne_groups
  
  
  long_df <- tfs_df %>%
    pivot_longer(
      cols = starts_with("Ne_"),
      names_to = c("Ne", "rep"),
      names_pattern = "Ne_(\\d+)_rep(\\d+)",
      values_to = "value"
    ) %>%
    mutate(
      Ne = as.integer(Ne),
      rep = as.integer(rep)
    )
  
  summary_df <- long_df %>%
    group_by(bin, Ne) %>%
    summarise(
      mean = mean(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      .groups = "drop"
    )
  summary_df
}

sub_df <- get_df_2_plot(bin_mat_list = bin_mat_list_sub)
sub_df$approach <- "sub"
sites_df <- get_df_2_plot(bin_mat_list = bin_mat_list_sites)
sites_df$approach <- "sites"
wind_df <- get_df_2_plot(bin_mat_list = bin_mat_list_wind)
wind_df$approach <- "wind"
uni_df <- get_df_2_plot(bin_mat_list = bin_mat_list_uni)
uni_df$approach <- "uni"

# get a big df storing all of them
# to facilitate the plotting in another scriptr
df2save <- rbind(sub_df, wind_df, uni_df, sites_df)
df2save$model <- sprintf("model%s", suffix)

# save the df
write_tsv(df2save, file = sprintf("tfs_df/tfs_model_%s.tsv", suffix))

# p <- ggplot(wind_df, aes(x = bin, y = mean, color = factor(Ne), fill = factor(Ne))) +
#   facet_wrap(~Ne)+
#   geom_line() +
#   # geom_point() +
#   geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.2, color = NA) +
#   labs(
#     x = "Bin",
#     y = "Density",
#     color = "Ne",
#     fill = "Ne",
#     title = "TFS "
#   ) +
#   theme_minimal()
# p
# 
# ggsave(plot = p, filename = "exxxx.png")




# extra -------------------------------------------------------------------


# plot_tfs_line <- function(df){
#   df$bin <- 1:nrow(df)
#   
#   summary_df <- df %>%
#     rowwise() %>%
#     mutate(mean = mean(c_across(starts_with("Ne_"))),
#            sd = sd(c_across(starts_with("Ne_")))) %>%
#     ungroup() %>%
#     select(bin, mean, sd)
#   # half sd cause to clutterd
#   ggplot(summary_df, aes(x = bin, y = mean)) +
#     geom_line(color = "red") +
#     # geom_point(color = "red") +
#     geom_ribbon(aes(ymin = mean - sd/2, ymax = mean + sd/2), alpha = 0.1, fill = "red") +
#     labs(x = "Bin", y = "Density", title = "TFS") +
#     theme_minimal()  
# }
# 
# plot_tfs_line(mat_list$Ne_2000)


