#########################
# In this script I explore the joint tract frequency spectrum
# aka frequency spectrum for shared introgressed tracts
# for a population that experience a single pulse of introgression and
# subsequently split into two populations

library(ggplotify)
library(dplyr)
library(ggplot2)
library(cowplot)
source("utils.R")

# devtools::install_github("bodkan/slendr")

library(slendr)
init_env()

if(!dir.exists("figs")) dir.create("figs")
set.seed(1234)

# define Ne receiving population
Ne <- 3000

anc <- population("ancestor",
                  time = 700e3,
                  N = 10000,
                  remove = 640e3)
afr <- population(
  "AFR",
  parent = anc,
  time = 650e3,
  N = 10000
)
nea <- population(
  "NEA",
  parent = anc,
  time = 650e3,
  N = 1000,
  remove = 39e3
)
ooa <- population(
  "OOA",
  parent = afr,
  time = 70e3,
  N = Ne,
  remove = 49e3
)
eur1 <- population(
  "EUR1",
  parent = ooa,
  time = 51e3,
  N = Ne
)

eur2 <- population(
  "EUR2",
  parent = ooa,
  time = 51e3,
  N = Ne
)

# 3% introgressino into OOA
gf <- gene_flow(
  from = nea,
  to = ooa,
  rate = 0.03,
  start = 60000,
  end = 55000,
  overlap = FALSE
)

model <- compile_model(
  populations = list(anc, afr, nea, ooa, eur1, eur2),
  gene_flow = gf,
  generation_time = 30,
  serialize = FALSE,
  time_units = "years before present"
)

schedule <- schedule_sampling(
    model,
    times = c(50, 40, 30, 20, 10, 0) * 1e3,
    list(eur1, 20),
    list(eur2, 20)
  )

plot_model(
  model,
  order = c("AFR", "ancestor",   "EUR2", "EUR1", "OOA", "NEA"),
  proportions = TRUE,
  file = "figs/model_1p.png"
)
# define number of rows and column for the final matrx
n_row <- unique(schedule$n)*2+1
# define list tp store results
time_list_subset <- lapply(unique(schedule$time), function(x) matrix(0, nrow = n_row, ncol = n_row))
time_list_sites <- lapply(unique(schedule$time), function(x) matrix(0, nrow = n_row, ncol = n_row))
time_list_windows <- lapply(unique(schedule$time), function(x) matrix(0, nrow = n_row, ncol = n_row))
time_list_unique <- lapply(unique(schedule$time), function(x) matrix(0, nrow = n_row, ncol = n_row))
names(time_list_sites) <- names(time_list_subset) <- names(time_list_unique) <- names(time_list_windows) <- unique(schedule$time)

i <- 1
while (i<6) {
ts <- msprime(
  model,
  sequence_length = 10e6,
  recombination_rate = 1e-8,
  samples = schedule
)

tracts <- ts_tracts(ts, census = 60000, squashed = TRUE)
# add chromosome column
tracts$chrom <- "chr1"

tracts <- tracts %>% mutate(name = paste0(name, "_hap_", haplotype))

tracts_gr <- makeGRangesFromDataFrame(
  tracts,
  keep.extra.columns = TRUE,
  seqnames.field = "chrom",
  start.field = "left",
  end.field = "right"
)

# build 2d TFS
# and store it in list according to time
for(time_point in names(time_list_subset)){
  samples <- tracts_gr %>% filter(time==as.numeric(time_point))
  bin_mat <- get_bin_mat_subset(samples)
  eur1_ids <- as_tibble(samples) %>% filter(pop=="EUR1") %>% .[["name"]] %>% unique()
  eur2_ids <- as_tibble(samples) %>% filter(pop=="EUR2") %>% .[["name"]] %>% unique()
  res <- get_2d_tfs(bin_mat = bin_mat, pop1_ids = eur1_ids, pop2_ids = eur2_ids)
  res[is.na(res)] <- 0 
  time_list_subset[[time_point]][1:nrow(res), 1:ncol(res)] <-time_list_subset[[time_point]][1:nrow(res), 1:ncol(res)]+res
}

for(time_point in names(time_list_subset)){
  samples <- tracts_gr %>% filter(time==as.numeric(time_point))
  bin_mat <- get_bin_mat_sites(samples)
  eur1_ids <- as_tibble(samples) %>% filter(pop=="EUR1") %>% .[["name"]] %>% unique()
  eur2_ids <- as_tibble(samples) %>% filter(pop=="EUR2") %>% .[["name"]] %>% unique()
  res <- get_2d_tfs(bin_mat = bin_mat, pop1_ids = eur1_ids, pop2_ids = eur2_ids)
  res[is.na(res)] <- 0 
  time_list_sites[[time_point]][1:nrow(res), 1:ncol(res)] <-time_list_sites[[time_point]][1:nrow(res), 1:ncol(res)]+res
}

for(time_point in names(time_list_windows)){
  samples <- tracts_gr %>% filter(time==as.numeric(time_point))
  bin_mat <- get_bin_mat_windows(tracts_gr = samples, window_size = 50e3, step_size = 50e3, len_chr = 10e6)
  eur1_ids <- as_tibble(samples) %>% filter(pop=="EUR1") %>% .[["name"]] %>% unique()
  eur2_ids <- as_tibble(samples) %>% filter(pop=="EUR2") %>% .[["name"]] %>% unique()
  res <- get_2d_tfs(bin_mat = bin_mat, pop1_ids = eur1_ids, pop2_ids = eur2_ids)
  res[is.na(res)] <- 0 
  time_list_windows[[time_point]][1:nrow(res), 1:ncol(res)] <-time_list_windows[[time_point]][1:nrow(res), 1:ncol(res)]+res
}

for(time_point in names(time_list_unique)){
  samples <- tracts_gr %>% filter(time==as.numeric(time_point))
  bin_mat <- get_bin_mat_unique(tracts_gr = samples)
  eur1_ids <- as_tibble(samples) %>% filter(pop=="EUR1") %>% .[["name"]] %>% unique()
  eur2_ids <- as_tibble(samples) %>% filter(pop=="EUR2") %>% .[["name"]] %>% unique()
  res <- get_2d_tfs(bin_mat = bin_mat, pop1_ids = eur1_ids, pop2_ids = eur2_ids)
  res[is.na(res)] <- 0 
  time_list_unique[[time_point]][1:nrow(res), 1:ncol(res)] <- time_list_unique[[time_point]][1:nrow(res), 1:ncol(res)]+res
}

i <- i+1
}

#labbeller 
my_labeller <- function(time){
  return(paste0("Sample time = ", time, " ybp"))
}


# subset approach ---------------------------------------------------------


# get long df for plotting of subset approach
df_subset <- purrr::imap_dfr(time_list_subset, function(mat, sample_time) {
  mat[mat==0] <- NA
  as_tibble(mat) %>%
    mutate(row = row_number()-1) %>%
    pivot_longer(-row, names_to = "col", values_to = "Count") %>% 
    mutate(col = as.integer(stringr::str_remove(col, "V")),
           time = sample_time)
})

df_subset$time <- factor(df_subset$time, levels = rev(unique(schedule$time))) 

saveRDS(df_subset, file = "df_subset_1p.rds")
df_subset <- readRDS("df_subset_1p.rds")

p_subset <- ggplot(df_subset, aes(x = col, y = row, fill = Count)) +
  geom_tile() +
  scale_fill_viridis_c() +
  facet_wrap(.~time, labeller = labeller(time = my_labeller)) +
  labs(x = "EUR_1", y = "EUR_2")+theme(aspect.ratio=1)+ggtitle("Model X - sub-tracts approach")

ggsave("figs/2d_tfs_1p_subset.png", plot = p_subset, width = 190, height = 143, units = "mm")

# recombination sites approach --------------------------------------------

df_sites <- purrr::imap_dfr(time_list_sites, function(mat, sample_time) {
  mat[mat==0] <- NA
  as_tibble(mat) %>%
    mutate(row = row_number()-1) %>%
    pivot_longer(-row, names_to = "col", values_to = "Count") %>% 
    mutate(col = as.integer(stringr::str_remove(col, "V")),
           time = sample_time)
})

df_sites$time <- factor(df_sites$time, levels = rev(unique(schedule$time)))

saveRDS(df_sites, file = "df_sites_1p.rds")
df_sites <- readRDS("df_sites_1p.rds")

p_sites <- ggplot(df_sites, aes(x = col, y = row, fill = Count)) +
  geom_tile() +
  scale_fill_viridis_c() +
  facet_wrap(.~time, labeller = labeller(time = my_labeller)) +
  labs(x = "EUR_1", y = "EUR_2")+theme(aspect.ratio=1)+ggtitle("Model X - recombination sites approach")

ggsave("figs/2d_tfs_1p_sites.png", plot = p_sites, width = 190, height = 143, units = "mm")


# windows approach --------------------------------------------------------

df_wind <- purrr::imap_dfr(time_list_windows, function(mat, sample_time) {
  mat[mat==0] <- NA
  as_tibble(mat) %>%
    mutate(row = row_number()-1) %>%
    pivot_longer(-row, names_to = "col", values_to = "Count") %>% 
    mutate(col = as.integer(stringr::str_remove(col, "V")),
           time = sample_time)
})

df_wind$time <- factor(df_wind$time, levels = rev(unique(schedule$time)))

saveRDS(df_wind, file = "df_wind_1p.rds")
df_wind <- readRDS("df_wind_1p.rds")

p_wind <- ggplot(df_wind, aes(x = col, y = row, fill = Count)) +
  geom_tile() +
  scale_fill_viridis_c() +
  facet_wrap(.~time, labeller = labeller(time = my_labeller)) +
  labs(x = "EUR_1", y = "EUR_2")+theme(aspect.ratio=1)+ggtitle("Model X - windows approach")


ggsave("figs/2d_tfs_1p_wind.png", plot = p_sites, width = 190, height = 143, units = "mm")


# unique ------------------------------------------------------------------

df_uni <- purrr::imap_dfr(time_list_unique, function(mat, sample_time) {
  mat[mat==0] <- NA
  as_tibble(mat) %>%
    mutate(row = row_number()-1) %>%
    pivot_longer(-row, names_to = "col", values_to = "Count") %>% 
    mutate(col = as.integer(stringr::str_remove(col, "V")),
           time = sample_time)
})

df_uni$time <- factor(df_uni$time, levels = rev(unique(schedule$time)))

saveRDS(df_uni, file = "df_uni_1p.rds")
df_uni <- readRDS("df_uni_1p.rds")

p_uni <- ggplot(df_uni, aes(x = col, y = row, fill = Count)) +
  geom_tile() +
  scale_fill_viridis_c() +
  facet_wrap(.~time, labeller = labeller(time = my_labeller)) +
  labs(x = "EUR_1", y = "EUR_2")+theme(aspect.ratio=1)+ggtitle("Model X - uniqueness approach")

ggsave("figs/2d_tfs_1p_uni.png", plot = p_uni, width = 190, height = 143, units = "mm")

plot_list <- list(uni = p_uni, wind= p_wind, sub=p_subset, sites=p_sites)
saveRDS(plot_list, "plots_1p.rds")
