library(data.table)
library(dplyr)
library(ggplot2)
library(readr)

if(!dir.exists("data_sfs")) dir.create("data_sfs")
if(!dir.exists("figs_sfs")) dir.create("figs_sfs")

# SFS analysis ------------------------------------------------------------

# Read table of Neanderthal-derived alleles in a given simulation and compute SFS
read_sfs <- function(model) {
  model_dir <- model

  sites_files <- list.files(model_dir, pattern = "*_sites.tsv.gz", full.names = TRUE)

  lapply(seq_along(sites_files), function(i) {
    f <- sites_files[i]

    cat(sprintf("Reading file %s [%s / %s]\n", f, i, length(sites_files)))

    model <- dirname(f)

    Ne_ooa <- as.integer(gsub("model._._(\\d+)_rep\\d+_sites.tsv.gz", "\\1", basename(f)))
    rep_i <- as.integer(gsub("model._._\\d+_rep(\\d+)_sites.tsv.gz", "\\1", basename(f)))

    df <- fread(f)

    sfs_counts <- table(rowSums(df))
    sfs_bins <- 1:100
    sfs <- rep(0, 100)
    if(!is.na(sfs_counts["0"])) sfs_counts <- sfs_counts[-1]
    sfs[as.numeric(names(sfs_counts))] <- sfs_counts
    # normalize to get proportion of tracts
    sfs <- sfs/sum(sfs)

    tibble(model = model, rep_i = rep_i, Ne_ooa = Ne_ooa, bin = sfs_bins, count = sfs)
  }) %>% do.call(rbind, .)
}

# Sanity check to compute the estimate of Neanderthal % in each EUR individual
# from the informative sites (we should get about 3% on average)
sim_sites <- fread("model1_A/model1_A_1000_rep1_sites.tsv.gz")
sim_props <- sim_sites[, sapply(.SD, mean)]
hist(sim_props, xlim = c(0, 0.1))

# Sanity check to visualize the average SFS for all Ne_ooa values (averaged
# across replicates), for a single model A here -- we should see the expected
# influence of drift

models <- c(paste0("model", 1:3, "_A"), paste0("model", 1:3, "_B"))
models_sfs <- do.call(rbind, lapply(models, read_sfs)) 

models_sfs_sd <- models_sfs %>%
  group_by(model, Ne_ooa, bin) %>%
  summarise(mean_count = mean(count), sd_count=sd(count)) %>% ungroup()

write_tsv(x = models_sfs_sd, file = "data_sfs/sfs_df.tsv")
models_sfs_sd <- read_tsv("data_sfs/sfs_df.tsv")

ne_labeller <- function(variable, value) {
  paste("Ne OOA:", value)
}

models_sfs_sd %>% ggplot(aes(bin, mean_count, color = factor(model))) +
  geom_line() + geom_ribbon(mapping = aes(
    ymin = mean_count - sd_count,
    ymax = mean_count + sd_count
  ), alpha = 0.2,
  color = NA) + facet_wrap( ~ Ne_ooa, ncol = 4, labeller = ne_labeller) + theme(aspect.ratio = 1) +
  labs(y = "Proportion", x = "Frequency class", color = "Ne OOA", fill="none")


models_sfs_sd %>% filter(model == "model2_A") %>% ggplot(aes(bin, mean_count, color = factor(Ne_ooa))) +
  geom_line() + geom_ribbon(mapping = aes(
    ymin = mean_count - sd_count,
    ymax = mean_count + sd_count
  ), alpha = 0.2,
  color = NA) + facet_wrap( ~ Ne_ooa, ncol = 4, labeller = ne_labeller) + theme(aspect.ratio = 1) +
  labs(y = "Proportion", x = "Frequency class", color = "Ne OOA", fill="none")


p <- models_sfs_sd %>% ggplot(aes(x = bin, y = mean_count, color = factor(model))) +
  facet_wrap(~Ne_ooa, labeller =  ne_labeller, ncol = 4)+
  geom_line() +
  # geom_point() +
  # geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.2, color = NA) +
  labs(
    x = "Frequency class",
    y = "Proportion of alleles",
    color = "Model",
    fill = "Model",
    title = "Site Frequency Spectrum of Introgressed Alleles"
  ) +theme(aspect.ratio = 1)
# theme_minimal()
p

ggsave(plot = p, filename = "figs_sfs//all_model_sfs.png", width = 8, height = 10)



lapply(unique(models_sfs_sd$model), function(mod) {
  p <- models_sfs_sd %>% filter(model == mod) %>% ggplot(aes(x = bin, y = mean_count)) +
    facet_wrap( ~ Ne_ooa,
                labeller = ne_labeller,
                ncol = 4) +
    geom_line() +
    # geom_point() +
    geom_ribbon(aes(ymin = mean_count - sd_count, ymax = mean_count + sd_count),
                alpha = 0.3,
                color = NA) +
    labs(
      x = "Frequency class",
      y = "Proportion of alleles",
      color = "Model",
      fill = "Ne OOA",
      title = paste0("Site Frequency Spectrum of Introgressed Alleles - ", mod),
    )+theme(aspect.ratio = 1)
  # theme_minimal()
  ggsave(plot = p, filename = sprintf("figs_sfs//ribbon_%s.png", mod), width = 8, height = 10)
})




p_models <- models_sfs_sd %>% ggplot(aes(
  x = bin,
  y = mean_count,
  color = factor(Ne_ooa)
)) +
  facet_wrap(~ model, ncol = 2) +
  geom_line() +
  # geom_point() +
  # geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.3, color = NA) +
  labs(
    x = "Frequency class",
    y = "Proportion of alleles",
    color = "Ne OOA",
    fill = "Model"
  ) + theme(aspect.ratio = 1)
p_models
# theme_minimal()
ggsave(
  plot = p_models,
  filename = "figs_sfs//full_models_sfs.png",
  width = 6,
  height = 8
)


# You should of course visualize this for all models, and present those figures for each



# Comparison with SFS and ancestry proportions in real data

# read equivalent of informative sites on real data (fixed African-vs-Vindija differences)
# (I will tell you more about the methods used to get these sites)
emp_sites <- fread("https://raw.githubusercontent.com/bodkan/ibdmix-dating/refs/heads/main/data/info_gt.tsv.gz")

# read metadata
metadata <- fread("https://raw.githubusercontent.com/bodkan/ibdmix-dating/refs/heads/main/data/neo.impute.1000g.sampleInfo_clusterInfo.txt")

# get all samples present in the info sites table
samples_all <- colnames(emp_sites[, .SD, .SDcols = !c("chrom", "pos", "ref", "alt", "AA")]) %>% gsub("_hap\\d$", "", .) %>% unique

# split samples into groups of interest
samples_1kg_afr <- intersect(samples_all, metadata[groupAge == "Modern" & grepl("Africa", region), sampleId])
samples_1kg_eur <- intersect(samples_all, metadata[groupAge == "Modern" & grepl("Europe", region), sampleId])
# unfortunately I messed up the generation of sites and only took "ancient" individuals
# as those from MesoNeo (the data set has much more, from other studies, but that's OK)
samples_emh <- intersect(samples_all, metadata[groupAge == "Ancient", sampleId])
samples_arch <- "neand" # we only have a single chromosome for Vindija, because it's fixed sites anyway

# compute proportions of Neanderthal informative sites in each set of samples
props_1kg_afr <- emp_sites[, sapply(.SD, mean), .SDcols = paste0(samples_1kg_afr, c("_hap1", "_hap2"))]
props_1kg_eur <- emp_sites[, sapply(.SD, mean), .SDcols = paste0(samples_1kg_eur, c("_hap1", "_hap2"))]
props_emh <- emp_sites[, sapply(.SD, mean), .SDcols = paste0(samples_emh, c("_hap1", "_hap2"))]
props_arch <- emp_sites[, sapply(.SD, mean), .SDcols = samples_arch]

# convert those proportions into data frames for plotting
props_1kg_afr_df <- tibble(name = names(props_1kg_afr), prop = props_1kg_afr, set = "present-day AFR")
props_1kg_eur_df <- tibble(name = names(props_1kg_eur), prop = props_1kg_eur, set = "present-day EUR")
props_emh_df <- tibble(name = names(props_emh), prop = props_emh, set = "EMH")
props_arch_df <- tibble(name = names(props_emh), prop = props_arch, set = "archaic")

# combine into a single data frame
props_df <- rbind(props_1kg_afr_df, props_1kg_eur_df, props_emh_df, props_arch_df)

props_df %>%
  filter(set != "archaic") %>%
  ggplot() +
  geom_density(aes(prop, fill = set), alpha = 0.8) +
  coord_cartesian(xlim = c(0, 0.1))

# Compute SFS for present-day modern humans+

sfs_emp_df <- lapply(1:10, function(rep_i) {
  samples_1kg_eur_hap <- samples_1kg_eur %>% sample(size = 50) %>% {c(paste0(., "_hap1"), paste0(., "_hap2"))}
  sites_1kg_eur <- emp_sites %>% select(any_of(samples_1kg_eur_hap))
  sfs <- rep(0, 100)
  sfs_counts <- table(rowSums(sites_1kg_eur))
  sfs_bins <- 1:100
  sfs <- rep(0, 100)
  if(!is.na(sfs_counts["0"])) sfs_counts <- sfs_counts[-1]
  sfs[as.numeric(names(sfs_counts))] <- sfs_counts
  # normalize to get proportion of tracts
  sfs <- sfs/sum(sfs)
  tibble(rep_i = rep_i, bin = sfs_bins, count = sfs)
}) %>% do.call(rbind, .) %>% group_by(bin) %>% summarise(mean_count = mean(count), sd_count=sd(count)) %>% ungroup()


pt <- sfs_emp_df %>% ggplot(aes(bin, mean_count)) +
  geom_line() + geom_ribbon(mapping = aes(
    ymin = mean_count - sd_count,
    ymax = mean_count + sd_count
  ), alpha = 0.3,
  color = NA) + theme(aspect.ratio = 1) +
  labs(y = "Proportion of alleles", x = "Frequency class", title = "Site Frequency Spectrum of Introgressed Alleles", subtitle = "Present-day Western Eurasians (n=50)")
pt

ggsave(plot = pt, filename = "figs_sfs//emp_sfs.png")

p_emp <- ggplot() + geom_line(data = models_sfs_sd, aes(
  x = bin,
  y = mean_count,
  color = factor(Ne_ooa)
)) +  geom_line(data = sfs_emp_df, aes(x = bin, y = mean_count, linetype="Empirical data"), alpha=0.7) +
  facet_wrap( ~ model, ncol = 2) +
  scale_linetype_manual(values = c("Empirical data" = "dashed"))+
  # geom_point() +
  # geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.3, color = NA) +
  labs(
    x = "Frequency class",
    y = "Proportion of alleles",
    color = "Ne OOA",
    fill = "Model",
    linetype = ""
  ) + theme(aspect.ratio = 1)
p_emp

ggsave(
  plot = p_emp,
  filename = "figs_sfs//emp_sfs_and_sim.png",
  width = 6,
  height = 8
)

