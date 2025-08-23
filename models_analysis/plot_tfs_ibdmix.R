# script to plot all the TFS

library(readr)
library(ggplot2)
library(cowplot)
library(BSgenome.Hsapiens.UCSC.hg19)
source("utils.R")

if(!dir.exists("figs_ibdmix")) dir.create("figs_ibdmix")

tfs_files <- list.files("tfs_df_ibdmix/", full.names = TRUE)
# get the df as a list
list_df <- lapply(tfs_files, read_tsv)
# join them together
def_df <- do.call(rbind, list_df) %>% rename("Ne"="Ne OOA")

def_df$approach <- factor(def_df$approach, levels = c("uni", "wind", "sub", "sites"))


p_sub <- def_df %>% filter(approach=="sub") %>% ggplot(aes(x = bin, y = mean, color = factor(model))) +
  facet_wrap(~`Ne OOA`, labeller = labeller(.rows = label_both), ncol = 4)+
  geom_line() +
  # geom_point() +
  # geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.2, color = NA) +
  labs(
    x = "Frequency class",
    y = "Proportion of tracts",
    color = "Model",
    fill = "Model",
    title = "Tract Frequency Spectrum",
    subtitle = "Sub-tracts approach"
  ) +theme(aspect.ratio = 1)
# theme_minimal()
p_sub

ggsave(plot = p_sub, filename = "figs_ibdmix/all_model_sub_ibdmix.png", width = 8, height = 10)


p_uni <- def_df %>% filter(approach=="uni") %>% ggplot(aes(x = bin, y = mean, color = factor(model))) +
  facet_wrap(~`Ne OOA`, labeller = labeller(.rows = label_both), ncol = 4)+
  geom_line() +
  # geom_point() +
  # geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.2, color = NA) +
  labs(
    x = "Frequency class",
    y = "Proportion of tracts",
    color = "Model",
    fill = "Model",
    title = "Tract Frequency Spectrum",
    subtitle = "Uniqueness approach"
  ) +theme(aspect.ratio = 1)
# theme_minimal()
p_uni

ggsave(plot = p_uni, filename = "figs_ibdmix/all_model_uni_ibdmix.png", width = 8, height = 10)


p_wind <- def_df %>% filter(approach=="wind") %>% ggplot(aes(x = bin, y = mean, color = factor(model))) +
  facet_wrap(~`Ne OOA`, labeller = labeller(.rows = label_both), ncol = 4)+
  geom_line() +
  # geom_point() +
  # geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.2, color = NA) +
  labs(
    x = "Frequency class",
    y = "Proportion of tracts",
    color = "Model",
    fill = "Model",
    title = "Tract Frequency Spectrum",
    subtitle = "Windows approach"
  ) +theme(aspect.ratio = 1)
# theme_minimal()
p_wind

ggsave(plot = p_wind, filename = "figs_ibdmix/all_model_wind_ibdmix.png", width = 8, height = 10)


p_sites <- def_df %>% filter(approach=="sites") %>% ggplot(aes(x = bin, y = mean, color = factor(model))) +
  facet_wrap(~`Ne OOA`, labeller = labeller(.rows = label_both), ncol = 4)+
  geom_line() +
  # geom_point() +
  # geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.2, color = NA) +
  labs(
    x = "Frequency class",
    y = "Proportion of tracts",
    color = "Model",
    fill = "Model",
    title = "Tract Frequency Spectrum",
    subtitle = "Recombination sites approach"
  ) +theme(aspect.ratio = 1)
# theme_minimal()
p_sites

ggsave(plot = p_sites, filename = "figs_ibdmix/all_model_sites_ibdmix.png", width = 8, height = 10)


# ribbon  -----------------------------------------------------------------

lapply(unique(def_df$model), function(mod) {
  p <- def_df %>% filter(model == mod, approach == "sub") %>% ggplot(aes(x = bin, y = mean)) +
    facet_wrap( ~ `Ne OOA`,
                labeller = labeller(.rows = label_both),
                ncol = 4) +
    geom_line() +
    # geom_point() +
    geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd),
                alpha = 0.2,
                color = NA) +
    labs(
      x = "Frequency class",
      y = "Proportion of tracts",
      color = "Model",
      fill = "Ne OOA",
      title = paste0("Tract Frequency Spectrum - ", mod),
      subtitle = "Sub-tracts approach"
    )+theme(aspect.ratio = 1)
  # theme_minimal()
  ggsave(plot = p, filename = sprintf("figs_ibdmix/ribbon_sub_%s_ibdmix.png", mod), width = 8, height = 10)
})


lapply(unique(def_df$model), function(mod) {
  p <- def_df %>% filter(model == mod, approach == "uni") %>% ggplot(aes(x = bin, y = mean)) +
    facet_wrap( ~ `Ne OOA`,
                labeller = labeller(.rows = label_both),
                ncol = 4) +
    geom_line() +
    # geom_point() +
    geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd),
                alpha = 0.2,
                color = NA) +
    labs(
      x = "Frequency class",
      y = "Proportion of tracts",
      color = "Model",
      fill = "Ne OOA",
      title = paste0("Tract Frequency Spectrum - ", mod),
      subtitle = "Uniqueness approach"
    )+theme(aspect.ratio = 1)
  # theme_minimal()
  ggsave(plot = p, filename = sprintf("figs_ibdmix/ribbon_uni_%s_ibdmix.png", mod), width = 8, height = 10)
})

lapply(unique(def_df$model), function(mod) {
  p <- def_df %>% filter(model == mod, approach == "sites") %>% ggplot(aes(x = bin, y = mean)) +
    facet_wrap( ~ `Ne OOA`,
                labeller = labeller(.rows = label_both),
                ncol = 4) +
    geom_line() +
    # geom_point() +
    geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd),
                alpha = 0.2,
                color = NA) +
    labs(
      x = "Frequency class",
      y = "Proportion of tracts",
      color = "Model",
      fill = "Ne OOA",
      title = paste0("Tract Frequency Spectrum - ", mod),
      subtitle = "Recombination sites approach"
    )+theme(aspect.ratio = 1)
  # theme_minimal()
  ggsave(plot = p, filename = sprintf("figs_ibdmix/ribbon_sites_%s_ibdmix.png", mod), width = 8, height = 10)
})


lapply(unique(def_df$model), function(mod) {
  p <- def_df %>% filter(model == mod, approach == "wind") %>% ggplot(aes(x = bin, y = mean)) +
    facet_wrap( ~ `Ne OOA`,
                labeller = labeller(.rows = label_both),
                ncol = 4) +
    geom_line() +
    # geom_point() +
    geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd),
                alpha = 0.2,
                color = NA) +
    labs(
      x = "Frequency class",
      y = "Proportion of tracts",
      color = "Model",
      fill = "Ne OOA",
      title = paste0("Tract Frequency Spectrum - ", mod),
      subtitle = "Windows approach"
    )+theme(aspect.ratio = 1)
  # theme_minimal()
  ggsave(plot = p, filename = sprintf("figs_ibdmix/ribbon_wind_%s_ibdmix.png", mod), width = 8, height = 10)
})



# all together ------------------------------------------------------------


app <- "sub"
p_models <- def_df %>% filter(approach == app) %>% ggplot(aes(
  x = bin,
  y = mean,
  color = factor(`Ne OOA`)
)) +
  facet_wrap( ~ model, ncol = 2) +
  geom_line() +
  # geom_point() +
  # geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.3, color = NA) +
  labs(
    x = "Frequency class",
    y = "Proportion of tracts",
    color = "Ne OOA",
    fill = "Model",
    title = "Tract Frequency Spectrum",
    subtitle = "Sub-tracts approach"
  ) + theme(aspect.ratio = 1)
# theme_minimal()
ggsave(plot = p_models, filename = sprintf("figs_ibdmix/full_model_%s_ibdmix.png", app), width = 6, height = 8)

app <- "sites"
p_models <- def_df %>% filter(approach == app) %>% ggplot(aes(
  x = bin,
  y = mean,
  color = factor(`Ne OOA`)
)) +
  facet_wrap( ~ model, ncol = 2) +
  geom_line() +
  # geom_point() +
  # geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.3, color = NA) +
  labs(
    x = "Frequency class",
    y = "Proportion of tracts",
    color = "Ne OOA",
    fill = "Model",
    title = "Tract Frequency Spectrum",
    subtitle = "Recombination sites approach"
  ) + theme(aspect.ratio = 1)
# theme_minimal()
ggsave(plot = p_models, filename = sprintf("figs_ibdmix/full_model_%s_ibdmix.png", app), width = 6, height = 8)

app <- "wind"
p_models <- def_df %>% filter(approach == app) %>% ggplot(aes(
  x = bin,
  y = mean,
  color = factor(`Ne OOA`)
)) +
  facet_wrap( ~ model, ncol = 2) +
  geom_line() +
  # geom_point() +
  # geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.3, color = NA) +
  labs(
    x = "Frequency class",
    y = "Proportion of tracts",
    color = "Ne OOA",
    fill = "Model",
    title = "Tract Frequency Spectrum",
    subtitle = "Windows approach"
  ) + theme(aspect.ratio = 1)
# theme_minimal()
ggsave(plot = p_models, filename = sprintf("figs_ibdmix/full_model_%s_ibdmix.png", app), width = 6, height = 8)

app <- "uni"
p_models <- def_df %>% filter(approach == app) %>% ggplot(aes(
  x = bin,
  y = mean,
  color = factor(`Ne OOA`)
)) +
  facet_wrap( ~ model, ncol = 2) +
  geom_line() +
  # geom_point() +
  # geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.3, color = NA) +
  labs(
    x = "Frequency class",
    y = "Proportion of tracts",
    color = "Ne OOA",
    fill = "Model",
    title = "Tract Frequency Spectrum",
    subtitle = "Uniqueness approach"
  ) + theme(aspect.ratio = 1)
# theme_minimal()
ggsave(plot = p_models, filename = sprintf("figs_ibdmix/full_model_%s_ibdmix.png", app), width = 6, height = 8)


# empirical data ----------------------------------------------------------


metadata <- read_metadata()
# tracts
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

n_sample <- 50

if(!file.exists("tfs_emp.tsv")) {
  tfs_emp_df <- mclapply(1:10, function(rep_i) {
    samples <- metadata$sampleId %>% sample(size = n_sample)
    tracts_samp_gr <- tracts_gr %>% filter(name %in% samples)
    bin_mat_subset <- get_bin_mat_subset_emp(tracts_samp_gr)
    bin_mat_windows <- get_bin_mat_windows_emp(tracts_samp_gr)
    bin_mat_uni <- get_bin_mat_unique(tracts_samp_gr)
    bin_mat_sites <- get_bin_mat_sites_emp(tracts_samp_gr)
    tfs_bins <- 1:n_sample
    tfs_subset <- compute_tfs(bin_mat_subset, n_sample = n_sample)
    sub_tbl <- tibble(
      rep_i = rep_i,
      bin = tfs_bins,
      count = tfs_subset,
      approach = "sub"
    )
    tfs_wind <- compute_tfs(bin_mat_windows, n_sample = n_sample)
    wind_tbl <- tibble(
      rep_i = rep_i,
      bin = tfs_bins,
      count = tfs_wind,
      approach = "wind"
    )
    tfs_uni <- compute_tfs(bin_mat_uni, n_sample = n_sample)
    uni_tbl <- tibble(
      rep_i = rep_i,
      bin = tfs_bins,
      count = tfs_uni,
      approach = "uni"
    )
    tfs_sites <- compute_tfs(bin_mat_sites, n_sample = n_sample)
    sites_tbl <- tibble(
      rep_i = rep_i,
      bin = tfs_bins,
      count = tfs_sites,
      approach = "sites"
    )
    print(rep_i)
    rbind(sub_tbl, sites_tbl, uni_tbl, wind_tbl)
  }, mc.cores = 3) %>% do.call(rbind, .) %>% group_by(bin, approach) %>% summarise(mean_count = mean(count), sd_count =
                                                                                     sd(count)) %>% ungroup()
  write_tsv(tfs_emp_df, file = "tfs_df_ibdmix/tfs_emp.tsv")
}

tfs_emp_df <- read_tsv("tfs_emp.tsv")
tfs_emp_df$approach <- factor(tfs_emp_df$approach, levels = c("uni", "wind", "sub", "sites"))

approach_labels <- c(
  "sites" = "Recombination breakpoints",
  "uni" = "Uniqueness approach",
  "wind" = "Windows approach",
  "sub" = "Sub-tracts approach"
)

p <- tfs_emp_df %>% ggplot(aes(x = bin, y = mean_count)) +
  geom_line() +
  facet_wrap(~approach, labeller = labeller(approach=approach_labels))+
  # geom_point() +
  geom_ribbon(aes(ymin = mean_count - sd_count, ymax = mean_count + sd_count),
              alpha = 0.2,
              color = NA) +
  labs(
    x = "Frequency class",
    y = "Proportion of tracts",
    color = "Model",
    fill = "Ne OOA",
  )+theme(aspect.ratio = 1)

p


ggsave(plot = p, filename = "figs_ibdmix/emp_data_only.png", width = 6, height = 6)



p_emp <- ggplot() +
  geom_line(
    data = def_df,
    aes(x = bin, y = mean, color = factor(`Ne OOA`))
  ) +
  geom_line(
    data = tfs_emp_df,
    aes(x = bin, y = mean_count, linetype = "Empirical data"),
    alpha = 0.7
  ) +
  facet_grid(model ~ approach, labeller = labeller(approach=approach_labels)) +
  scale_linetype_manual(values = c("Empirical data" = "dashed")) +
  labs(
    x = "Frequency class",
    y = "Proportion of tracts",
    color = "Ne OOA",
    fill = "Model",
    title = "Tract Frequency Spectrum",
    subtitle = "Models vs. Empirical Data",
    linetype = ""
  ) +
  theme(aspect.ratio = 1)
p_emp


ggsave(plot = p_emp, filename = "figs_ibdmix/emp_data.png", width = 10, height = 13)

p_emp <- ggplot() +
  geom_line(
    data = def_df,
    aes(x = bin, y = mean, color = factor(`Ne OOA`))
  ) +
  geom_line(
    data = tfs_emp_df,
    aes(x = bin, y = mean_count, linetype = "Empirical data"),
    alpha = 0.7
  ) +
  facet_grid(approach~model, labeller = labeller(approach=approach_labels), scales = "free_y") +
  scale_linetype_manual(values = c("Empirical data" = "dashed")) +
  labs(
    x = "Frequency class",
    y = "Proportion of tracts",
    color = "Ne OOA",
    fill = "Model",
    linetype = ""
  ) +
  theme(aspect.ratio = 1)
p_emp

ggsave(plot = p_emp, filename = "figs_ibdmix/models_vs_approach_hor_emp.png", width = 14, height = 11)

# only some of the ne OOA
p_emp <- ggplot() +
  geom_line(
    data = filter(def_df, `Ne OOA` %in% c(500, 1500, 2500, 3500, 4500, 5500)),
    aes(x = bin, y = mean, color = model)
  ) +
  geom_line(
    data = tfs_emp_df,
    aes(x = bin, y = mean_count, linetype = "Empirical data"),
    alpha = 0.7
  ) +
  facet_grid(approach~`Ne OOA`, labeller = labeller(approach=approach_labels), scales = "free_y") +
  scale_linetype_manual(values = c("Empirical data" = "dashed")) +
  labs(
    x = "Frequency class",
    y = "Proportion of tracts",
    color = "Model",
    fill = "Model",
    linetype = ""
  ) +
  theme(aspect.ratio = 1)
p_emp

ggsave(plot = p_emp, filename = "figs_ibdmix/models_vs_approach_hor_emp_ne.png", width = 14, height = 11)


# all together ------------------------------------------------------------


app <- "sub"
p_models <- def_df %>% filter(approach == app) %>% ggplot() +  geom_line(aes(
  x = bin,
  y = mean,
  color = factor(`Ne OOA`)
))+
  facet_wrap( ~ model, ncol = 2) +
  geom_line(
    data = filter(tfs_emp_df, approach==app),
    aes(x = bin, y = mean_count, linetype = "Empirical data"),
    alpha = 0.7
  ) +scale_linetype_manual(values = c("Empirical data" = "dashed")) +
  # geom_point() +
  # geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.3, color = NA) +
  labs(
    x = "Frequency class",
    y = "Proportion of tracts",
    color = "Ne OOA",
    fill = "Model",
    title = "Tract Frequency Spectrum",
    subtitle = "Sub-tracts approach", 
    linetype=""
  ) + theme(aspect.ratio = 1)
p_models
# theme_minimal()
ggsave(plot = p_models, filename = sprintf("figs_ibdmix/full_model_%s_ibdmix_emp.png", app), width = 6, height = 8)

app <- "sites"
p_models <- def_df %>% filter(approach == app) %>% ggplot() +  geom_line(aes(
  x = bin,
  y = mean,
  color = factor(`Ne OOA`)
))+
  facet_wrap( ~ model, ncol = 2) +
  geom_line(
    data = filter(tfs_emp_df, approach==app),
    aes(x = bin, y = mean_count, linetype = "Empirical data"),
    alpha = 0.7
  ) + scale_linetype_manual(values = c("Empirical data" = "dashed")) +
  labs(
    x = "Frequency class",
    y = "Proportion of tracts",
    color = "Ne OOA",
    fill = "Model",
    title = "Tract Frequency Spectrum",
    subtitle = "Recombination sites approach",
    linetype = ""
  ) + theme(aspect.ratio = 1)
p_models
ggsave(plot = p_models, filename = sprintf("figs_ibdmix/full_model_%s_ibdmix_emp.png", app), width = 6, height = 8)


app <- "wind"
p_models <- def_df %>% filter(approach == app) %>% ggplot() +  geom_line(aes(
  x = bin,
  y = mean,
  color = factor(`Ne OOA`)
))+
  facet_wrap( ~ model, ncol = 2) +
  geom_line(
    data = filter(tfs_emp_df, approach==app),
    aes(x = bin, y = mean_count, linetype = "Empirical data"),
    alpha = 0.7
  ) + scale_linetype_manual(values = c("Empirical data" = "dashed")) +
  labs(
    x = "Frequency class",
    y = "Proportion of tracts",
    color = "Ne OOA",
    fill = "Model",
    title = "Tract Frequency Spectrum",
    subtitle = "Windows approach",
    linetype = ""
  ) + theme(aspect.ratio = 1)
p_models
ggsave(plot = p_models, filename = sprintf("figs_ibdmix/full_model_%s_ibdmix_emp.png", app), width = 6, height = 8)

app <- "uni"
p_models <- def_df %>% filter(approach == app) %>% ggplot() +  geom_line(aes(
  x = bin,
  y = mean,
  color = factor(`Ne OOA`)
))+
  facet_wrap( ~ model, ncol = 2) +
  geom_line(
    data = filter(tfs_emp_df, approach==app),
    aes(x = bin, y = mean_count, linetype = "Empirical data"),
    alpha = 0.7
  ) + scale_linetype_manual(values = c("Empirical data" = "dashed")) +
  labs(
    x = "Frequency class",
    y = "Proportion of tracts",
    color = "Ne OOA",
    fill = "Model",
    title = "Tract Frequency Spectrum",
    subtitle = "Uniqueness approach",
    linetype = ""
  ) + theme(aspect.ratio = 1)
p_models
ggsave(plot = p_models, filename = sprintf("figs_ibdmix/full_model_%s_ibdmix_emp.png", app), width = 6, height = 8)

