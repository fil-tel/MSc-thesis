# script to plot all the TFS

library(readr)
library(ggplot2)
library(cowplot)
library(dplyr)
if(!dir.exists("figs")) dir.create("figs")

tfs_files <- list.files("tfs_df/", full.names = TRUE)
# get the df as a list
list_df <- lapply(tfs_files, read_tsv)
# join them together
def_df <- do.call(rbind, list_df) %>% rename("Ne"="Ne OOA")

def_df$approach <- factor(def_df$approach, levels = c("uni", "wind", "sub", "sites"))


p_sub <- def_df %>% filter(approach=="sub") %>% ggplot(aes(x = bin, y = mean, color = factor(model))) +
  facet_wrap(~`Ne OOA`, labeller = labeller(.rows = label_both), ncol = 4)+
  geom_line() +
  # geom_point() +
  # geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.3, color = NA) +
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

ggsave(plot = p_sub, filename = "figs/all_model_sub.png", width = 8, height = 10)


p_uni <- def_df %>% filter(approach=="uni") %>% ggplot(aes(x = bin, y = mean, color = factor(model))) +
  facet_wrap(~`Ne OOA`, labeller = labeller(.rows = label_both), ncol = 4)+
  geom_line() +
  # geom_point() +
  # geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.3, color = NA) +
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

ggsave(plot = p_uni, filename = "figs/all_model_uni.png", width = 8, height = 10)


p_wind <- def_df %>% filter(approach=="wind") %>% ggplot(aes(x = bin, y = mean, color = factor(model))) +
  facet_wrap(~`Ne OOA`, labeller = labeller(.rows = label_both), ncol = 4)+
  geom_line() +
  # geom_point() +
  # geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.3, color = NA) +
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

ggsave(plot = p_wind, filename = "figs/all_model_wind.png", width = 8, height = 10)


p_sites <- def_df %>% filter(approach=="sites") %>% ggplot(aes(x = bin, y = mean, color = factor(model))) +
  facet_wrap(~`Ne OOA`, labeller = labeller(.rows = label_both), ncol = 4)+
  geom_line() +
  # geom_point() +
  # geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.3, color = NA) +
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

ggsave(plot = p_sites, filename = "figs/all_model_sites.png", width = 8, height = 10)


# ribbon  -----------------------------------------------------------------

lapply(unique(def_df$model), function(mod) {
  p <- def_df %>% filter(model == mod, approach == "sub") %>% ggplot(aes(x = bin, y = mean)) +
    facet_wrap( ~ `Ne OOA`,
                labeller = labeller(.rows = label_both),
                ncol = 4) +
    geom_line() +
    # geom_point() +
    geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd),
                alpha = 0.3,
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
  ggsave(plot = p, filename = sprintf("figs/ribbon_sub_%s.png", mod), width = 8, height = 10)
})


lapply(unique(def_df$model), function(mod) {
  p <- def_df %>% filter(model == mod, approach == "uni") %>% ggplot(aes(x = bin, y = mean)) +
    facet_wrap( ~ `Ne OOA`,
                labeller = labeller(.rows = label_both),
                ncol = 4) +
    geom_line() +
    # geom_point() +
    geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd),
                alpha = 0.3,
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
  ggsave(plot = p, filename = sprintf("figs/ribbon_uni_%s.png", mod), width = 8, height = 10)
})

lapply(unique(def_df$model), function(mod) {
  p <- def_df %>% filter(model == mod, approach == "sites") %>% ggplot(aes(x = bin, y = mean)) +
    facet_wrap( ~ `Ne OOA`,
                labeller = labeller(.rows = label_both),
                ncol = 4) +
    geom_line() +
    # geom_point() +
    geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd),
                alpha = 0.3,
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
  ggsave(plot = p, filename = sprintf("figs/ribbon_sites_%s.png", mod), width = 8, height = 10)
})


lapply(unique(def_df$model), function(mod) {
  p <- def_df %>% filter(model == mod, approach == "wind") %>% ggplot(aes(x = bin, y = mean)) +
    facet_wrap( ~ `Ne OOA`,
                labeller = labeller(.rows = label_both),
                ncol = 4) +
    geom_line() +
    # geom_point() +
    geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd),
                alpha = 0.3,
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
  ggsave(plot = p, filename = sprintf("figs/ribbon_wind_%s.png", mod), width = 8, height = 10)
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
ggsave(plot = p_models, filename = sprintf("figs/full_model_%s.png", app), width = 6, height = 8)

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
ggsave(plot = p_models, filename = sprintf("figs/full_model_%s.png", app), width = 6, height = 8)

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
ggsave(plot = p_models, filename = sprintf("figs/full_model_%s.png", app), width = 6, height = 8)

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
ggsave(plot = p_models, filename = sprintf("figs/full_model_%s.png", app), width = 6, height = 8)

approach_labels <- c(
  "sites" = "Recombination breakpoints",
  "uni" = "Uniqueness approach",
  "wind" = "Windows approach",
  "sub" = "Sub-tracts approach"
)



p <- ggplot() +
  geom_line(
    data = def_df,
    aes(x = bin, y = mean, color = factor(`Ne OOA`))
  ) +
    facet_grid(approach~model, labeller = labeller(approach=approach_labels), scales = "free_y") +
  labs(
    x = "Frequency class",
    y = "Proportion of tracts",
    color = "Ne OOA",
    fill = "Model"
  ) +
  theme(aspect.ratio = 1)
p


ggsave(plot = p, filename = "figs/models_vs_approach_hor.png", width = 14, height = 11)


# only some of the ne OOA
p <- ggplot() +
  geom_line(
    data = filter(def_df, `Ne OOA` %in% c(500, 1500, 2500, 3500, 4500, 5500)),
    aes(x = bin, y = mean, color = model)
  ) +
  facet_grid(approach~`Ne OOA`, labeller = labeller(approach=approach_labels), scales = "free_y") +
  labs(
    x = "Frequency class",
    y = "Proportion of tracts",
    color = "Model",
    fill = "Model"
  ) +
  theme(aspect.ratio = 1)
p

ggsave(plot = p, filename = "figs/models_vs_approach_hor_ne.png", width = 14, height = 11)
