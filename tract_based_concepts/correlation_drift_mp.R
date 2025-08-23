#########################
# In this script I explore the correlation between individuals sample
# indiviuals of tract and snp sharing

library(ggplotify)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(comperes)
source("utils.R")


# devtools::install_github("bodkan/slendr")

library(slendr)
init_env()

if(!dir.exists("figs")) dir.create("figs")
# set.seed(1234)

# just for saving plots, e..g 1 pulse, 2 pulse bla
mode <- "mp"

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
  N = 1000
)
ooa <- population(
  "OOA",
  parent = afr,
  time = 70e3,
  N = Ne,
  remove = 59e3
)
pop1 <- population(
  "POP1",
  parent = ooa,
  time = 62e3,
  N = Ne,
  remove = 29e3
)

pop1a <- population(
  "POP1a",
  parent = pop1,
  time = 30e3,
  N = Ne
)

pop1b <- population(
  "POP1b",
  parent = pop1,
  time = 30e3,
  N = Ne,
  remove = 14e3
)

pop1bx <- population(
  "POP1bx",
  parent = pop1b,
  time = 15e3,
  N = Ne
)
pop1by <- population(
  "POP1by",
  parent = pop1b,
  time = 15e3,
  N = Ne
)

pop2 <- population(
  "POP2",
  parent = ooa,
  time = 62e3,
  N = Ne,
  remove = 50e3
)

pop3 <- population(
  "POP3",
  parent = pop2,
  time = 51e3,
  N = Ne,
  remove = 19e3
)

pop3a <- population(
  "POP3a",
  parent = pop3,
  time = 20e3,
  N = Ne
)

pop3b <- population(
  "POP3b",
  parent = pop3,
  time = 20e3,
  N = Ne
)

pop4 <- population(
  "POP4",
  parent = pop2,
  time = 51e3,
  N = Ne,
  remove = 11e3
)

pop4a <- population(
  "POP4a",
  parent = pop4,
  time = 12e3,
  N = Ne
)

pop4b <- population(
  "POP4b",
  parent = pop4,
  time = 12e3,
  N = Ne
)

# 3% introgressino into each EUR pop
gf <- list(gene_flow(
  from = nea,
  to = pop1a,
  rate = 0.03,
  start = 12000-30,
  end = 8000,
  overlap = FALSE
), gene_flow(
  from = nea,
  to = pop1bx,
  rate = 0.03,
  start = 12000-30,
  end = 8000,
  overlap = FALSE
),
gene_flow(
  from = nea,
  to = pop1by,
  rate = 0.03,
  start = 12000-30,
  end = 8000,
  overlap = FALSE
),gene_flow(
  from = nea,
  to = pop3a,
  rate = 0.03,
  start = 12000-30,
  end = 8000,
  overlap = FALSE
),gene_flow(
  from = nea,
  to = pop3b,
  rate = 0.03,
  start = 12000-30,
  end = 8000,
  overlap = FALSE
), gene_flow(
  from = nea,
  to = pop4a,
  rate = 0.03,
  start = 12000-30,
  end = 8000,
  overlap = FALSE
), gene_flow(
  from = nea,
  to = pop4b,
  rate = 0.03,
  start = 12000-30,
  end = 8000,
  overlap = FALSE
)
)

# gf <- list(gene_flow(
#   from = nea,
#   to = ooa,
#   rate = 0.03,
#   start = 66000,
#   end = 61000,
#   overlap = FALSE
# ))

model <- compile_model(
  populations = list(anc, afr, nea, ooa,
                     pop1, pop1a, pop1b, pop1bx, pop1by,
                     pop2, pop3, pop3a, pop3b, pop4, pop4a, pop4b),
  gene_flow = gf,
  generation_time = 30,
  serialize = FALSE,
  time_units = "years before present"
)

# schedule <- schedule_sampling(
#   model,
#   times = c(50, 40, 30, 20, 10, 0) * 1e3,
#   list(afr, 20),
#   list(pop1, 20),
#   list(pop3, 20),
#   list(pop4, 20)
# )

schedule <- schedule_sampling(
  model,
  times = 0,
  list(afr, 5),
  list(pop1, 5),
  list(pop1a, 5),
  list(pop1b, 5),
  list(pop1bx, 5),
  list(pop1by, 5),
  list(pop2, 5),
  list(pop3, 5),
  list(pop3a, 5),
  list(pop3b, 5),
  list(pop4, 5),
  list(pop4a, 5),
  list(pop4b, 5)
)

p_model <- plot_model(
  model,
  order = c("AFR", "ancestor",  "POP1a", "POP1",  "POP1bx", "POP1b", "POP1by",
            "OOA","POP3a", "POP3",  "POP3b", "POP2", "POP4a", "POP4", "POP4b", "NEA"),
  proportions = TRUE,
  samples = schedule
  # file = "figs/model_correlation.png"
)+scale_y_break(c(70e3, 6.3e5), scales = c(0.1,1))+theme(axis.text.y.right =  element_blank(),
                                                         axis.ticks.y.right =  element_blank(),
                                                         axis.line.y.right =  element_blank())

p_model

# define number of rows and column for the final matrx
# n_row <- unique(schedule$n)*2+1
# define list tp store results
# time_list_subset <- lapply(unique(schedule$time), function(x) matrix(0, nrow = n_row, ncol = n_row))
# time_list_sites <- lapply(unique(schedule$time), function(x) matrix(0, nrow = n_row, ncol = n_row))
# names(time_list_sites) <- names(time_list_subset) <- unique(schedule$time)

ts <- msprime(
  model,
  sequence_length = 50e6,
  recombination_rate = 1e-8,
  samples = schedule
) %>% ts_mutate(mutation_rate = 1e-8)

ts_samples(ts) %>% group_by(pop) %>% tally()

pairs_df <- ts_samples(ts) %>%
  filter(pop != "AFR") %>%
  { expand.grid(x = .$name, y = .$name) } %>%
  filter(x != y) %>%
  as_tibble

# f3 computation ----------------------------------------------------------

pairs_df <-
  pairs_df %>%
  rowwise() %>%
  mutate(f3 = ts_f3(ts, A = list(paste0("AFR_", 1:5)), B = as.character(x), C = as.character(y))$f3)

p_f3 <- pairs_df %>% ggplot(aes(x, y)) + geom_tile(aes(fill = f3))

cowplot::plot_grid(p_model, p_f3)

# save <- pairs_df
# pairs_df <- save

# mat_f3[lower.tri(mat_f3)] <- NA

# tracts ------------------------------------------------------------------


tracts <- ts_tracts(ts, census = 12000-30, squashed = TRUE)
# add chromosome column
tracts$chrom <- "chr1"

# tracts <- tracts %>% mutate(name = paste0(name, "_hap_", haplotype))

tracts_gr <- makeGRangesFromDataFrame(
  tracts,
  keep.extra.columns = TRUE,
  seqnames.field = "chrom",
  start.field = "left",
  end.field = "right"
)

bin_mat_sub <- get_bin_mat_subset(tracts_gr)
bin_mat_wind <- get_bin_mat_windows(tracts_gr, len_chr = 50e6)
bin_mat_uni <- get_bin_mat_unique(tracts_gr)
bin_mat_sites <- get_bin_mat_sites(tracts_gr)

get_df_2_plot <- function(bin_mat, pairs_df){
  pairs_df <- pairs_df %>% rowwise() %>% mutate(cor = cor(bin_mat[,x], bin_mat[,y])) %>% ungroup()
  # pairs_df <- pairs_df %>% rowwise() %>% mutate(cor = cor(bin_mat_wind[,x], bin_mat_wind[,y])) %>% ungroup()
  # pairs_df <- pairs_df %>% rowwise() %>% mutate(cor = cor(bin_mat_uni[,x], bin_mat_uni[,y])) %>% ungroup()
  pairs_df <- pairs_df %>% mutate(f3=(f3-min(f3))/(max(f3)-min(f3)), cor=(cor-min(cor)/(max(cor)-min(cor))))
  # get matrix for f3
  mat_f3 <- long_to_mat(pairs_df, "x", "y", "f3")
  # get matrix for tracts corrr
  mat_cor <- long_to_mat(pairs_df, "x", "y", "cor")
  # join them, uooer is f3, lower tract corr
  mat_f3[lower.tri(mat_f3)] <- mat_cor[lower.tri(mat_cor)]  
  mat_to_long(mat_f3, "x", "y", "value")
}

get_df_2_plot <- function(bin_mat, pairs_df){
  # bin_mat[bin_mat==0] <- -1  
  pairs_df <- pairs_df %>% rowwise() %>% mutate(cor = cor(bin_mat[,x], bin_mat[,y])) %>% ungroup()
  # pairs_df <- pairs_df %>% rowwise() %>% mutate(cor = 1-e1071::hamming.distance(bin_mat[,x], bin_mat[,y])/nrow(bin_mat)) %>% ungroup()
  # pairs_df <- pairs_df %>% rowwise() %>% mutate(cor = 1-c(dist(t(bin_mat[,c(x, y)]), method="binary"))) %>% ungroup()
  # pairs_df <- pairs_df %>% rowwise() %>% mutate(cor = cor(bin_mat_wind[,x], bin_mat_wind[,y])) %>% ungroup()
  # pairs_df <- pairs_df %>% rowwise() %>% mutate(cor = cor(bin_mat_uni[,x], bin_mat_uni[,y])) %>% ungroup()
  # pairs_df <- pairs_df %>% mutate(f3=(f3-min(f3))/(max(f3)-min(f3)), cor=(cor-min(cor)/(max(cor)-min(cor))))
  # get matrix for f3
  mat_f3 <- long_to_mat(pairs_df, "x", "y", "f3")
  # get matrix for tracts corrr
  mat_cor <- long_to_mat(pairs_df, "x", "y", "cor")
  # join them, uooer is f3, lower tract corr
  mat_f3[lower.tri(mat_f3)] <- mat_cor[lower.tri(mat_cor)]  
  mat_long <- mat_to_long(mat_f3, "x", "y", "value")
  levels_order <- sort(unique(c(mat_long$x, mat_long$y)))
  mat_long$x <- factor(mat_long$x, levels = levels_order)
  mat_long$y <- factor(mat_long$y, levels = levels_order)
  mat_long$triangle <- with(mat_long, ifelse(
    as.integer(x) < as.integer(y), "upper", "lower"
  ))
  mat_long
}

# get_df_2_plot <- function(bin_mat, pairs_df){
#   bin_mat[bin_mat==0] <- -1  
#   pairs_df <- pairs_df %>% rowwise() %>% mutate(cor = (bin_mat[,x]%*%bin_mat[,y])/nrow(bin_mat)) %>% ungroup()
#   # pairs_df <- pairs_df %>% rowwise() %>% mutate(cor = 1-e1071::hamming.distance(bin_mat[,x], bin_mat[,y])/nrow(bin_mat)) %>% ungroup()
#   # pairs_df <- pairs_df %>% rowwise() %>% mutate(cor = 1-c(dist(t(bin_mat[,c(x, y)]), method="binary"))) %>% ungroup()
#   # pairs_df <- pairs_df %>% rowwise() %>% mutate(cor = cor(bin_mat_wind[,x], bin_mat_wind[,y])) %>% ungroup()
#   # pairs_df <- pairs_df %>% rowwise() %>% mutate(cor = cor(bin_mat_uni[,x], bin_mat_uni[,y])) %>% ungroup()
#   # pairs_df <- pairs_df %>% mutate(f3=(f3-min(f3))/(max(f3)-min(f3)), cor=(cor-min(cor)/(max(cor)-min(cor))))
#   # get matrix for f3
#   mat_f3 <- long_to_mat(pairs_df, "x", "y", "f3")
#   # get matrix for tracts corrr
#   mat_cor <- long_to_mat(pairs_df, "x", "y", "cor")
#   # join them, uooer is f3, lower tract corr
#   mat_f3[lower.tri(mat_f3)] <- mat_cor[lower.tri(mat_cor)]  
#   mat_long <- mat_to_long(mat_f3, "x", "y", "value")
#   levels_order <- sort(unique(c(mat_long$x, mat_long$y)))
#   mat_long$x <- factor(mat_long$x, levels = levels_order)
#   mat_long$y <- factor(mat_long$y, levels = levels_order)
#   mat_long$triangle <- with(mat_long, ifelse(
#     as.integer(x) < as.integer(y), "upper", "lower"
#   ))
#   mat_long
# }



# get_df_2_plot <- function(bin_mat, pairs_df){
#   bin_mat[bin_mat==0] <- -1  
#   pairs_df <- pairs_df %>% rowwise() %>% mutate(cor = (bin_mat[,x]%*%bin_mat[,y])/nrow(bin_mat)) %>% ungroup()
#   # pairs_df <- pairs_df %>% rowwise() %>% mutate(cor = 1-e1071::hamming.distance(bin_mat[,x], bin_mat[,y])/nrow(bin_mat)) %>% ungroup()
#   # pairs_df <- pairs_df %>% rowwise() %>% mutate(cor = 1-c(dist(t(bin_mat[,c(x, y)]), method="binary"))) %>% ungroup()
#   # pairs_df <- pairs_df %>% rowwise() %>% mutate(cor = cor(bin_mat_wind[,x], bin_mat_wind[,y])) %>% ungroup()
#   # pairs_df <- pairs_df %>% rowwise() %>% mutate(cor = cor(bin_mat_uni[,x], bin_mat_uni[,y])) %>% ungroup()
#   # pairs_df <- pairs_df %>% mutate(f3=(f3-min(f3))/(max(f3)-min(f3)), cor=(cor-min(cor)/(max(cor)-min(cor))))
#   # get matrix for f3
#   mat_f3 <- long_to_mat(pairs_df, "x", "y", "f3")
#   # get matrix for tracts corrr
#   mat_cor <- long_to_mat(pairs_df, "x", "y", "cor")
#   # join them, uooer is f3, lower tract corr
#   mat_f3[lower.tri(mat_f3)] <- mat_cor[lower.tri(mat_cor)]  
#   mat_long <- mat_to_long(mat_f3, "x", "y", "value")
#   levels_order <- sort(unique(c(mat_long$x, mat_long$y)))
#   mat_long$x <- factor(mat_long$x, levels = levels_order)
#   mat_long$y <- factor(mat_long$y, levels = levels_order)
#   mat_long$triangle <- with(mat_long, ifelse(
#     as.integer(x) < as.integer(y), "upper", "lower"
#   ))
#   mat_long
# }

df_sub <- get_df_2_plot(bin_mat_sub, pairs_df)
df_wind <- get_df_2_plot(bin_mat_wind, pairs_df)
df_uni <- get_df_2_plot(bin_mat_uni, pairs_df)
df_sites <- get_df_2_plot(bin_mat_sites, pairs_df)

saveRDS(df_sub, file = "df_sub_mp.rds")
plot_corr <- function(df){
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_tile(data = subset(df, triangle == "lower"),
              aes(fill = value)) +
    scale_fill_gradient2(low = "lightblue", high = "blue", name = "Tracts sharing") +
    new_scale_fill() +
    geom_tile(data = subset(df, triangle == "upper"),
              aes(fill = value)) +
    scale_fill_gradient2(low = "blue", high = "red", name = "F3 statistics") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio = 1, axis.title.x = element_blank(), axis.title.y = element_blank())+
    labs(title = "Correlation between pairwise F3 statistic and introgressed tracts sharing")
  p
}


# mat_f3 <- (mat_f3-mean(mat_f3, na.rm=TRUE))/sd(mat_f3, na.rm = TRUE)
# bleurgh <- mat_to_long(mat_f3, "x", "y", "value")
# p_f3_sub <- df_sub %>% ggplot(aes(x, y)) + geom_tile(aes(fill = value))+scale_fill_viridis_c()+theme(aspect.ratio = 1, axis.text.x = element_text(angle = 45, hjust = 1))+
#   ggtitle("Sub-tracts approach")
p_f3_sub <- df_sub %>% plot_corr() + labs(subtitle = "Sub-tracts approach")
p_sub <- cowplot::plot_grid(p_model, p_f3_sub, nrow = 1,  labels = "AUTO")
p_sub
ggsave(plot = p_sub, filename = sprintf("figs/corr_sub_%s.png", mode), width = 14, height = 10)

# p_f3_wind <- df_wind %>% ggplot(aes(x, y)) + geom_tile(aes(fill = value))+scale_fill_viridis_c()+theme(aspect.ratio = 1, axis.text.x = element_text(angle = 45, hjust = 1))+
#   ggtitle("Windows approach")
p_f3_wind <- df_wind %>% plot_corr() + labs(subtitle = "Windows approach")
p_wind <- cowplot::plot_grid(p_model, p_f3_wind, nrow = 1, labels = "AUTO")
ggsave(plot = p_wind, filename = sprintf("figs/corr_wind_%s.png", mode), width = 14, height = 10)

# p_f3_sites <- df_sites %>% ggplot(aes(x, y)) + geom_tile(aes(fill = value))+scale_fill_viridis_c()+theme(aspect.ratio = 1, axis.text.x = element_text(angle = 45, hjust = 1))+
#   ggtitle("Recombination sites approach")
p_f3_sites <- df_sites %>% plot_corr() + labs(subtitle = "Recombination sites approach")
p_sites <- cowplot::plot_grid(p_model, p_f3_sites, nrow = 1, labels = "AUTO")
ggsave(plot = p_sites, filename = sprintf("figs/corr_sites_%s.png", mode), width = 14, height = 10)

# p_f3_uni <- df_uni %>% ggplot(aes(x, y)) + geom_tile(aes(fill = value))+scale_fill_viridis_c()+theme(aspect.ratio = 1, axis.text.x = element_text(angle = 45, hjust = 1))+
#   ggtitle("Uniqueness approach")
p_f3_uni <- df_uni %>% plot_corr() + labs(subtitle = "Uniqueness approach")
p_uni <- cowplot::plot_grid(p_model, p_f3_uni, nrow = 1, labels = "AUTO")
ggsave(plot = p_uni, filename = sprintf("figs/corr_uni_%s.png", mode), width = 14, height = 10)


# grid_plot <- cowplot::plot_grid(p_f3_uni, p_f3_wind, p_f3_sub, p_f3_sites, labels = LETTERS[2:5])
# model_and_grid <- cowplot::plot_grid(p_model, grid_plot, labels = "A", ncol = 1)

df_sub$mode <- "sub"
df_sites$mode <- "sites"
df_uni$mode <- "uni"
df_wind$mode <- "wind"

full_df <- rbind(df_sub, df_sites, df_uni, df_wind)
full_df$mode <- factor(full_df$mode, levels = c("uni", "wind", "sub", "sites"))

label_vec <- c(
  sub = "Sub-tracts approach",
  sites = "Recombination sites approach",
  uni = "Uniqueness approach",
  wind = "Windows appraoch"
)


p <- full_df %>% ggplot(aes(x, y)) +   geom_tile(data = subset(full_df, triangle == "lower"), aes(fill = value)) +
  scale_fill_gradient2(low = "lightblue",
                       high = "blue",
                       name = "Tracts sharing") +
  new_scale_fill() +
  geom_tile(data = subset(full_df, triangle == "upper"), aes(fill = value)) +
  scale_fill_gradient2(low = "blue",
                       high = "red",
                       name = "F3 statistics") +facet_wrap(~ mode, labeller = labeller(mode = label_vec), ncol =
                                                             2) + 
  theme(
    aspect.ratio = 1,
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0
    ),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

p
# last_p <- cowplot::plot_grid(p_model, p, labels = "AUTO", ncol = 2)
last_p <- aplot::plot_list(p_model, p, labels = c("A", "B"), ncol = 2)

ggsave(plot = last_p, filename = sprintf("figs/grid_plot_%s.png", mode), width = 15, height = 10)

