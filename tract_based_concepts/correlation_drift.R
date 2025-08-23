#########################
# In this script I explore the correlation between individuals sample
# indiviuals of tract and snp sharing

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
  remove = 59e3
)
pop1 <- population(
  "POP1",
  parent = ooa,
  time = 62e3,
  N = Ne
)

pop2 <- population(
  "POP2",
  parent = pop1,
  time = 55e3+30,
  N = Ne,
  remove = 50e3
)

pop3 <- population(
  "POP3",
  parent = pop2,
  time = 51e3,
  N = Ne
)

pop4 <- population(
  "POP4",
  parent = pop2,
  time = 51e3,
  N = Ne
)

# 3% introgressino into each EUR pop
gf <- list(gene_flow(
  from = nea,
  to = pop1,
  rate = 0.03,
  start = 55000,
  end = 51000,
  overlap = FALSE
), gene_flow(
  from = nea,
  to = pop2,
  rate = 0.03,
  start = 55000,
  end = 51000,
  overlap = FALSE
))

model <- compile_model(
  populations = list(anc, afr, nea, ooa, pop1, pop2, pop3, pop4),
  gene_flow = gf,
  generation_time = 30,
  serialize = FALSE,
  time_units = "years before present"
)

schedule <- schedule_sampling(
  model,
  times = c(50, 40, 30, 20, 10, 0) * 1e3,
  list(afr, 20),
  list(pop1, 20),
  list(pop3, 20),
  list(pop4, 20)
)

plot_model(
  model,
  order = c("AFR", "ancestor",    "OOA","POP1",  "POP3","POP2", "POP4", "NEA"),
  proportions = TRUE,
  file = "figs/model_correlation.png"
)
# define number of rows and column for the final matrx
# n_row <- unique(schedule$n)*2+1
# define list tp store results
# time_list_subset <- lapply(unique(schedule$time), function(x) matrix(0, nrow = n_row, ncol = n_row))
# time_list_sites <- lapply(unique(schedule$time), function(x) matrix(0, nrow = n_row, ncol = n_row))
# names(time_list_sites) <- names(time_list_subset) <- unique(schedule$time)

ts <- msprime(
  model,
  sequence_length = 10e6,
  recombination_rate = 1e-8,
  samples = schedule
) %>% ts_mutate(mutation_rate = 1e-8)

tracts <- ts_tracts(ts, census = 55000, squashed = TRUE)
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

# get pairwise f3
get_pairwise_f3 <- function(ts){
  samples <- ts_samples(ts)
}
  
samples <- ts_samples(ts) %>% filter(time==30000)

ts_f3(ts, A = c("POP1_2", "POP1_3"), B = "POP3_1", C = "AFR_1")

list <- group_split(samples, pop)

pops_df <- purrr::map2_dfc(list, "name", ~ select(.x, all_of(.y)))
colnames(pops_df) <- unique(samples$pop)

t <- pops_df %>% tidyr::expand(POP1, POP3)
t_h <- head(t)

# purrr::map(t, ts_f3(A = .x[[1]], B = , C 0))
ts_f3(ts, A = pops_df$POP3, B = pops_df$POP1, C = pops_df$AFR)

do.call(rbind, apply(t, MARGIN = 1, function(x) ts_f3(ts = ts, A = x[1], B =x[2], C = pops_df$AFR[1])))


