library(dplyr)
library(ggplot2)
library(readr)
library(parallel)
library(slendr)
init_env()

if(!dir.exists("heterozygosity")) dir.create("heterozygosity")

# Apparenlty defining it only once work the same,
# without having to iterate over all the models
Ne_ooa <- 3000

anc <- population("ancestor",
                  time = 700e3,
                  N = 15000,
                  remove = 640e3)
afr <- population("AFR",
                  parent = anc,
                  time = 650e3,
                  N = 15000)
nea <- population(
  "NEA",
  parent = anc,
  time = 650e3,
  N = 2000,
  remove = 39e3
)


ooa <- population(
  "OOA",
  parent = afr,
  time = 70e3,
  N = Ne_ooa,
  remove = 39e3
)

eur <- population("EUR",
                  parent = ooa,
                  time = 42e3,
                  N = 5000)
# 3% introgressino into OOA
gf <- gene_flow(
  from = nea,
  to = ooa,
  rate = 0.03,
  start = 49000,
  end = 45000,
  overlap = FALSE
)

model <- compile_model(
  populations = list(anc, afr, nea, ooa, eur),
  gene_flow = gf,
  generation_time = 30,
  serialize = FALSE,
  time_units = "years before present"
)

schedule <- rbind(
  schedule_sampling(model, times = 50e3, list(nea, 2)),
  schedule_sampling(model, times = 0, list(eur, 50)),
  schedule_sampling(model, times = 0, list(afr, 5))
)

suffices <- c(paste0(1:3, c("_A")), paste0(1:3, c("_B")))


if (!file.exists("heterozygosity/heterozygosity.tsv")) {
  het_df <- mclapply(suffices, function(suffix) {
    ts_files <- list.files(sprintf("model%s/", suffix),
                           pattern = "*.trees",
                           full.names = TRUE)
    lapply(ts_files, function(f) {
      model_name <- dirname(f)
      
      Ne_ooa <- as.integer(gsub("model._._(\\d+)_rep\\d+.trees", "\\1", basename(f)))
      rep_i <- as.integer(gsub("model._._\\d+_rep(\\d+).trees", "\\1", basename(f)))
      ts <- ts_read(file = f, model = model)
      sample_sets <- ts_samples(ts) %>%
        split(., .$pop) %>%
        lapply(function(pop)
          pop$name) %>% .[c("AFR", "EUR")] %>% unlist() %>% unname
      # print(rep_i)
      ts_diversity(ts, sample_sets = sample_sets) %>% mutate(
        pop = substr(set, 1, 3),
        model = model_name,
        rep = rep_i,
        Ne_ooa = Ne_ooa
      )
    }) %>% do.call(rbind, .)
  }, mc.cores = 3) %>% do.call(rbind, .)
  write_tsv(het_df, "heterozygosity/heterozygosity.tsv")
}else{
  het_df <- read_tsv("heterozygosity/heterozygosity.tsv")
}


het_df %>%
  ggplot() +
  geom_density(aes(diversity, fill = factor(Ne_ooa)), alpha = 0.8) + facet_grid(pop~model)
  coord_cartesian(xlim = c(0, 0.005))

  

# empirical diversity -----------------------------------------------------

emp_df <- readxl::read_xlsx("data/41586_2016_BFnature18964_MOESM205_ESM.xlsx")
emp_df <- emp_df %>% filter(Region %in% c("Africa", "WestEurasia"))
emp_het_df <- emp_df %>% select(Region, `samtools autosomal heterozygosity`) %>% mutate(diversity = as.numeric(`samtools autosomal heterozygosity`)) %>% 
  na.omit() %>% mutate(pop=ifelse(Region=="WestEurasia", "EUR", "AFR")) %>%  select(-one_of("samtools autosomal heterozygosity", "Region"))


# library(ggpattern)

p <- ggplot() +
  geom_density(
    data = het_df,
    aes(diversity, fill = factor(Ne_ooa)),
    alpha = 0.5,
    stat = "density"
  ) + facet_grid(model ~ pop) + geom_density(
    data = emp_het_df,
    mapping=aes(diversity, colour = "Empirical data"),
    alpha = 0.8,
    stat = "density",
    fill = "red"
  ) + scale_colour_manual(values = c("Empirical data" = NA))+
  coord_cartesian(xlim = c(0, 0.0013))+labs(x="Heterozygosity", y="Density", fill="Ne OOA", colour="")
p

ggsave(plot = p, "heterozygosity/het_plot.png", width = 10, height = 8)
