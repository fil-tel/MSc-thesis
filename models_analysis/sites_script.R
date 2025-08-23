library(dplyr)
library(data.table)
library(readr)
library(parallel)
library(slendr)
init_env()

get_info_sites <- function(ts) {
  gt <- as.data.table(ts_genotypes(ts))

  afr <- grep("^AFR", colnames(gt), value = TRUE)
  neand <- grep("^NEA", colnames(gt), value = TRUE)

  # filter for African-Neanderthal fixed differences
  # (Africans = fixed ancestral, Neanderthals = fixed derived)
  ancestral <- 0
  afr_anc <- rowMeans(gt[, lapply(.SD, function(x) x == ancestral), .SDcols = afr]) == 1
  neand_der <- rowMeans(gt[, lapply(.SD, function(x) x != ancestral), .SDcols = neand]) == 1
  fixed_sites <- afr_anc & neand_der

  # select genotypes at fixed African-Neanderthal differences in Eurasians
  info_gt <- gt[fixed_sites, .SD, .SDcols = !c("pos", afr, neand)]

  info_gt
}

# bakcbine structure
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

# parameters
Ne_ooa_grid <- seq(500, 10000, 500)
seq_length <- 100e6
mc.cores <- detectCores()-1
n_runs <- 10


# model 1 A ---------------------------------------------------------------

dir_m1a <- "model1_A/"

# 10 runs
runs <- mclapply(1:n_runs, function(x) {
  lapply(Ne_ooa_grid, function(Ne_ooa) {
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

    prefix <- sprintf("model1_A_%s_rep%s", Ne_ooa, x)
    ts <- ts_read(file = paste0(dir_m1a, prefix, ".trees"), model = model)

    info_sites <- get_info_sites(ts)
    write_tsv(info_sites, file = paste0(dir_m1a, prefix, "_sites.tsv.gz"))

  })
}, mc.cores=mc.cores)


# model 1 B ---------------------------------------------------------------

dir_m1b <- "model1_B/"

runs <- mclapply(1:n_runs, function(x) {
  lapply(Ne_ooa_grid, function(Ne_ooa) {
    ooa <- population(
      "OOA",
      parent = afr,
      time = 70e3,
      N = Ne_ooa,
      remove = 39e3
    )
    # first vary Ne of the OOA receiving population (let's say 1k...20k Ne)
    # and compare the obtained TFS against the empirical IBDmix TFS

    ########################################
    # model 1
    #

    eur <- population("EUR",
                      parent = ooa,
                      time = 42e3,
                      N = 5000)

    # 2% introgressino into OOA
    gf <- list(
      gene_flow(
        from = nea,
        to = ooa,
        rate = 0.02,
        start = 49000,
        end = 45000,
        overlap = FALSE
      ),
      gene_flow(
        from = nea,
        to = eur,
        rate = 0.01,
        start = 42000 - 30,
        end = 40000,
        overlap = FALSE
      )
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
    #
    # plot_model(
    #   model,
    #   order = c("AFR", "ancestor", "EUR", "OOA", "ASN", "NEA"),
    #   proportions = TRUE,
    #   log = FALSE
    # )

    prefix <- sprintf("model1_B_%s_rep%s", Ne_ooa, x)
    ts <- ts_read(file = paste0(dir_m1b, prefix, ".trees"), model = model)

    info_sites <- get_info_sites(ts)
    write_tsv(info_sites, file = paste0(dir_m1b, prefix, "_sites.tsv.gz"))

  })
}, mc.cores=mc.cores)


# model 2 A ---------------------------------------------------------------

dir_m2a <- "model2_A/"

runs <- mclapply(1:n_runs, function(x) {
  lapply(Ne_ooa_grid, function(Ne_ooa) {
    ooa <- population(
      "OOA",
      parent = afr,
      time = 70e3,
      N = Ne_ooa,
      remove = 39e3
    )

    n_demes <- 9
    eur_demes <- lapply(seq_len(n_demes), function(i) {
      population(
        paste0("EUR_", i),
        parent = ooa,
        time = 42e3,
        N = 2000,
        remove = 29e3
      )
    })

    ana <- population(
      "ANA",
      parent = eur_demes[[1]],
      time = 30e3,
      N = 5000,
      remove = 3000 - 30
    )

    whg <- population(
      "WHG",
      parent = eur_demes[[4]],
      time = 30e3,
      N = 5000,
      remove = 3000 - 30
    )

    yam <- population(
      "YAM",
      parent = eur_demes[[7]],
      time = 30e3,
      N = 5000,
      remove = 3000 - 30
    )

    eur <- population("EUR",
                      parent = whg,
                      time = 3e3,
                      N = 15000)

    gf <- c(
      list(
        gene_flow(
          from = nea,
          to = ooa,
          rate = 0.03,
          start = 49000,
          end = 45000,
          overlap = FALSE
        )
      ),
      lapply(eur_demes[2:3], function(x) {
        gene_flow(
          from = x,
          to = ana,
          rate = 0.33,
          start = 30e3 - 30,
          end = 29e3 + 30,
          overlap = FALSE
        )
      }),
      lapply(eur_demes[5:6], function(x) {
        gene_flow(
          from = x,
          to = whg,
          rate = 0.33,
          start = 30e3 - 30,
          end = 29e3 + 30,
          overlap = FALSE
        )
      }),
      lapply(eur_demes[8:9], function(x) {
        gene_flow(
          from = x,
          to = yam,
          rate = 0.33,
          start = 30e3 - 30,
          end = 29e3 + 30,
          overlap = FALSE
        )
      }),
      list(
        gene_flow(
          from = ana,
          to = whg,
          rate = 0.5,
          start = 8500,
          end = 3000
        ),
        gene_flow(
          from = yam,
          to = whg,
          rate = 0.5,
          start = 5500,
          end = 3000
        )
      )
    )



    model <- compile_model(
      populations = append(list(anc, afr, nea, ooa, eur, yam, whg, ana), eur_demes),
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
    #
    # plot_model(
    #   model,
    #   order = c("AFR", "ancestor", "EUR", "OOA", "ASN", "NEA"),
    #   proportions = TRUE,
    #   log = FALSE
    # )

    prefix <- sprintf("model2_A_%s_rep%s", Ne_ooa, x)
    ts <- ts_read(file = paste0(dir_m2a, prefix, ".trees"), model = model)

    info_sites <- get_info_sites(ts)
    write_tsv(info_sites, file = paste0(dir_m2a, prefix, "_sites.tsv.gz"))

  })
}, mc.cores = mc.cores)



# model 2 B ---------------------------------------------------------------



dir_m2b <- "model2_B/"

runs <- mclapply(1:n_runs, function(x) {
  lapply(Ne_ooa_grid, function(Ne_ooa) {
    ooa <- population(
      "OOA",
      parent = afr,
      time = 70e3,
      N = Ne_ooa,
      remove = 39e3
    )

    n_demes <- 9
    eur_demes <- lapply(seq_len(n_demes), function(i) {
      population(
        paste0("EUR_", i),
        parent = ooa,
        time = 42e3,
        N = 2000,
        remove = 29e3
      )
    })

    ana <- population(
      "ANA",
      parent = eur_demes[[1]],
      time = 30e3,
      N = 5000,
      remove = 3000 - 30
    )

    whg <- population(
      "WHG",
      parent = eur_demes[[4]],
      time = 30e3,
      N = 5000,
      remove = 3000 - 30
    )

    yam <- population(
      "YAM",
      parent = eur_demes[[7]],
      time = 30e3,
      N = 5000,
      remove = 3000 - 30
    )

    eur <- population("EUR",
                      parent = whg,
                      time = 3e3,
                      N = 15000)

    gf <- c(
      list(
        gene_flow(
          from = nea,
          to = ooa,
          rate = 0.02,
          start = 49000,
          end = 45000,
          overlap = FALSE
        )
      ),
      lapply(eur_demes, function(x) {
        gene_flow(
          from = nea,
          to = x,
          rate = 0.01,
          start = 42e3 - 30,
          end = 40e3,
          overlap = FALSE
        )
      }),
      lapply(eur_demes[2:3], function(x) {
        gene_flow(
          from = x,
          to = ana,
          rate = 0.33,
          start = 30e3 - 30,
          end = 29e3 + 30,
          overlap = FALSE
        )
      }),
      lapply(eur_demes[5:6], function(x) {
        gene_flow(
          from = x,
          to = whg,
          rate = 0.33,
          start = 30e3 - 30,
          end = 29e3 + 30,
          overlap = FALSE
        )
      }),
      lapply(eur_demes[8:9], function(x) {
        gene_flow(
          from = x,
          to = yam,
          rate = 0.33,
          start = 30e3 - 30,
          end = 29e3 + 30,
          overlap = FALSE
        )
      }),
      list(
        gene_flow(
          from = ana,
          to = whg,
          rate = 0.5,
          start = 8500,
          end = 3000
        ),
        gene_flow(
          from = yam,
          to = whg,
          rate = 0.5,
          start = 5500,
          end = 3000
        )
      )
    )


    model <- compile_model(
      populations = append(list(anc, afr, nea, ooa, eur, yam, whg, ana), eur_demes),
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
    #
    # plot_model(
    #   model,
    #   order = c("AFR", "ancestor", "EUR", "OOA", "ASN", "NEA"),
    #   proportions = TRUE,
    #   log = FALSE
    # )

    prefix <- sprintf("model2_B_%s_rep%s", Ne_ooa, x)
    ts <- ts_read(file = paste0(dir_m2b, prefix, ".trees"), model = model)

    info_sites <- get_info_sites(ts)
    write_tsv(info_sites, file = paste0(dir_m2b, prefix, "_sites.tsv.gz"))

  })
}, mc.cores = mc.cores)


# model 3 a ---------------------------------------------------------------


dir_m3a <- "model3_A/"

runs <- mclapply(1:n_runs, function(x) {
  lapply(Ne_ooa_grid, function(Ne_ooa) {
    ooa <- population(
      "OOA",
      parent = afr,
      time = 70e3,
      N = Ne_ooa,
      remove = 39e3
    )

    n_demes <- 9
    eur_demes <- lapply(seq_len(n_demes), function(i) {
      population(
        paste0("EUR_", i),
        parent = ooa,
        time = 42e3,
        N = 2000,
        remove = 3e3 - 60
      )
    })

    eur <- population("EUR",
                      parent = eur_demes[[1]],
                      time = 3e3 + 30,
                      N = 15000)

    gf <- c(
      list(
        gene_flow(
          from = nea,
          to = ooa,
          rate = 0.03,
          start = 49000,
          end = 45000,
          overlap = FALSE
        )
      ),
      lapply(1:(length(eur_demes) - 1), function(x) {
        gene_flow(
          from = eur_demes[[x]],
          to = eur_demes[[x + 1]],
          rate = 1 / length(eur_demes),
          start = 30e3 - 30,
          end = 3e3 - 30,
          overlap = FALSE
        )
      }),
      lapply(length(eur_demes):2, function(x) {
        gene_flow(
          from = eur_demes[[x]],
          to = eur_demes[[x - 1]],
          rate = 1 / length(eur_demes),
          start = 30e3 - 30,
          end = 3e3 - 30,
          overlap = FALSE
        )
      }),
      lapply(eur_demes[2:length(eur_demes)], function(x) {
        gene_flow(
          from = x,
          to = eur,
          rate = 1 / length(eur_demes),
          start = 3e3,
          end = 3e3 - 30,
          overlap = FALSE
        )
      })
    )

    model <- compile_model(
      populations = append(list(anc, afr, nea, ooa, eur), eur_demes),
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

    prefix <- sprintf("model3_A_%s_rep%s", Ne_ooa, x)
    ts <- ts_read(file = paste0(dir_m3a, prefix, ".trees"), model = model)

    info_sites <- get_info_sites(ts)
    write_tsv(info_sites, file = paste0(dir_m3a, prefix, "_sites.tsv.gz"))

  })}, mc.cores = mc.cores)


# model 3 b ---------------------------------------------------------------


dir_m3b <- "model3_B/"

runs <- mclapply(1:n_runs, function(x) {
  lapply(Ne_ooa_grid, function(Ne_ooa) {
    ooa <- population(
      "OOA",
      parent = afr,
      time = 70e3,
      N = Ne_ooa,
      remove = 39e3
    )

    n_demes <- 9
    eur_demes <- lapply(seq_len(n_demes), function(i) {
      population(
        paste0("EUR_", i),
        parent = ooa,
        time = 42e3,
        N = 2000,
        remove = 3e3 - 60
      )
    })

    eur <- population("EUR",
                      parent = eur_demes[[1]],
                      time = 3e3 + 30,
                      N = 15000)

    gf <- c(
      list(
        gene_flow(
          from = nea,
          to = ooa,
          rate = 0.02,
          start = 49000,
          end = 45000,
          overlap = FALSE
        )
      ),
      lapply(eur_demes, function(x) {
        gene_flow(
          from = nea,
          to = x,
          rate = 0.01,
          start = 42e3 - 30,
          end = 40e3,
          overlap = FALSE
        )
      }),
      lapply(1:(length(eur_demes) - 1), function(x) {
        gene_flow(
          from = eur_demes[[x]],
          to = eur_demes[[x + 1]],
          rate = 1 / length(eur_demes),
          start = 30e3 - 30,
          end = 3e3 - 30,
          overlap = FALSE
        )
      }),
      lapply(length(eur_demes):2, function(x) {
        gene_flow(
          from = eur_demes[[x]],
          to = eur_demes[[x - 1]],
          rate = 1 / length(eur_demes),
          start = 30e3 - 30,
          end = 3e3 - 30,
          overlap = FALSE
        )
      }),
      lapply(eur_demes[2:length(eur_demes)], function(x) {
        gene_flow(
          from = x,
          to = eur,
          rate = 1 / length(eur_demes),
          start = 3e3,
          end = 3e3 - 30,
          overlap = FALSE
        )
      })
    )

    model <- compile_model(
      populations = append(list(anc, afr, nea, ooa, eur), eur_demes),
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

    prefix <- sprintf("model3_B_%s_rep%s", Ne_ooa, x)
    ts <- ts_read(file = paste0(dir_m3b, prefix, ".trees"), model = model)

    info_sites <- get_info_sites(ts)
    write_tsv(info_sites, file = paste0(dir_m3b, prefix, "_sites.tsv.gz"))

  })}, mc.cores = mc.cores)
