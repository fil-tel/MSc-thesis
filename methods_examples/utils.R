library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(parallel)

library(GenomicRanges)
library(plyranges)


# i/o functions -------------------------------------------------------------------------------

read_metadata <- function() {
  raw_info <- read_tsv("data/neo.impute.1000g.sampleInfo_clusterInfo.txt")

  info <-
    raw_info %>%
    filter(region %in% c("SouthernEurope", "WesternEurope", "WesternAsia", "NorthernEurope", "CentralEasternEurope"),
           country != "Greenland") %>%
    filter(groupAge != "Archaic") %>%
    mutate(ageAverage = ifelse(groupAge == "Modern", 0, ageAverage)) %>%
    mutate(coverage = ifelse(groupAge == "Modern", Inf, coverage))

  info
}

read_tracts <- function(set = NULL) {
  metadata <- read_metadata()
  
  if(!is.null(set)){
    info <- filter(metadata, groupAge == set) 
  }

  raw_tracts <- read_tsv(here::here("data/Vindija33.19_raw_eurasian_wModern"))

  tracts <- raw_tracts %>%
    mutate(chrom = paste0("chr", chrom)) %>%
    filter(ID %in% unique(info$sampleId)) %>%
    dplyr::select(ID, chrom, start, end) %>%
    dplyr::mutate(length = end - start, set = set)

  tracts
}

# archaic deserts -----------------------------------------------------------------------------

generate_windows <- function(gaps_gr, window_size, step_size) {
  autosomes_gr <- GenomeInfoDb::getChromInfoFromUCSC("hg19") %>%
    filter(grepl("chr\\d+$", chrom)) %>%
    mutate(start = 1) %>%
    dplyr::rename(end = size) %>%
    makeGRangesFromDataFrame()
  seqinfo(autosomes_gr) <- seqinfo(tracts_gr)

  windows_grl <- slidingWindows(autosomes_gr, width = window_size, step = step_size)
  for (i in seq_along(windows_grl)) {
    windows_grl[[i]]$midpoint <- (start(windows_grl[[i]]) + end(windows_grl[[i]])) / 2
    windows_grl[[i]]$gap <- FALSE
    to_remove <- queryHits(findOverlaps(windows_grl[[i]], gaps_gr))
    if (length(to_remove) > 0)
      windows_grl[[i]][to_remove]$gap <- TRUE
  }

  unlist(windows_grl)
}

generate_windows_no_gap <- function(window_size, step_size) {
  autosomes_gr <- GenomeInfoDb::getChromInfoFromUCSC("hg19") %>%
    filter(grepl("chr\\d+$", chrom)) %>%
    mutate(start = 1) %>%
    dplyr::rename(end = size) %>%
    makeGRangesFromDataFrame()
  seqinfo(autosomes_gr) <- seqinfo(tracts_gr)
  
  windows_grl <- slidingWindows(autosomes_gr, width = window_size, step = step_size)
  for (i in seq_along(windows_grl)) {
    windows_grl[[i]]$midpoint <- (start(windows_grl[[i]]) + end(windows_grl[[i]])) / 2
  }
  
  unlist(windows_grl)
}

generate_windows_sim <- function(n_chr=1, len_chr=100e6, window_size, step_size) {
  
  seq_info <- Seqinfo(paste0("chr", 1:n_chr), seqlengths = len_chr, isCircular = FALSE, genome = "sim")
  gr <- GRanges(paste0("chr", 1:n_chr), IRanges(rep(1, n_chr), len_chr), seqinfo = seq_info)
  
  windows_grl <- slidingWindows(gr, width = window_size, step = step_size)
  for (i in seq_along(windows_grl)) {
    windows_grl[[i]]$midpoint <- (start(windows_grl[[i]]) + end(windows_grl[[i]])) / 2
  }
  
  unlist(windows_grl)
}

compute_ancestry <- function(tracts_gr, windows_gr, keep_gaps = FALSE) {
  autosomes <- getChromInfoFromUCSC("hg19") %>%
    filter(grepl("chr\\d+$", chrom)) %>%
    mutate(start = 1) %>%
    dplyr::rename(end = size) %>%
    makeGRangesFromDataFrame()

  # first compute coverage...
  cov <- coverage(tracts_gr)
  # ... then convert that to proportions
  cov <- lapply(cov, function(x) x / length(unique(tracts_gr$sampleId)))

  ancestry_list <- mclapply(as.character(unique(seqnames(windows_gr))),
                            function(chrom) {
    chrom_coverage <- cov[[chrom]]
    chrom_gr <- windows_gr[seqnames(windows_gr) == chrom]

    # count overlaps between windows and tracts
    average_coverage_per_window <- sapply(
      seq_along(chrom_gr),
      function(i) {
        start_idx <- start(chrom_gr[i])
        end_idx <- end(chrom_gr[i])
        mean(chrom_coverage[start_idx:end_idx])
    })

    mcols(chrom_gr)$coverage <- average_coverage_per_window
    if (!keep_gaps)
      mcols(chrom_gr)$coverage[mcols(chrom_gr)$gap] <- NA

    chrom_gr
  }, mc.cores = detectCores()-2)

  ancestry_grl <- GRangesList(ancestry_list)


  unlist(ancestry_grl)
}

plot_desert_ancestry <- function(ancestry_gr, deserts_gr, chrom, full = FALSE) {
  desert_df <- deserts_gr %>% filter(seqnames == chrom) %>% as_tibble

  ancestry_df <- as_tibble(ancestry_gr) %>%
    filter(seqnames == chrom) %>%
    select(chrom = seqnames, start, end, ancient, modern, gap, midpoint)

  if (!full) {
    ancestry_df <- ancestry_df %>%
      filter(start >= (desert_df$start * 0.9) & end <= (desert_df$end * 1.1))
  }

  p <- ancestry_df %>%
    {
      ggplot(data = .) +
        geom_rect(data = desert_df, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = "red", alpha = 0.05) +

        geom_line(aes(midpoint, ancient), color = "blue") +

        geom_line(aes(midpoint, modern), color = "orange") +

        geom_vline(data = desert_df, aes(xintercept = start, linetype = "desert boundary"), color = "red") +
        geom_vline(data = desert_df, aes(xintercept = end, linetype = "desert boundary"), color = "red") +

        geom_hline(yintercept = mean(.$ancient, na.rm = TRUE), color = "blue", linetype = "dashed", alpha = 0.5) +
        geom_hline(yintercept = mean(.$modern, na.rm = TRUE), color = "orange", linetype = "dashed", alpha = 0.5) +

        scale_color_manual(values = c("blue", "orange")) +
        guides(color = guide_legend("", override.aes = list(size = 5)),
               linetype = guide_legend("")) +
        labs(x = "genomic coordinate [bp]", y = "proportion of Neanderthal ancestry") +
        scale_x_continuous(labels = scales::comma) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 0.01)) +
        coord_cartesian(ylim = c(0, 0.5)) +
        scale_linetype_manual(values = "dashed") +
        theme_minimal() +
        theme(legend.position = "bottom", text = element_text(size = 13)) +
        ggtitle(paste("Archaic ancestry desert on chromosome", gsub("chr", "", .$chrom[1])))
    }
    if (!full) {
      p <- p +
        geom_point(data = filter(ancestry_df, ancient > 0), aes(midpoint, ancient, color = "ancient individuals"), size = 0.8) +
        geom_point(data = filter(ancestry_df, modern > 0), aes(midpoint, modern, color = "present-day individuals"), size = 0.8)
    }
    p
}

plot_desert_correlation <- function(ancestry_gr, chrom) {
  ancestry_df <- as_tibble(ancestry_gr) %>% filter(seqnames == chrom)

  rho <- ancestry_gr %>% filter(within_desert, seqnames == chrom) %>% { cor(.$modern, .$ancient) }

  ggplot() +
    geom_smooth(data = filter(ancestry_df, within_desert, ancient > 0, modern > 0), aes(modern, ancient, color = within_desert),
                formula = y ~ x, color = "red", fill = "black", method = "lm", linetype = "dashed", linewidth = 0.8, alpha = 0.35) +

    geom_point(data = filter(ancestry_df, !within_desert, ancient > 0, modern > 0), aes(modern, ancient, color = desert, shape = "outside desert"),
               color = "lightgray", alpha = 0.5) +
    geom_point(data = filter(ancestry_df, within_desert, ancient > 0, modern > 0), aes(modern, ancient, color = within_desert, shape = "within desert"), color = "black") +

    geom_abline(slope = 1, linetype = "dashed") +

    geom_vline(aes(color = "modern", xintercept = mean(ancestry_df$modern, na.rm = TRUE)), linetype = "dashed", color = "blue") +
    geom_hline(aes(color = "ancient", yintercept = mean(ancestry_df$ancient, na.rm = TRUE)), linetype = "dashed", color = "orange") +

    scale_x_log10(breaks = c(0.0001, mean(ancestry_df$modern, na.rm = TRUE), 1), labels = scales::percent_format(accuracy = 0.01), limits = c(0.00001, 1)) +
    scale_y_log10(breaks = c(0.0001, mean(ancestry_df$ancient, na.rm = TRUE), 1), labels = scales::percent_format(accuracy = 0.01), limits = c(0.00001, 1)) +
    labs(x = "Neanderthal ancestry proportion\nin present-day Eurasians [log scale]",
         y = "Neanderthal ancestry proportion\nin ancient Eurasians [log scale]") +
    coord_fixed() +
    theme_minimal() +
    theme(legend.position = "bottom", text = element_text(size = 13)) +
    guides(shape = guide_legend("window", override.aes = list(alpha = 1, size = 3)),
           linetype = "none") +
    scale_shape_manual(values = c(4, 20)) +
    ggtitle("", subtitle = paste0("Pearson correlation within desert = ",
                                  formatC(rho, format = "f", digits = 3)))
}

plot_desert_ancestry2 <- function(ancestry_gr, deserts_gr, chrom) {
  desert_df <- deserts_gr %>% filter(seqnames == chrom) %>% as_tibble

  ancestry_df <- as_tibble(ancestry_gr) %>%
    filter(seqnames == chrom) %>%
    select(chrom = seqnames, start, end, modern, chen, gap, midpoint)

  ancestry_df %>%
    filter(start >= (desert_df$start * 0.9) & end <= (desert_df$end * 1.1)) %>%
    {
      ggplot(data = .) +
        geom_rect(data = desert_df, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = "red", alpha = 0.1) +

        geom_line(aes(midpoint, chen), color = "orange") +
        geom_point(data = filter(., chen > 0), aes(midpoint, chen, color = "Chen et al."), size = 0.8) +

        geom_line(aes(midpoint, modern), color = "blue") +
        geom_point(data = filter(., modern > 0), aes(midpoint, modern, color = "Alba et al."), size = 0.8) +

        geom_vline(data = desert_df, aes(xintercept = start, linetype = "desert boundary"), color = "red") +
        geom_vline(data = desert_df, aes(xintercept = end, linetype = "desert boundary"), color = "red") +

        geom_hline(yintercept = mean(.$chen), color = "orange", linetype = "dashed", alpha = 0.5) +
        geom_hline(yintercept = mean(.$modern), color = "blue", linetype = "dashed", alpha = 0.5) +

        scale_color_manual(values = c("Chen et al." = "orange", "Alba et al." = "blue")) +
        guides(color = guide_legend("", override.aes = list(size = 5)), linetype = guide_legend("")) +

        labs(x = "genomic coordinate [bp]", y = "proportion of Neanderthal ancestry") +
        scale_x_continuous(labels = scales::comma) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 0.01)) +
        # coord_cartesian(ylim = c(0, 0.1)) +
        scale_linetype_manual(values = "dashed") +
        theme_minimal() +
        theme(legend.position = "bottom") +
        ggtitle(paste("Archaic ancestry desert on chromosome", gsub("chr", "", .$chrom[1])))
    }

}

# LD-based admixture dating -------------------------------------------------------------------

# Define positions of archaic ancestry informative sites
generate_info_sites <- function(tracts_gr, interval) {
  sites_grl <- lapply(seqlevels(tracts_gr), function(chrom) {
    positions <- seq(from = 1, to = seqlengths(tracts_gr)[chrom], by = interval)

    gr <- GRanges(seqnames = chrom, ranges = IRanges(start = positions, end = positions))
    mcols(gr)$index <- seq_len(length(gr))

    gr
  }) %>% GRangesList()
  seqlevels(sites_grl) <- seqlevels(tracts_gr)
  seqlengths(sites_grl) <- seqlengths(tracts_gr)

  sites_grl
}

# Define list of pairs of sites at given distances
# (one element of the list for each distance bin)
collect_pairs <- function(sites_grl, distances, recmap = NULL, ncores = parallel::detectCores()/2) {

  if (is.null(recmap)) {
    chrom_lengths <- seqlengths(sites_grl)
  } else {
    chrom_lengths <- split(recmap, recmap$chrom) %>% sapply(function(recmap_chr) max(recmap_chr$posg))
  }

  chroms <- sapply(sites_grl, function(x) as.character(unique(seqnames(x))))

  chr_pairs <- lapply(chroms, function(chrom) {

    sites_gr <- sites_grl[chroms == chrom, ] %>% unlist
    if (is.null(recmap)) {
      site_pos <- start(sites_gr)
    } else {
      site_pos <- convert_genetic(recmap, as.data.table(sites_gr), "start", chrom = "seqnames")$start_gen
    }

    pairs <- parallel::mclapply(distances, function(distance) {

      pair1 <- c()
      pair2 <- c()

      # iterate through each site one by one...
      for (i in sites_gr$index) {
        index1 <- i
        # ... and find the index of the first site that is at a given distance
        index2 <- sites_gr[site_pos >= site_pos[i] + distance]$index[1]

        if (is.na(index2)) {
          if (chrom_lengths[chrom] < site_pos[i] + distance  + distance / 10)
            break
          else
            next
        }

        # otherwise record the indices of the pair of sites and proceed with searching
        # for the next pair
        pair1 <- c(pair1, index1)
        pair2 <- c(pair2, index2)
      }

      list(pair1 = pair1, pair2 = pair2)

    }, mc.cores = ncores)

    pairs

  })

  names(chr_pairs) <- chroms
  chr_pairs
}

# Compute covariances of allele states at pairs of sites
compute_tract_covariances <- function(tracts_gr, sites_grl, pairs) {
  lapply(seqlevels(sites_grl), function(chrom) {

    chrom_sites_gr <- sites_grl[seqlevels(sites_grl) == chrom, ] %>% unlist

    parallel::mclapply(unique(tracts_gr$name), function(name) {

      ind_tracts_gr <- tracts_gr %>% filter(name == !!name, seqnames == chrom)
      ind_sites_gr <- chrom_sites_gr

      # mark sites falling within an introgressed tract
      tract_overlaps <- queryHits(findOverlaps(ind_sites_gr, ind_tracts_gr))
      mcols(ind_sites_gr)$neand <- FALSE
      if (length(tract_overlaps) > 0)
        mcols(ind_sites_gr[tract_overlaps])$neand <- TRUE
      mcols(ind_sites_gr)$neand <- as.integer(mcols(ind_sites_gr)$neand)

      covariances <- sapply(seq_along(distances), function(i) {
        sites1 <- ind_sites_gr[pairs[[chrom]][[i]]$pair1]$neand
        sites2 <- ind_sites_gr[pairs[[chrom]][[i]]$pair2]$neand
        cov(sites1, sites2)
      })

      tibble(
        chrom = chrom,
        name = name,
        sample_age = unique(ind_tracts_gr$sample_age),
        distance = distances,
        covariance = covariances
      )
    }, mc.cores = detectCores()) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .)
}


# Compute covariances of allele states at pairs of sites
compute_match_covariances <- function(info_gt, pairs, metadata) {
  archaic_name <- "NEA_1"
  samples <- setdiff(colnames(info_gt), c("chrom", "pos", archaic_name))

  lapply(unique(info_gt$chrom), function(chrom) {

    chrom_info_gt <- info_gt[, .SD[chrom %in% ..chrom, ], .SDcols = !c("chrom", "pos")]

    parallel::mclapply(samples, function(name) {

      ind_matches <- chrom_info_gt[, .(match = .SD[, get(name) == .SD[, get(archaic_name)]])]

      covariances <- sapply(seq_along(distances), function(i) {
        sites1 <- ind_matches[pairs[[chrom]][[i]]$pair1]$match
        sites2 <- ind_matches[pairs[[chrom]][[i]]$pair2]$match
        cov(sites1, sites2)
      })

      tibble(
        chrom = chrom,
        name = name,
        sample_age = filter(metadata, name == gsub("_hap\\d", "", !!name))$sample_age,
        distance = distances,
        covariance = covariances
      )
    }, mc.cores = detectCores()) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .)
}

fit_exponential <- function(cov_df, distance) {
  cov_df <- cov_df %>% mutate(name = gsub("_hap\\d", "", name))

  distance <- match.arg(distance, choices = c("physical", "genetic"))

  grid_df <- expand_grid(name = unique(cov_df$name))

  fit_df <- lapply(1:nrow(grid_df), function(i) {
    name <- grid_df[i, ]$name

    data_df <- filter(cov_df, name == !!name) %>%
      group_by(name, sample_age, distance) %>%
      summarise(covariance = mean(covariance), .groups = "keep")

    lambda <- tryCatch({
      nls_res <- nls(covariance ~ SSasymp(distance, Asym, R0, lrc), data = data_df)
      exp(coef(nls_res)["lrc"])
    }, error = function(e) NA
    )

    if (!is.na(lambda)) {
      df <- tibble(
        distance = data_df$distance,
        covariance = predict(nls_res, newdata = data_df[, "distance"])
      )
    } else
      df <- NULL

    tibble(
      name = name,
      sample_age = data_df$sample_age[1],
      lambda = lambda,
      t_gens_before = ifelse(distance == "physical", lambda / 1e-8, lambda * 100),
      t_admix = t_gens_before * gen_time + sample_age,
      fit = list(df)
    )
  }) %>%
    do.call(rbind, .) %>%
    unnest(fit)
}


# conversion of physical distances to genetic distances ---------------------------------------

read_recmap <- function(path) {
  if (!dir.exists(path))
    stop("Path '", path, "' does not exist", call. = FALSE)

  files <- list.files(path, "plink.chr\\d+.*.map", full.names = TRUE)
  lapply(files, function(f) readr::read_tsv(f, col_names = c("chrom", "_", "posg", "pos"),
                                            show_col_types = FALSE)) %>%
    do.call(rbind, .) %>%
    mutate(chrom = paste0("chr", chrom))
}

read_recmap2 <- function(path, cumulative = TRUE) {
  if (!dir.exists(path))
    stop("Path '", path, "' does not exist", call. = FALSE)

  files <- list.files(path, ".*_recombination_map_hg19_chr_\\d+.bed", full.names = TRUE)
  lapply(files, function(f) {
    df <- read_tsv(f, show_col_types = FALSE) %>% setNames(c("chrom", "start", "end", "rate"))
    if (cumulative) {
      df <- rbind(tibble(chrom = df[1, ]$chrom, start = 0, end = df[1, ]$start, rate = 0), df) %>%
        mutate(posg = cumsum((end - start) * rate * 100)) %>%
        select(chrom, posg, pos = end)
    } else
      df <- mutate(df, length = end - start) %>% select(chrom, start, end, length, rate)
    df
  }) %>% do.call(rbind, .)
}

convert_genetic <- function(recmap, df, cols, chrom = "chrom") {
  interpolators <- recmap %>%
    split(.$chrom) %>%
    lapply(function(chrom_map) approxfun(chrom_map$pos, chrom_map$posg, rule = 2))

  df %>%
    { split(., .[[chrom]]) } %>%
    lapply(function(chrom_df) {
      if (!nrow(chrom_df)) return(NULL)
      chrom <- as.integer(gsub("chr", "", chrom_df[[chrom]][1]))
      for (c in cols) {
        chrom_df[[paste0(c, "_gen")]] <- interpolators[[!!chrom]](chrom_df[[c]])
      }
      chrom_df
    }) %>%
    do.call(rbind, .)
}

# Functions to compute tract frequency ------------------------------------

# This function computes the frequency of Neanderthal tract across windows
# NOTE: it returns them as a vector, ideally to be added at the windows_gr object as an metadata column

compute_tract_freq <- function(tracts_gr, windows_gr, age_group = NULL) {
  if(!is.null(age_group)){
    tracts_gr <- tracts_gr[tracts_gr$age_group==age_group]
  }
  autosomes <- getChromInfoFromUCSC("hg19") %>%
    filter(grepl("chr\\d+$", chrom)) %>%
    mutate(start = 1) %>%
    dplyr::rename(end = size) %>%
    makeGRangesFromDataFrame()
  
  # first compute coverage...
  cov <- coverage(tracts_gr)
  # ... then convert that to proportions
  cov <- lapply(cov, function(x) x / length(unique(tracts_gr$sampleId)))
  
  ancestry_list <- mclapply(as.character(unique(seqnames(windows_gr))),
                            function(chrom) {
                              chrom_coverage <- cov[[chrom]]
                              chrom_gr <- windows_gr[seqnames(windows_gr) == chrom]
                              
                              # count overlaps between windows and tracts
                              # check Views object
                              cov_view <- Views(chrom_coverage, ranges(chrom_gr))
                              average_coverage_per_window <- viewMeans(cov_view)
                              average_coverage_per_window
                            }, mc.cores = 1)
  
  # ancestry_grl <- GRangesList(ancestry_list)
  unlist(ancestry_list)
}

# Taken from Signac
StringToGRanges <- function(regions, sep = c("-", "-"), ...) {
  ranges.df <- data.frame(ranges = regions)
  ranges.df <- separate(
    data = ranges.df,
    col = "ranges",
    sep = paste0(sep[[1]], "|", sep[[2]]),
    into = c("chr", "start", "end")
  )
  granges <- makeGRangesFromDataFrame(df = ranges.df, ...)
  return(granges)
}


# Plot frequencies over time ----------------------------------------------

# Take in inout a GR object

plot_frequencies <- function(gr_obj){
  # Get frequencies as df
  meta_df <- as.matrix(mcols(gr_obj))
  rownames(meta_df) <- paste0(seqnames(gr_obj), ":",ranges(gr_obj))
  meta_df <- meta_df %>% as_tibble(rownames="bin")
  # Long format for plotting
  meta_long <- meta_df %>% pivot_longer(cols = colnames(meta_df)[-c(1,ncol(meta_df))]) %>% transform(name=factor(name, levels=rev(c("present-day", "(0,2e+03]", "(2e+03,5e+03]", "(5e+03,1e+04]", "(1e+04,1.2e+04]" ,"(1.2e+04,3e+04]", "(3e+04,Inf]"))))
  p <- ggplot()+geom_line(data = meta_long, aes(x = name, y = as.numeric(value), group=bin, color=gene))+theme_minimal()
  print(p)
}


# Topic modeling ----------------------------------------------------------

compute_bin_vec <- function(tracts_gr, windows_gr, sample_id = NULL) {
  if(is.null(sample_id)){
    stop("Forgot to specify sampleID.")
  }
  
  # Get sample
  sample_gr <- tracts_gr %>% filter(name==sample_id) 
  seqinfo(sample_gr) <- seqinfo(windows_gr)
  # first compute coverage...
  cov <- coverage(sample_gr)
  
  ancestry_list <- mclapply(as.character(unique(seqnames(windows_gr))),
                            function(chrom) {
                              chrom_coverage <- cov[[chrom]]
                              chrom_gr <- windows_gr[seqnames(windows_gr) == chrom]
                              
                              # count overlaps between windows and tracts
                              # check Views object
                              cov_view <- Views(chrom_coverage, ranges(chrom_gr))
                              average_coverage_per_window <- viewMeans(cov_view)
                              average_coverage_per_window
                            }, mc.cores = 1)
  
  # ancestry_grl <- GRangesList(ancestry_list)
  unlist(ancestry_list)
}

# for simulation
# only changes the name
compute_bin_vec_sim <- function(tracts_gr, windows_gr, sample_id = NULL) {
  if(is.null(sample_id)){
    stop("Forgot to specify sampleID.")
  }
  
  # Get sample
  sample_gr <- tracts_gr %>% filter(name==sample_id) 
  seqinfo(sample_gr) <- seqinfo(windows_gr)
  # first compute coverage...
  cov <- coverage(sample_gr)
  
  ancestry_list <- mclapply(as.character(unique(seqnames(windows_gr))),
                            function(chrom) {
                              chrom_coverage <- cov[[chrom]]
                              chrom_gr <- windows_gr[seqnames(windows_gr) == chrom]
                              
                              # count overlaps between windows and tracts
                              # check Views object
                              cov_view <- Views(chrom_coverage, ranges(chrom_gr))
                              average_coverage_per_window <- viewMeans(cov_view)
                              average_coverage_per_window
                            }, mc.cores = 1)
  
  # ancestry_grl <- GRangesList(ancestry_list)
  unlist(ancestry_list)
}

# function that returns a vector of which exact tracts are in that sample
get_exact_match <- function(target, uni_tracts_gr){
  # binary vector of results
  vector <- rep(0, length(uni_tracts_gr))
  names(vector) <- as.character(uni_tracts_gr)
  # pairs <- findOverlapPairs(target, uni_tracts_gr, type="equal")
  res <- findOverlaps(target, uni_tracts_gr, type = "equal")
  vector[subjectHits(res)] <- 1
  vector
}

plot_on_map <- function(metadata){
  world <- ne_countries(scale = "medium", returnclass = "sf")
  sf::st_agr(world) <- "constant"
  bbox <- st_as_sfc(st_bbox(c(xmin = -25, xmax = 65, ymin = 25, ymax = 70), crs = st_crs(world)))
  western_eurasia <- st_crop(st_make_valid(world), bbox)
  
  metadata %>%
    filter(!is.na(latitude) & !is.na(longitude)) %>%
    st_as_sf(coords = c("longitude", "latitude")) %>%
    st_set_crs(4326) %>%
    ggplot() +
    geom_sf(data = western_eurasia) +
    geom_sf(aes(color = clusterAlias, alpha=ageAverage)) +
    coord_sf(crs = 3035)
}


# tfs  --------------------------------------------------------------------

get_coords <- function(ts) {
  model <- attr(ts, "model")
  lapply(model$populations, function(pop) {
    if (slendr:::has_map(pop))
      sf::st_centroid(pop[, c("pop", "geometry")])
    else
      NULL
  }) %>%
    do.call(rbind, .) %>%
    as_tibble %>%
    sf::st_as_sf()
}

# function to plot it 
plot_tfs <- function(tfs){
  ggplot(data.frame(tfs))+geom_bar(stat="identity", aes(seq_along(tfs),tfs))+theme_minimal()+xlab("Frequency class")+ylab("Proportion of tracts")
}

# function to compute the tract frequency spectrum
# note that I will consider the samples rather than the chromosome,
# I will have to talk with Martin for this
# I find really weird that IBDmix returns only the tracts, but not the chromosome
# I will not use the hist function
compute_tfs <- function(bin_mat){
  n_sample <- ncol(bin_mat)
  tfs <- rep(0, n_sample)
  # find how many individual share a tract
  counts <- table(rowSums(bin_mat))
  # remove the 0 counts if there are
  if(!is.na(counts["0"])) counts <- counts[-1]
  tfs[as.numeric(names(counts))] <- counts
  # normalize to get proportion of tracts
  tfs <- tfs/sum(tfs)
  tfs
}

# tract matrices approaches -----------------------------------------------

# uniqueness approach

get_bin_mat_unique <- function(tracts_gr){
  tracts_unique_gr <- unique(tracts_gr)
  if(length(tracts_unique_gr)>1){
    bin_mat <- sapply(unique(tracts_gr$name), function(x) get_exact_match(filter(tracts_gr, name==x), tracts_unique_gr))
  }else{
    bin_mat <- matrix(1, nrow = 1, ncol = length(unique(tracts_gr$name)))
    rownames(bin_mat) <- as.character(tracts_unique_gr)
    colnames(bin_mat) <- unique(tracts_gr$name)
  }
  return(bin_mat)
}

# windows approach

get_bin_mat_windows <- function(tracts_gr, window_size = 50e3, step_size = 50e3, len_chr = NULL){
  if (is.null(len_chr)) {
    stop("Specify length of chromosome.")
  }
  windows_gr <- generate_windows_sim(window_size = window_size, step_size = step_size, len_chr = len_chr)  
  bin_mat <- sapply(unique(tracts_gr$name), function(sample_id){
    compute_bin_vec_sim(tracts_gr, windows_gr, sample_id)
  })
  rownames(bin_mat) <- paste0(seqnames(windows_gr), ":",ranges(windows_gr))
  bin_mat <- ifelse(bin_mat>0, 1, 0)
  bin_mat <- bin_mat[rowSums(bin_mat)!=0,]
  return(bin_mat)
}

get_bin_mat_windows_emp <- function(tracts_gr, window_size = 50e3, step_size = 50e3){
  windows_gr <- generate_windows_no_gap(window_size = window_size, step_size = step_size)  
  bin_mat <- sapply(unique(tracts_gr$name), function(sample_id){
    compute_bin_vec(tracts_gr, windows_gr, sample_id)
  })
  rownames(bin_mat) <- paste0(seqnames(windows_gr), ":",ranges(windows_gr))
  bin_mat <- ifelse(bin_mat>0, 1, 0)
  bin_mat <- bin_mat[rowSums(bin_mat)!=0,]
  return(bin_mat)
}

# function to get which subtracts are contained in an inidvidual
# NOTE: the function is meant to work with one sample at a time
get_tract_on <- function(tracts_gr, sample_id){
  tract_df <- tracts_gr %>% as.data.frame() %>% select(start, end)
  # convert to matrix
  tract_mat <- tract_df %>% as.matrix()
  # get vector of starts and ends
  starts_ends <- tract_mat %>% c() %>% unique() %>% sort()
  # it is meant to be for only one chrom, so it the gr has to be only one chr
  seqname <- seqnames(tracts_gr) %>% as.character() %>% unique()
  # The start end vector is sorted, so for the start I will exclude the last element, and fot he ends of the ranges the first one
  mask_gr <- GRanges(seqnames = seqname, ranges = IRanges(start = starts_ends[-length(starts_ends)], end = starts_ends[-1]))
  names <- as.character(mask_gr)
  vec <- rep(0, length(names))
  # get which ranges are in that sample
  samp_tracts_gr <- tracts_gr %>% filter(name==sample_id)
  hit <- subjectHits(findOverlaps(samp_tracts_gr, mask_gr, minoverlap = 2L))
  vec[hit] <- 1
  names(vec) <- names
  return(vec)
}

# function to get the subset of the tracts
# 3rd approach
get_bin_mat_subset <- function(tracts_gr){
  bin_mat <- sapply(unique(tracts_gr$name), function(sample_id){
    get_tract_on(tracts_gr, sample_id)
  })
  # return NULL if the bin mat has nothing inisde
  # this to avoid the problem of rowSums over a NULL elements
  if(is.null(dim(bin_mat))) return(NULL)
  bin_mat <- bin_mat[rowSums(bin_mat)!=0,]
  return(bin_mat)
}

get_bin_mat_subset_emp <- function(tracts_gr){
  # iterate over each chromosome
  # to get one matrix for each chromosome
  col_names <- unique(tracts_gr$name) 
  bin_mat_l <- lapply(as.character(unique(seqnames(tracts_gr))),
                        function(chrom) {
                          chrom_gr <- filter(tracts_gr, seqnames==chrom)
                          mat <- get_bin_mat_subset(chrom_gr)
                          # This is because the some inds might not have tracts 
                          # for some chromosomes and in order to create one big matrix it is important
                          # to have the dimensions consistent
                          full_matrix <- matrix(0, nrow = nrow(mat), ncol = length(col_names), dimnames = list(rownames(mat), col_names))
                          full_matrix[rownames(mat), colnames(mat)] <- mat
                          full_matrix
                        })
  bin_mat <- do.call(rbind, bin_mat_l)
  return(bin_mat)
}

# 4th approach
get_bin_mat_sites <- function(tracts_gr){
  bin_mat <- sapply(unique(tracts_gr$name), function(sample_id){
    get_rec_sites(tracts_gr, sample_id)
  })
  bin_mat <- bin_mat[rowSums(bin_mat)!=0,]
  return(bin_mat)
}

get_rec_sites <- function(tracts_gr, sample_id){
  tract_df <- tracts_gr %>% as.data.frame() %>% select(start, end)
  # convert to matrix
  tract_mat <- tract_df %>% as.matrix()
  # get vector of starts and ends
  starts_ends <- tract_mat %>% c() %>% unique() %>% sort()
  # it is meant to be for only one chrom, so it the gr has to be only one chr
  seqname <- seqnames(tracts_gr) %>% as.character() %>% unique()
  # The start end vector is sorted, so for the start I will exclude the last element, and fot he ends of the ranges the first one
  names <- paste0(seqname, ":", as.character(starts_ends))
  vec <- rep(0, length(names))
  # get which ranges are in that sample
  samp_tracts_gr <- tracts_gr %>% filter(name==sample_id)
  samp_start_ends <- samp_tracts_gr %>% as.data.frame() %>% select(start, end) %>% as.matrix() %>% c() %>% unique() %>% sort()
  names(vec) <- names
  vec[paste0(seqname, ":", as.character(samp_start_ends))] <- 1
  return(vec)
}
# plot functions ----------------------------------------------------------

##########################
# allele frequency spectrum

plot_afs <- function(afs, title){
  # title <- unlist(strsplit(title, "_"))
  p <- ggplot()+geom_bar(stat="identity", aes(seq_along(afs),afs/sum(afs)))+theme_minimal()+ylab("Proportion of SNPs")+xlab("Derived allele counts")+ggtitle(title)
  p  
}

# afs_pop is a string with a pattern like this
# afs_eur
get_afs_df <- function(model_list, afs_pop){
  afs_df <- sapply(names(model_list), function(x) {
    afs <- model_list[[x]][[afs_pop]]
    # https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.allele_frequency_spectrum
    # check this to understand why I am removing the first and the last
    afs <- afs[-1]
    afs <- afs[-length(afs)]
    afs
  })
  # names(afs_list) <- names(model_list)  
  sweep(afs_df, MARGIN = 2, FUN = "/", colSums(afs_df))
}

###############################

## heterozygosity

get_het_df <- function(model_list){
  pops <- sort(model_list$Ne_OOA_100$het$set)
  het_df <- do.call(cbind, lapply(model_list, function(x){
    het <- x[["het"]] %>% arrange(set) %>%  .[["diversity"]]
    het}))
  rownames(het_df) <- pops
  het_df
}


###########################

# fst

get_fst_df <- function(model_list){
  pops_df <- model_list$Ne_OOA_100$fst %>% arrange(x, y) %>% select(x, y)
  fst_df <- do.call(cbind, lapply(model_list, function(x){ 
    fst <- x[["fst"]] %>% arrange(x, y) %>% .[["Fst"]] 
    fst}))
  cbind(pops_df, fst_df)
}

##########################

# tract based statistics


# uniqueness approach
# return a list od matrices where each column in the matrices
# correspond to a Ne of the OOA and each eleemnt of the list
# to a run
get_unique_tfs <- function(model_runs, n_samples = 100) {
  unique_tfs_df <- lapply(model_runs, function(x) {
    sapply(x, function(y) {
      tracts <- y$tracts
      if (nrow(tracts) == 0 |
          (summarise(tracts, n_ind = n_distinct(name)) %>% .[["n_ind"]] < 2)) {
        return(rep(0, n_samples))
      }
      # add chrom for simplicity
      tracts$chrom <- "chr1"
      tracts_gr <- makeGRangesFromDataFrame(
        tracts,
        keep.extra.columns = TRUE,
        seqnames.field = "chrom",
        start.field = "left",
        end.field = "right"
      )
      bin_mat <- get_bin_mat_unique(tracts_gr)
      if(is.null(dim(bin_mat))) return(rep(0, n_samples))
      tfs_unique <- compute_tfs(bin_mat)
      c(tfs_unique, rep(0, n_samples-length(tfs_unique)))
    })
  })
  unique_tfs_df
}

# windows approach
get_windows_tfs <- function(model_runs, n_samples = 100, len_chr = 1e6) {
  windows_tfs_df <- lapply(model_runs, function(x) {
    sapply(x, function(y) {
      tracts <- y$tracts
      if (nrow(tracts) == 0 | (summarise(tracts, n_ind = n_distinct(name)) %>% .[["n_ind"]] < 2)) {
        return(rep(0, n_samples))
      }
      # add chrom for simplicity
      tracts$chrom <- "chr1"
      tracts_gr <- makeGRangesFromDataFrame(
        tracts,
        keep.extra.columns = TRUE,
        seqnames.field = "chrom",
        start.field = "left",
        end.field = "right"
      )
      bin_mat <- get_bin_mat_windows(tracts_gr, len_chr = len_chr)
      if(is.null(dim(bin_mat))) return(rep(0, n_samples))
      tfs_windows <- compute_tfs(bin_mat)
      c(tfs_windows, rep(0, n_samples - length(tfs_windows)))
    })
  })
  windows_tfs_df
}

# subset approach
get_subset_tfs <- function(model_runs, n_samples = 100){
  subset_tfs_df <- lapply(model_runs, function(x) {
    sapply(x, function(y){
      tracts <- y$tracts
      if(nrow(tracts)==0 | (summarise(tracts, n_ind=n_distinct(name)) %>% .[["n_ind"]]<2)){
        return(rep(0, n_samples))
      }
      # add chrom for simplicity
      tracts$chrom <- "chr1"
      tracts_gr <- makeGRangesFromDataFrame(
        tracts,
        keep.extra.columns = TRUE,
        seqnames.field = "chrom",
        start.field = "left",
        end.field = "right"
      )
      bin_mat <- get_bin_mat_subset(tracts_gr)
      if(is.null(dim(bin_mat))) return(rep(0, n_samples))
      # print(bin_mat)
      tfs_subset <- compute_tfs(bin_mat)
      c(tfs_subset, rep(0, n_samples-length(tfs_subset)))
    })
  })
  subset_tfs_df
}

plot_tfs <- function(tfs, title){
  # title <- unlist(strsplit(title, "_"))
  p <- ggplot()+geom_bar(stat="identity", aes(seq_along(tfs),tfs/sum(tfs)))+theme_minimal()+ylab("Proportion of tracts")+xlab("Frequency class")+ggtitle(title)
  p  
}

# Function to average over the runs
get_avg_tfs <- function(tfs_l){
  Reduce("+", tfs_l)/length(tfs_l)
}

get_avg_afs <- function(afs_l) {
  afs_avg_df <- Reduce("+", afs_l) / length(afs_l)
  afs_avg_df %>% as_tibble() %>%
    mutate(bin = 1:nrow(.)) %>% pivot_longer(cols = starts_with("Ne"),
                                             names_to = "Ne",
                                             values_to = "Freq") %>%  mutate(across(Ne, ~ factor(
                                               .,
                                               levels = c(
                                                 "Ne_OOA_100",
                                                 "Ne_OOA_500",
                                                 "Ne_OOA_1000",
                                                 "Ne_OOA_2000",
                                                 "Ne_OOA_5000",
                                                 "Ne_OOA_10000"
                                               )
                                             )))
  
}

get_fst_df_long <- function(model_runs) {
  fst_df_list <- lapply(model_runs, function(x) {
    get_fst_df(x)
  })
  tmp_list <- lapply(seq_along(fst_df_list), function(i)
    fst_df_list[[i]] %>% as.data.frame() %>% pivot_longer(
      cols = starts_with("Ne"),
      names_to = "Ne",
      values_to = "fst"
    )%>% mutate(run = i)) 
  
  do.call(rbind, tmp_list)
}

# function to convert a list into a long df
convert_list_to_long_df <- function(df_list){
  tmp_list <- lapply(seq_along(df_list), function(i) df_list[[i]] %>% as_tibble() %>% mutate(bin = 1:nrow(.), run = i))
  do.call(rbind, tmp_list) %>% pivot_longer(cols = starts_with("Ne"), names_to = "Ne", values_to = "Freq")
}

get_fst_avg_df <- function(model_runs) {
  fst_df_list <- lapply(model_runs, function(x) {
    get_fst_df(x)
  })
  fst_avg <- cbind(select(fst_df_list[[1]], x, y), Reduce("+", lapply(fst_df_list, function(df) select(df, !c("x", "y"))))/length(fst_df_list))
  fst_avg
}

get_het_df_long <- function(model_runs) {
  het_df_list <- lapply(model_runs, function(x) {
    get_het_df(x)
  })
  tmp_list <- lapply(seq_along(het_df_list), function(i)
    het_df_list[[i]] %>% as.data.frame() %>% rownames_to_column("pop") %>% pivot_longer(
      cols = starts_with("Ne"),
      names_to = "Ne",
      values_to = "het"
    )%>% mutate(run = i)) 
  
  do.call(rbind, tmp_list)
}

get_het_avg_df <- function(model_runs) {
  het_df_list <- lapply(model_runs, function(x) {
    get_het_df(x)
  })
  Reduce("+", het_df_list)/length(het_df_list)
}

# function to simplify the tracts from simulations
# to mimick the ones of IBDmix
simplify_tracts <- function(tracts_df, size_cutoff = 50e3){
  tracts$chrom <- "chr1"
  tracts_gr <- makeGRangesFromDataFrame(
    tracts,
    keep.extra.columns = TRUE,
    seqnames.field = "chrom",
    start.field = "left",
    end.field = "right"
  )
  tracts_grl <- lapply(unique(tracts_gr$name), function(id){
    # extract tract for individual
    ind_gr <- tracts_gr %>% filter(name == id)
    # create a df for the metacolumns, they will be lost otherwise
    # and remove some columns because not needed or wrong
    df <- mcols(ind_gr) %>% as.data.frame() %>% select(-one_of("node_id", "length", "haplotype"))
    # reduce tracts to simuate what IBDmix would do
    # since we do not have the haplotype in output
    res_gr <- reduce(ind_gr)
    # add metadata colummns
    mcols(res_gr) <- df[length(res_gr),]
    res_gr$length <- width(res_gr)
    # filter out tracts less than cutoff, gain to mimick IBDmix
    filter(res_gr, length>=size_cutoff)
  })
  do.call(c, tracts_grl)
}


# function to get a 2d tfs
get_2d_tfs <- function(bin_mat, pop1_ids, pop2_ids){
  bin_mat_pop1 <- bin_mat[, pop1_ids]
  bin_mat_pop2 <- bin_mat[, pop2_ids]
  # get matrix with frequencies in each pop per bin
  freq_pop1 <- rowSums(bin_mat_pop1)
  freq_pop2 <- rowSums(bin_mat_pop2)
  mat <- matrix(c(freq_pop1, freq_pop2), ncol = 2)
  res <- matrix(0, nrow = length(pop1_ids)+1, ncol = length(pop2_ids)+1)
  for (row in 1:nrow(mat)) {
    ind <- mat[row,]
    res[ind[1]+1,ind[2]+1] <- res[ind[1]+1,ind[2]+1]+1
  }
  # mat
  colnames(res) <- 0:(ncol(res)-1)
  rownames(res) <- 0:(nrow(res)-1)
  res[res==0] <- NA
  res
}