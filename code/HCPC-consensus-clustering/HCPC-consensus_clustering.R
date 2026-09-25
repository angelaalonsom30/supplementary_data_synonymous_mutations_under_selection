
# Consensus clustering + PCA factor map + pheatmap heatmap (K=2-7)

library(here)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(FactoMineR)
library(pheatmap)
library(ggrepel)
library(readxl)


df <- read_excel("4_fold_sites_in_gtsB.xlsx") ##here the file available as supplementary Table S3 (for only four-fold sites
                                            # or Table S8 (including four-fold sites and neighboring nonsynonymous sites)

gene_col      <- "gene"        
site_id_col   <- "nucleotide_position"        
codon_pos_col <- "order_codon"   

gene_name <- "gtsB"  #the gene to evaluate

active_cols <- c(
  "substitution rate" ,"transitions_avg" ,"transversions_avg","T-ended synonymous codon frequency" ,"A-ended synonymous codon frequency", "G-ended synonymous codon frequency" ,"C-ended synonymous codon frequency"
)

##Parameters
ks <- 2:7   #the number of Clusters to test
B <- 5000   #iterations   
frac <- 0.8  #here define the resampling fraction: 0.8, 0.7 or 0.9
seed <- 1
scale_unit <- TRUE

out_dir <- "out_gtsB_consensus_final_final_4fold"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

#Filtering codon positions 2 & 3
df2 <- df %>% filter(.data[[codon_pos_col]] %in% c(3))  ## c(3) only when analyzing four-fold synonymous sites, for neighboring and 4-fols sites
                                                          ##replace for c(2,3)
df_g <- df2 %>% filter(.data[[gene_col]] == gene_name)
stopifnot(nrow(df_g) >= 10)

#storing outputs
site_out <- df_g %>% select(all_of(c(gene_col, site_id_col, codon_pos_col)))
metrics_rows <- list()
cdf_rows <- list()
A_by_k <- numeric(length(ks)); names(A_by_k) <- ks


###Functions

# PCA on X (sites x variables), keep PC1-PC2, Ward clustering
cluster_pc12_ward <- function(X, k, scale_unit = TRUE) {
  pca <- FactoMineR::PCA(X, ncp = 2, scale.unit = scale_unit, graph = FALSE)
  coords <- pca$ind$coord[, 1:2, drop = FALSE]
  hc <- hclust(dist(coords), method = "ward.D2")
  cutree(hc, k = k)
}

#Stratified sampling when 2nd and 3rd codon positions are compared
stratified_sample <- function(strata, frac) {
  strata <- as.factor(strata)
  idx <- integer(0)
  for (lv in levels(strata)) {
    ids <- which(strata == lv)
    m <- max(2, floor(frac * length(ids)))
    m <- min(m, length(ids))
    idx <- c(idx, sample(ids, m))
  }
  unique(idx)
}

#Build consensus matrix 
build_consensus <- function(df_g, active_cols, k, B, frac, seed, scale_unit, strata_col) {
  set.seed(seed)
  X <- as.matrix(df_g[, active_cols, drop = FALSE])
  n <- nrow(X)
  strata <- df_g[[strata_col]]
  
  num <- matrix(0, n, n)  # how often i,j co-clustered
  den <- matrix(0, n, n)  # how often i,j co-occurred in same resample
  
  for (b in seq_len(B)) {
    idx <- stratified_sample(strata, frac)
    cl  <- cluster_pc12_ward(X[idx, , drop = FALSE], k = k, scale_unit = scale_unit)
    
    den[idx, idx] <- den[idx, idx] + 1
    
    for (c in seq_len(k)) {
      mem <- idx[which(cl == c)]
      if (length(mem) < 2) next
      num[mem, mem] <- num[mem, mem] + 1
    }
  }
  
  C <- num / pmax(den, 1)
  diag(C) <- 1
  C
}

#Final clustering + metrics from consensus
analyze_consensus <- function(C, k) {
  hcC <- hclust(as.dist(1 - C), method = "ward.D2")
  cl_final <- cutree(hcC, k = k)
  
  within_means <- sapply(seq_len(k), function(c) {
    members <- which(cl_final == c)
    if (length(members) < 2) return(NA_real_)
    sub <- C[members, members, drop = FALSE]
    mean(sub[upper.tri(sub)])
  })
  overall_within <- mean(within_means, na.rm = TRUE)
  
  vals <- C[upper.tri(C, diag = FALSE)]
  pac <- mean(vals > 0.1 & vals < 0.9)
  
  stability <- sapply(seq_along(cl_final), function(i) {
    mates <- which(cl_final == cl_final[i])
    mates <- setdiff(mates, i)
    if (length(mates) == 0) return(NA_real_)
    mean(C[i, mates])
  })
  
  list(cluster = cl_final,
       within_means = within_means,
       overall_within = overall_within,
       pac = pac,
       stability = stability)
}

#Monti A(K): area under CDF of consensus values (off-diagonal)
cdf_area_monti <- function(C) {
  v <- C[upper.tri(C, diag = FALSE)]
  1 - mean(v)
}

cdf_curve <- function(C, grid = seq(0, 1, by = 0.01)) {
  v <- C[upper.tri(C, diag = FALSE)]
  data.frame(x = grid, F = ecdf(v)(grid))
}

#Plots

#pheatmap consensus matrix with dendrogram with cluster bars
plot_consensus_heatmap_pheatmap <- function(C, cluster, outfile, title,
                                            hc_method = "average",
                                            fontsize = 8,
                                            fontsize_title = 10,
                                            legend_breaks = c(0, 0.5, 1),
                                            width_in = 8, height_in = 6,
                                            res = 600) {
  
  C <- (C + t(C)) / 2
  diag(C) <- 1
  if (is.null(rownames(C))) rownames(C) <- seq_len(nrow(C))
  if (is.null(colnames(C))) colnames(C) <- seq_len(ncol(C))
  
  hc <- hclust(as.dist(1 - C), method = hc_method)
  
  cl <- factor(cluster)
  ann <- data.frame(Cluster = cl)
  rownames(ann) <- rownames(C)
  
  lev <- levels(cl)
  cl_cols <- setNames(grDevices::rainbow(length(lev)), lev)
  ann_cols <- list(Cluster = cl_cols)
  
  cols <- grDevices::colorRampPalette(c("white", "#08306B"))(100)
  
  grDevices::png(outfile, width = width_in, height = height_in, units = "in", res = res)
  pheatmap::pheatmap(
    C,
    color = cols,
    breaks = seq(0, 1, length.out = 101),
    legend_breaks = legend_breaks,
    legend_labels = legend_breaks,
    cluster_rows = hc,
    cluster_cols = hc,
    show_rownames = FALSE,
    show_colnames = FALSE,
    annotation_row = ann,
    annotation_col = ann,
    annotation_colors = ann_cols,
    border_color = NA,
    main = title,
    fontsize = fontsize,
    fontsize_main = fontsize_title
  )
  grDevices::dev.off()
}

#Factor map PC1-PC2 (PCA on original X), color by consensus cluster, ellipses + centroids
plot_factor_map_pc12 <- function(df_g, active_cols, cluster, outfile, title,
                                 scale_unit = TRUE,
                                 ellipse_level = 0.80,
                                 point_size = 2.0,
                                 centroid_size = 4.0,
                                 centroid_stroke = 1.4,
                                 width_in = 8, height_in = 5, dpi = 800) {
  
  pca <- FactoMineR::PCA(df_g[, active_cols, drop = FALSE],
                         ncp = 2, scale.unit = scale_unit, graph = FALSE)
  
  coords <- as.data.frame(pca$ind$coord[, 1:2, drop = FALSE])
  colnames(coords) <- c("PC1", "PC2")
  coords$cluster <- factor(cluster)
  
  var_exp <- pca$eig[1:2, 2]
  xlab <- sprintf("Dim 1 (%.2f%%)", var_exp[1])
  ylab <- sprintf("Dim 2 (%.2f%%)", var_exp[2])
  
  centroids <- aggregate(cbind(PC1, PC2) ~ cluster, coords, mean)
  
  p <- ggplot(coords, aes(PC1, PC2, color = cluster)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6, color = "grey70") +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.6, color = "grey70") +
    geom_point(size = point_size, alpha = 0.9) +
    stat_ellipse(aes(group = cluster),
                 type = "norm", level = ellipse_level,
                 linetype = "dashed", linewidth = 0.7,
                 color = "grey40",
                 show.legend = FALSE) +
    geom_point(data = centroids, aes(PC1, PC2),
               shape = 4, size = centroid_size, stroke = centroid_stroke,
               inherit.aes = FALSE, color = "black") +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(),
          legend.position = "right",
          legend.title = element_blank()) +
    labs(title = title, x = xlab, y = ylab)
  
  ggsave(outfile, p, width = width_in, height = height_in, dpi = dpi)
}

# Stability histogram
plot_stability_hist <- function(stability, outfile, title) {
  p <- ggplot(data.frame(stability = stability), aes(stability)) +
    geom_histogram(bins = 30) +
    theme_minimal(base_size = 12) +
    labs(title = title, x = "Mean consensus with cluster mates", y = "Count")
  ggsave(outfile, p, width = 6.5, height = 4.5, dpi = 800)
}


#Running the analysis
for (k in ks) {
  message("gtsB: k = ", k)
  
  #consensus matrix
  C <- build_consensus(df_g, active_cols, k, B, frac, seed, scale_unit, strata_col = codon_pos_col)
  
  #clustering + metrics
  res <- analyze_consensus(C, k)
  
  #Monti A(K) + CDF curve
  A_by_k[as.character(k)] <- cdf_area_monti(C)
  cdf_rows[[paste0("k", k)]] <- cdf_curve(C) %>% mutate(k = k)
  
  #site-level outputs
  site_out[[paste0("consensus_cluster_k", k)]] <- res$cluster
  site_out[[paste0("stability_k", k)]] <- res$stability
  
  # plots: pheatmap consensus + factor map + stability hist
  plot_consensus_heatmap_pheatmap(
    C, res$cluster,
    file.path(out_dir, paste0("consensus_pheatmap_k", k, ".png")),
    paste0("Consensus matrix (k=", k, ")")
  )
  
  plot_factor_map_pc12(
    df_g, active_cols, res$cluster,
    file.path(out_dir, paste0("factor_map_pc12_like_example_k", k, ".png")),
    paste0("(k=", k, ")"),
    scale_unit = scale_unit,
    ellipse_level = 0.80
  )
  
  plot_stability_hist(
    res$stability,
    file.path(out_dir, paste0("stability_hist_k", k, ".png")),
    paste0("Stability histogram (k=", k, ")")
  )
  
  #sites per cluster list
  cl_df <- data.frame(site_id = df_g[[site_id_col]], cluster = res$cluster)
  cl_lists <- cl_df %>%
    group_by(cluster) %>%
    summarise(
      gene_id = gene_name,
      k = k,
      n_sites = n(),
      sites = paste(site_id, collapse = ";"),
      .groups = "drop"
    ) %>%
    arrange(cluster)
  
  write_csv(cl_lists, file.path(out_dir, paste0("cluster_lists_k", k, ".csv")))
  
  #metrics row
  metrics_rows[[paste0("k", k)]] <- tibble(
    gene_id = gene_name,
    k = k,
    n_sites = nrow(df_g),
    B = B,
    frac = frac,
    overall_within_mean = res$overall_within,
    pac = res$pac,
    within_cluster_mean = paste(round(res$within_means, 3), collapse = ","),
    stability_mean = mean(res$stability, na.rm = TRUE),
    stability_sd = sd(res$stability, na.rm = TRUE),
    A_cdf = A_by_k[as.character(k)]
  )
}


#Exporting the outputs

metrics_df <- dplyr::bind_rows(metrics_rows) %>%
  dplyr::arrange(k)

metrics_df <- metrics_df %>%
  dplyr::arrange(k) %>%
  dplyr::mutate(
    deltaK = dplyr::case_when(
      k == 2 ~ A_cdf,
      TRUE   ~ (A_cdf - dplyr::lag(A_cdf)) / dplyr::lag(A_cdf)
    )
  )

write_csv(site_out, file.path(out_dir, "site_assignments_gtsB_4fold_k2_k7_fraction08.csv"))
write_csv(metrics_df, file.path(out_dir, "consensus_08_fraction_metrics_gtsB_4fold_k2_k7.csv"))


#Monti plots: A(K), ΔK, CDF curves

pA <- ggplot(metrics_df, aes(x = k, y = A_cdf)) +
  geom_line() + geom_point() +
  theme_minimal(base_size = 12) +
  labs(title = "gtsB: Area under consensus CDF A(K) (Monti et al. 2003)",
       x = "K", y = "A(K)")
ggsave(file.path(out_dir, "monti_A_of_K_gtsB_4fold_.png"), pA, width = 6.5, height = 4.5, dpi = 600)


pD <- ggplot(metrics_df, aes(x = k, y = deltaK)) +
  geom_line() +
  geom_point() +
  theme_minimal(base_size = 12) +
  labs(title = "gtsB: Proportional increase Δ(K) in A(K)",
       x = "K", y = "Δ(K)")
ggsave(file.path(out_dir, "monti_deltaKgtsB_4fold_.png"), pD, width = 6.5, height = 4.5, dpi = 600)

cdf_df <- bind_rows(cdf_rows)
pCDF <- ggplot(cdf_df, aes(x = x, y = F, color = factor(k))) +
  geom_line(linewidth = 1) +
  theme_minimal(base_size = 12) +
  labs(title = "gtsB: Consensus CDF curves by K",
       x = "Consensus value", y = "Empirical CDF",
       color = "K")
ggsave(file.path(out_dir, "monti_CDF_curvesgtsB_4fold_.png"), pCDF, width = 6.5, height = 4.5, dpi = 600)



