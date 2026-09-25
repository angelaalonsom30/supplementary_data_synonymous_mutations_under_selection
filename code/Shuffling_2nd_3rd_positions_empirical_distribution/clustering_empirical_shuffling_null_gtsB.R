library(tidyverse)
library(readxl)
library(FactoMineR)  
library(cluster)   


observed_data <- read_excel("all_observed_4fold_2_3_gtsB_clustervariable.xlsx")

shuffled_data <- read.csv("shuffled_null_rates_profiles_4-fold_2_3post_gtsB.csv",
                          header = TRUE)

rate_col     <- "substitution_rate"
profile_cols <- c("A_C", "A_G", "A_T",
                  "C_A", "C_G", "C_T",
                  "G_A", "G_C", "G_T",
                  "T_A", "T_C", "T_G")
ts_col       <- "transitions"
tv_col       <- "transversions"

all_15_cols  <- c(rate_col, profile_cols, ts_col, tv_col)

#PCA on observed data


obs_matrix <- observed_data %>%
  select(all_of(all_15_cols)) %>%
  as.data.frame()

obs_pca <- PCA(
  obs_matrix,
  scale.unit = TRUE,
  ncp        = 2,
  graph      = FALSE
)

#Extract parameters for projecting null data
#Use svd$V for raw eigenvectors — critical for correct projection
pca_eigenvectors <- obs_pca$svd$V
pca_center       <- obs_pca$call$centre
pca_scale        <- obs_pca$call$ecart.type

cat("Variance explained by first 2 PCs:\n")
print(obs_pca$eig[1:2, ])

#HCPC on observed data with fixed k

optimal_k <- 5  # consensus clustering result

obs_hcpc <- HCPC(
  obs_pca,
  nb.clust = optimal_k,
  graph    = FALSE,
  consol   = FALSE
)

#Add cluster assignments to observed data
observed_data$cluster <- obs_hcpc$data.clust$clust

#Identify low substitution cluster
cluster_means <- observed_data %>%
  group_by(cluster) %>%
  summarise(
    mean_rate = mean(.data[[rate_col]]),
    n_pos2    = sum(position == 2),
    n_pos3    = sum(position == 3),
    n_total   = n(),
    .groups   = "drop"
  ) %>%
  arrange(mean_rate)


#Set low substitution cluster explicitly
# based on observed HCPC interpretation
low_sub_cluster <- 1
cat("\nLow substitution cluster set to:", 
    low_sub_cluster, "\n")

#Verify cluster 1 indeed has the lowest mean rate
cat("Mean rate of cluster 1:", 
    cluster_means %>% 
      filter(cluster == 1) %>% 
      pull(mean_rate), "\n")

#Observed constrained position 3 sites
obs_pos3_constrained <- observed_data %>%
  filter(position == 3,
         cluster  == low_sub_cluster)

cat("Observed constrained position 3 sites:",
    nrow(obs_pos3_constrained), "\n")

#clustering function using hclust

cluster_null_robust <- function(null_df,
                                all_15_cols,
                                pca_center,
                                pca_scale,
                                pca_eigenvectors,
                                optimal_k,
                                rate_col,
                                n_components = 2) {
  
  tryCatch({
    
    #Extract and scale null data using observed PCA parameters
    null_matrix <- null_df %>%
      select(all_of(all_15_cols)) %>%
      as.matrix()
    
    #Center and scale using observed parameters
    null_scaled <- sweep(null_matrix, 2, pca_center, "-")
    null_scaled <- sweep(null_scaled, 2, pca_scale,  "/")
    
    #Project onto observed PCA space using raw eigenvectors
    null_pca_scores <- null_scaled %*%
      pca_eigenvectors[, 1:n_components]
    
    null_scores_df <- as.data.frame(null_pca_scores)
    colnames(null_scores_df) <- paste0("Dim.", 1:n_components)
    
    #Hierarchical clustering using ward.D2
    #Same linkage as before
    null_dist    <- dist(null_scores_df, method = "euclidean")
    null_hclust  <- hclust(null_dist,    method = "ward.D2")
    
    #Cut tree at optimal k from consensus clustering
    null_clusters <- cutree(null_hclust, k = optimal_k)
    
    #Add cluster assignments to null data
    null_df$null_cluster <- null_clusters
    
    #Identify low substitution cluster in this null
    null_cluster_means <- null_df %>%
      group_by(null_cluster) %>%
      summarise(
        mean_rate = mean(.data[[rate_col]]),
        n_pos2    = sum(shuffled_position == 2),
        n_pos3    = sum(shuffled_position == 3),
        .groups   = "drop"
      ) %>%
      arrange(mean_rate)
    
    null_low_cluster <- null_cluster_means$null_cluster[1]
    
    #Count shuffled position 3 sites in low cluster
    n_constrained <- sum(
      null_df$shuffled_position == 3 &
        null_df$null_cluster      == null_low_cluster
    )
    
    return(list(
      success       = TRUE,
      n_constrained = n_constrained,
      site_data     = null_df
    ))
    
  }, error = function(e) {
    warning(paste("Clustering failed:", e$message))
    return(list(
      success       = FALSE,
      n_constrained = NA,
      site_data     = null_df
    ))
  })
}


#Apply to all 1000 shuffled nulls


null_ids       <- unique(shuffled_data$null_id)
null_results   <- numeric(length(null_ids))
null_site_data <- vector("list", length(null_ids))
failed_nulls   <- c()

for (i in seq_along(null_ids)) {
  
  null_subset <- shuffled_data %>%
    filter(null_id == null_ids[i])
  
  result <- cluster_null_robust(
    null_df          = null_subset,
    all_15_cols      = all_15_cols,
    pca_center       = pca_center,
    pca_scale        = pca_scale,
    pca_eigenvectors = pca_eigenvectors,
    optimal_k        = optimal_k,
    rate_col         = rate_col,
    n_components     = 2
  )
  
  if (result$success) {
    null_results[i]     <- result$n_constrained
    null_site_data[[i]] <- result$site_data
  } else {
    null_results[i] <- NA
    failed_nulls    <- c(failed_nulls, null_ids[i])
  }
  
  if (i %% 100 == 0) {
    cat("Processed null", i, "of", length(null_ids), "\n")
  }
}

#Check for failed nulls
if (length(failed_nulls) > 0) {
  cat("\nWarning:", length(failed_nulls),
      "nulls failed clustering:\n")
  print(failed_nulls)
}

#Remove failed nulls
null_results_clean <- null_results[!is.na(null_results)]
cat("\nSuccessful nulls:", length(null_results_clean),
    "of", length(null_ids), "\n")

#Comparison between observed vs null distribution


obs_n_constrained <- nrow(obs_pos3_constrained)

empirical_upper_PPP <- mean(null_results_clean >= obs_n_constrained)
empirical_lower_PPP <- mean(null_results_clean <= obs_n_constrained)
z_score <- (obs_n_constrained - mean(null_results_clean)) /
  sd(null_results_clean)

empirical_PPP_lower_corrected <- (sum(null_results_clean <= obs_n_constrained) + 1) /
  (length(null_results_clean) + 1)

empirical_PPP_upper_corrected <- (sum(null_results_clean >= obs_n_constrained) + 1) /
  (length(null_results_clean) + 1)



#Site level PP

#Pre-compute null_low for each null replicate once
null_low_lookup <- map_dbl(
  null_site_data[!sapply(null_site_data, is.null)],
  function(nd) {
    nd %>%
      group_by(null_cluster) %>%
      summarise(mr = mean(.data[[rate_col]]),
                .groups = "drop") %>%
      arrange(mr) %>%
      pull(null_cluster) %>%
      first()
  }
)

#Valid null data and lookup
valid_null_data   <- null_site_data[!sapply(null_site_data, is.null)]
valid_null_lookup <- null_low_lookup

# Site level PP with +1 correction
site_level_PPP <- obs_pos3_constrained %>%
  select(nucleotide_position) %>%
  mutate(
    site_PPP = map_dbl(nucleotide_position, function(pos) {
      
      site_null_counts <- map2_dbl(
        valid_null_data,
        valid_null_lookup,
        function(nd, null_low) {
          as.numeric(
            any(nd$nucleotide_position == pos &
                  nd$shuffled_position   == 3   &
                  nd$null_cluster        == null_low)
          )
        }
      )
      
      # Apply +1 correction
      n_nulls <- length(site_null_counts)
      (sum(site_null_counts) + 1) / (n_nulls + 1)
    })
  ) %>%
  arrange(site_PPP)

#Saving outputs


write.csv(site_level_PPP,
          "site_level_PPP_results_2_3_4fold_gtsB.csv",
          row.names = FALSE)

write.csv(
  data.frame(
    null_id       = null_ids[!is.na(null_results)],
    n_constrained = null_results_clean
  ),
  "null_distribution_summary_2_3_4fold_gtsB.csv",
  row.names = FALSE
)

#Plots

ggplot(
  data.frame(n_constrained = null_results_clean),
  aes(x = n_constrained)
) +
  geom_histogram(
    binwidth = 1,
    fill     = "steelblue",
    color    = "white",
    alpha    = 0.8
  ) +
  geom_vline(
    xintercept = obs_n_constrained,
    color      = "red",
    linewidth  = 1.2,
    linetype   = "dashed"
  ) +
  annotate(
    "text",
    x     = obs_n_constrained + 0.3,
    y     = Inf,
    label = paste("Observed =", obs_n_constrained),
    color = "red",
    hjust = 0,
    vjust = 2,
    size  = 4
  ) +
  labs(
    title    = paste("gtsB Null distribution — hclust ward.D2,",
                     "k =", optimal_k),
    subtitle = paste("lower-tail empirical p-value =",
                     round(empirical_PPP_lower_corrected, 3),
                     "| Z score =",
                     round(z_score, 3)),
    x        = "Number of third codon positions in constrained cluster",
    y        = "Count across 1100 null replicates"
  ) +
  theme_bw()

ggsave("null_distribution_hclust_gtsB.pdf",
       width = 8, height = 6)
ggsave("null_distribution_hclust_gtsB.png",
       width = 8, height = 6,dpi=1000 )



