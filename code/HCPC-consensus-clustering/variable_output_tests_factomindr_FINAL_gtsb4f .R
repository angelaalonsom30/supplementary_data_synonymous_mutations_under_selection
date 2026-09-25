##After evaluating the best cluster partition based on the consensus clustering metrics, comes the analysis of the variables on
#the selected cluster partition

library(FactoMineR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(readxl)

df <- read_excel("4_fold_sites_in_gtsB_clusterv.xlsx") #this data frame has all the active and supplementary quantitative and qualitative 
###variables per each site, but includes the cluster variable which was determined via consensus metrics

df_all<-df
#cluster (factor)
#quantitative vars (active + supplementary)
#categorical vars (supplementary, as factors)

out_dir <- "out_gtsB_consensus_final_final_4fold_variablesf"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

get_var_types <- function(df_test, cluster_col = "cluster") {
  stopifnot(cluster_col %in% names(df_test))
  num_vars <- names(df_test)[sapply(df_test, is.numeric)]
  num_vars <- setdiff(num_vars, cluster_col)
  
  fac_vars <- names(df_test)[sapply(df_test, is.factor)]
  fac_vars <- setdiff(fac_vars, cluster_col)
  
  list(num_vars = num_vars, fac_vars = fac_vars)
}

#Global contributions of each variable to the cluster construction

#Quantitative variables: ANOVA + eta^2 + p-value per variable
compute_quanti_link <- function(df_test, cluster_col = "cluster") {
  y <- factor(df_test[[cluster_col]])
  num_vars <- names(df_test)[sapply(df_test, is.numeric)]
  num_vars <- setdiff(num_vars, cluster_col)
  if (length(num_vars) == 0) return(tibble())
  
  bind_rows(lapply(num_vars, function(v) {
    dat <- df_test[, c(cluster_col, v)]
    dat <- dat[complete.cases(dat), , drop = FALSE]
    dat[[cluster_col]] <- factor(dat[[cluster_col]])
    if (nrow(dat) < 3) return(NULL)
    
    fit <- stats::aov(dat[[v]] ~ dat[[cluster_col]])
    sm  <- summary(fit)[[1]]
    
    ss_between <- sm[1, "Sum Sq"]
    ss_within  <- sm[2, "Sum Sq"]
    ss_total   <- ss_between + ss_within
    eta2 <- if (ss_total > 0) as.numeric(ss_between / ss_total) else NA_real_
    
    tibble(
      variable = v,
      eta2 = eta2,
      F = as.numeric(sm[1, "F value"]),
      df1 = as.integer(sm[1, "Df"]),
      df2 = as.integer(sm[2, "Df"]),
      p_value = as.numeric(sm[1, "Pr(>F)"])
    )
  }))
}

#Categorical variables: chi-square test per factor variable (cluster x variable)
compute_categorical_link <- function(df_test, cluster_col = "cluster") {
  fac_vars <- names(df_test)[sapply(df_test, is.factor)]
  fac_vars <- setdiff(fac_vars, cluster_col)
  if (length(fac_vars) == 0) return(tibble())
  
  bind_rows(lapply(fac_vars, function(v) {
    tab <- table(df_test[[cluster_col]], df_test[[v]])
    if (all(dim(tab) >= 2)) {
      tst <- suppressWarnings(stats::chisq.test(tab))
      tibble(
        variable = v,
        chisq = unname(tst$statistic),
        df = unname(tst$parameter),
        p_value = unname(tst$p.value)
      )
    } else {
      tibble(variable = v, chisq = NA_real_, df = NA_real_, p_value = NA_real_)
    }
  }))
}

#Contribution by cluster

#Quantitative: v-test + means/sds per cluster per variable like factormineR
compute_quanti_description <- function(df_test, cluster_col = "cluster") {
  df_test[[cluster_col]] <- factor(df_test[[cluster_col]])
  
  num_vars <- names(df_test)[sapply(df_test, is.numeric)]
  num_vars <- setdiff(num_vars, cluster_col)
  if (length(num_vars) == 0) return(tibble())
  
  N <- nrow(df_test)
  clusters <- levels(df_test[[cluster_col]])
  
  out <- lapply(num_vars, function(v) {
    x <- df_test[[v]]
    keep <- !is.na(x) & !is.na(df_test[[cluster_col]])
    x <- x[keep]
    g <- df_test[[cluster_col]][keep]
    N2 <- length(x)
    
    mu <- mean(x)
    sd_all <- stats::sd(x)
    
    bind_rows(lapply(clusters, function(cl) {
      idx <- which(g == cl)
      ncl <- length(idx)
      if (ncl == 0) return(NULL)
      
      mu_cl <- mean(x[idx])
      sd_cl <- if (ncl > 1) stats::sd(x[idx]) else NA_real_
      
      # Finite population correction (without replacement)
      # se = sd_all * sqrt( (N2 - ncl)/(N2 - 1) / ncl )
      se <- if (!is.na(sd_all) && N2 > 1 && ncl > 0) sd_all * sqrt(((N2 - ncl) / (N2 - 1)) / ncl) else NA_real_
      vtest <- if (!is.na(se) && se > 0) (mu_cl - mu) / se else NA_real_
      pval <- if (!is.na(vtest)) 2 * (1 - stats::pnorm(abs(vtest))) else NA_real_
      
      tibble(
        cluster = cl,
        variable = v,
        v.test = vtest,
        mean_in_cluster = mu_cl,
        overall_mean = mu,
        sd_in_cluster = sd_cl,
        overall_sd = sd_all,
        p.value = pval,
        n_in_cluster = ncl
      )
    }))
  })
  
  bind_rows(out)
}

#Categorical: Cla/Mod, Mod/Cla, Global, v.test, p.value for each level in each cluster
compute_categorical_levelwise <- function(df_test, cluster_col = "cluster") {
  df_test[[cluster_col]] <- factor(df_test[[cluster_col]])
  fac_vars <- names(df_test)[sapply(df_test, is.factor)]
  fac_vars <- setdiff(fac_vars, cluster_col)
  if (length(fac_vars) == 0) return(tibble())
  
  N <- nrow(df_test)
  clusters <- levels(df_test[[cluster_col]])
  
  bind_rows(lapply(fac_vars, function(v) {
    f <- df_test[[v]]
    g <- df_test[[cluster_col]]
    keep <- !is.na(f) & !is.na(g)
    f <- droplevels(f[keep])
    g <- droplevels(g[keep])
    
    N2 <- length(f)
    if (N2 == 0) return(NULL)
    
    # counts
    n_level <- table(f)            #level totals
    n_cluster <- table(g)          #cluster totals
    tab <- table(g, f)             #cluster x level counts
    
    bind_rows(lapply(clusters, function(cl) {
      if (!(cl %in% names(n_cluster))) return(NULL)
      ncl <- as.numeric(n_cluster[cl])
      
      bind_rows(lapply(levels(f), function(lv) {
        nlv <- as.numeric(n_level[lv])
        ncl_lv <- if (cl %in% rownames(tab) && lv %in% colnames(tab)) as.numeric(tab[cl, lv]) else 0
        
        # Definitions
        ClaMod <- if (nlv > 0) ncl_lv / nlv else NA_real_         #P(cluster/level)
        ModCla <- if (ncl > 0) ncl_lv / ncl else NA_real_         #P(level/cluster)
        Global <- if (N2 > 0) nlv / N2 else NA_real_              #P(level)
        
        # v-test for over/under-representation in cluster
        # v = (ModCla - Global) / sqrt(Global*(1-Global)/ncl)
        se <- if (!is.na(Global) && ncl > 0) sqrt(Global * (1 - Global) / ncl) else NA_real_
        vtest <- if (!is.na(se) && se > 0) (ModCla - Global) / se else NA_real_
        pval <- if (!is.na(vtest)) 2 * (1 - stats::pnorm(abs(vtest))) else NA_real_
        
        tibble(
          cluster = cl,
          variable = v,
          level = paste0(v, "=", lv),
          `Cla/Mod` = ClaMod,
          `Mod/Cla` = ModCla,
          Global = Global,
          p.value = pval,
          v.test = vtest,
          n_level = nlv,
          n_in_cluster = ncl,
          n_cluster_level = ncl_lv
        )
      }))
    }))
  }))
}


#functions and saving

run_full_variable_cluster_reports <- function(df_all,
                                              out_prefix,
                                              cluster_col = "cluster",
                                              id_cols = c("nucleotide_position"),
                                              force_factor_cols = NULL,
                                              convert_char_to_factor = TRUE,
                                              write_tsv = TRUE) {
  
  stopifnot(cluster_col %in% names(df_all))
  dir.create(dirname(out_prefix), showWarnings = FALSE, recursive = TRUE)
  
  df_test <- df_all
  
  #Drop ID columns (do not test them)
  drop_ids <- intersect(id_cols, names(df_test))
  if (length(drop_ids) > 0) df_test <- df_test %>% select(-all_of(drop_ids))
  
  #Ensure cluster is factor
  df_test[[cluster_col]] <- factor(df_test[[cluster_col]])
  
  #Optionally convert character -> factor
  if (isTRUE(convert_char_to_factor)) {
    char_cols <- names(df_test)[sapply(df_test, is.character)]
    char_cols <- setdiff(char_cols, cluster_col)
    if (length(char_cols) > 0) df_test[char_cols] <- lapply(df_test[char_cols], as.factor)
  }
  
  #Force specified categorical variables to factor
  if (!is.null(force_factor_cols)) {
    force_factor_cols <- intersect(force_factor_cols, names(df_test))
    df_test[force_factor_cols] <- lapply(df_test[force_factor_cols], as.factor)
  }
  
  #Drop constants
  is_constant <- sapply(df_test, function(x) {
    if (all(is.na(x))) return(TRUE)
    ux <- unique(x[!is.na(x)])
    length(ux) <= 1
  })
  if (any(is_constant)) df_test <- df_test[, !is_constant, drop = FALSE]
  
  #write
  writeLines(names(df_test), con = paste0(out_prefix, "_variables_tested.txt"))
  
  #Global results
  part1_quanti <- compute_quanti_link(df_test, cluster_col)
  part1_categ  <- compute_categorical_link(df_test, cluster_col)
  
  #per cluster
  part2_quanti <- compute_quanti_description(df_test, cluster_col)
  part2_categ  <- compute_categorical_levelwise(df_test, cluster_col)
  
  #export csv and txt
  write_both <- function(df, stem) {
    csv_path <- paste0(out_prefix, "_", stem, ".csv")
    write_csv(df, csv_path)
    if (isTRUE(write_tsv)) write_tsv(df, paste0(out_prefix, "_", stem, ".txt"))
  }
  
  if (nrow(part1_quanti) > 0) write_both(part1_quanti, "Global_quanti_link_eta2_F_df_pvalue")
  if (nrow(part1_categ)  > 0) write_both(part1_categ,  "Global_categorical_link_chisq_df_pvalue")
  if (nrow(part2_quanti) > 0) write_both(part2_quanti, "Per_cluster_quanti_description_vtests_means_sds_ALL")
  if (nrow(part2_categ)  > 0) write_both(part2_categ,  "Per_cluster_categorical_levelwise_ClaMod_ModCla_Global_vtest_ALL")
  
  invisible(list(
    df_test = df_test,
    part1_quanti = part1_quanti,
    part1_categ = part1_categ,
    part2_quanti = part2_quanti,
    part2_categ = part2_categ
  ))
}


#running analysis
out <- run_full_variable_cluster_reports(
   df_all,
   out_prefix = file.path(out_dir, "gtsB_consensus_K3"), 
   cluster_col = "cluster",
   id_cols = c("nucleotide_position"),
   force_factor_cols = c("amino_acid_P. fluorescens SBW25", "Topological structure", "region"),
   convert_char_to_factor = TRUE,
   write_tsv = TRUE
 )















