
library(tidyverse)
library(here)

#data
input_data1 <- read.csv("PP_distrib_per1000nulls_profiles_gtsB.csv", header = TRUE)

positions_2_3_4fold<-c(11,12,14,15,17,18,29,30,35,36,44,45,62,63,71,72,77,78,80,81,95,96,101,102,104,105,113,114,128,129,134,135,146,147,149,150,155,156,164,165,179,180,182,183,188,189,227,228,230,231,245,246,254,255,257,258,269,270,275,276,281,282,287,288,290,291,293,294,299,300,332,333,344,345,359,360,365,366,380,381,383,384,386,387,389,390,392,393,413,414,416,417,422,423,449,450,458,459,485,486,494,495,497,498,509,510,515,516,518,519,521,522,530,531,536,537,548,549,560,561,563,564,572,573,575,576,596,597,599,600,611,612,614,615,623,624,641,642,644,645,650,651,662,663,665,666,677,678,680,681,692,693,701,702,722,723,725,726,728,729,734,735,737,738,740,741,743,744,746,747,749,750,767,768,770,771,791,792,803,804,812,813,818,819,824,825,830,831,845,846,848,849,857,858,866,867,869,870,893,894
)

input_data<-input_data1%>%
    filter(nucleotide_position %in% positions_2_3_4fold)    
    

rate_col    <- "substitution_rate"
profile_cols <- c("A_C", "A_G", "A_T",
                  "C_A", "C_G", "C_T",
                  "G_A", "G_C", "G_T",
                  "T_A", "T_C", "T_G",
                  "transitions", "transversions")

#All data columns that should be together as a unit
data_cols <- c(rate_col, profile_cols)


#Shuffling function


shuffle_single_null <- function(null_df, 
                                rate_col,
                                profile_cols,
                                seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  #Separate position 2 and position 3
  pos2 <- null_df %>% filter(position == 2)
  pos3 <- null_df %>% filter(position == 3)
  
  #checking is working
  if (nrow(pos2) != nrow(pos3)) {
    warning(paste("Unequal position 2 and 3 sites:",
                  nrow(pos2), "vs", nrow(pos3)))
  }
  
  # Verify all 12 profile columns are present
  missing_cols <- setdiff(profile_cols, names(null_df))
  if (length(missing_cols) > 0) {
    stop(paste("Missing profile columns:", 
               paste(missing_cols, collapse = ", ")))
  }
  
  #Pool position 2 and 3 together
  #Rate and all 12 profile columns travel together as a unit
  #Position label is dropped from the pool
  pooled <- bind_rows(pos2, pos3) %>%
    select(null_id, codon_id, 
           all_of(rate_col), 
           all_of(profile_cols))  # 13 data columns travel together
  
  n_total <- nrow(pooled)
  n_each  <- nrow(pos2)
  
  # Randomly assign new position labels
  shuffled_labels <- sample(
    c(rep(2, n_each), rep(3, n_each)),
    size    = n_total,
    replace = FALSE
  )
  
  # Add shuffled position labels
  pooled$shuffled_position <- shuffled_labels
  
  return(pooled)
}

#shuffling across all 1000 nulls

null_ids        <- unique(input_data$null_id)
shuffled_results <- vector("list", length(null_ids))

for (i in seq_along(null_ids)) {
  
  null_subset <- input_data %>%
    filter(null_id == null_ids[i])
  
  shuffled_results[[i]] <- shuffle_single_null(
    null_df      = null_subset,
    rate_col     = rate_col,
    profile_cols = profile_cols,
    seed         = i * 42
  )
  
  if (i %% 100 == 0) {
    cat("Processed null", i, "of", length(null_ids), "\n")
  }
}

shuffled_data <- bind_rows(shuffled_results)


#Saving output

write.csv(shuffled_data,
          "shuffled_null_rates_profiles_4-fold_2_3post_gtsB.csv",
          row.names = FALSE)

