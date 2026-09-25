library(here)
library(tidyverse)

#Reading the csv file with the mappings parsed from phylobayes
df <- read.csv("observed_predicted_posterior_mappings_gtsB_positions.csv", header = TRUE)


df <- df %>%
  mutate(original_row_index = row_number())

#Calculate the total number of mutations for each row
df$substitution_rate <- rowSums(df[, c(1:12)])

PP_dataframe<-df%>%
  filter(
    original_row_index  %% 2 == 0
  )%>%
  mutate(nucleotide_position= as.numeric(str_extract(name, "(?<=_)[0-9]+")),
  )%>%
  arrange(nucleotide_position)


n_replicates  <- 1100   # number of nulls per site
n_codons      <- nrow(PP_dataframe) / (3 * n_replicates)  

#checking id working
cat("Total rows:", nrow(df), "\n")
cat("Expected rows:", n_codons * 3 * n_replicates, "\n")
cat("Number of codons:", n_codons, "\n")

#Assign null_id, codon_id and position
PP_dataframef <- PP_dataframe %>%
  mutate(
    # null replicate ID cycles 1 to 1000 repeatedly
    null_id  = rep(1:n_replicates, 
                   times = n_codons * 3),
    
    # position cycles 1, 2, 3 — each repeated n_replicates times
    position = rep(rep(1:3, each = n_replicates), 
                   times = n_codons),
    
    # codon ID — each codon has 3 * n_replicates rows
    codon_id = rep(1:n_codons, 
                   each = 3 * n_replicates)
  )

#Checking
cat("\nCodon ID counts:\n")
print(table(PP_dataframef$codon_id) %>% head())
#Each codon should have 3000 or more rows (3 positions x 1000 replicates)

cat("\nNull ID counts:\n")
print(table(PP_dataframef$null_id) %>% head())
#Each null should have n_codons * 3 rows


PP_dataframef$transitions <- rowSums(PP_dataframef[, c("A_G","G_A","C_T","T_C")])


PP_dataframef$transversions <- rowSums(PP_dataframef[, c("A_C","C_A","A_T","T_A","T_G","G_T","C_G","G_C" )])  


PP_dataframef$total_type<- rowSums(PP_dataframef[, c("A_G","G_A","C_T","T_C","A_C","C_A","A_T","T_A","T_G","G_T","C_G","G_C" )])  


write_csv(PP_dataframef, here("PP_distrib_per1000nulls_profiles_gtsB.csv"), append = FALSE)

  
#Group data by nucleotide position replicate and index even number as observed and uneven as predicted, compared them while maintaining order
pair_comparison <- df %>%
  group_by(name) %>% ##name here has all the position and replicate information
  arrange(original_row_index)%>% # To ensure rows are in their original order
  mutate(pair_group = (row_number() + 1) %/% 2) %>% # Create a group for each pair
  group_by(name, pair_group) %>%
  summarise(
    row1_index = original_row_index[1],
    row2_index = original_row_index[2],
    row1_substitutions = substitution_rate[1],
    row2_substitutions = substitution_rate[2],
    pp_high1 = as.integer(substitution_rate[2] >= substitution_rate[1]),
    pp_low1 = as.integer(substitution_rate[2] <= substitution_rate[1]),  
    .groups = "drop"
  )

comparison_summary <- pair_comparison %>%
  group_by(name) %>% 
  summarise(
    count_pp_1_observed = sum(pp_high1 == 1, na.rm = TRUE),
    total_pairs = n(), # Total number of pairs in the category
    pp_final = sum(pp_high1 == 1, na.rm = TRUE) / total_pairs,
    PP_high = mean(pp_high1),
    PP_low  = mean(pp_low1),
    PP_two  = min(1, 2 * min(PP_high, PP_low)),
    substitution_rate = mean(row1_substitutions, na.rm = TRUE),
    sd_observed_substitution_rate = sd(row1_substitutions, na.rm = TRUE),
    predicted_substitution_rate = mean(row2_substitutions, na.rm = TRUE),
    sd_predicted_substitution_rate = sd(row2_substitutions, na.rm = TRUE),
    median_predictive = median(row2_substitutions, na.rm = TRUE),  
    Effect_size = abs(substitution_rate - median_predictive) / sd_predicted_substitution_rate
  )

comparison_summary_plus_zscore<-comparison_summary%>%
  mutate(
    z_score=(substitution_rate - predicted_substitution_rate)/sd_predicted_substitution_rate,
    pval=(1-pnorm(abs(z_score))),
    nucleotide_position= as.numeric(str_extract(name, "(?<=_)[0-9]+"))
  )%>%
  arrange(nucleotide_position) %>% # Order by position from lowest to highest
  mutate(order_codon = rep(1:3, length.out = n()),
         codon_number = ceiling(row_number() / 3) 
  )
write_csv(comparison_summary_plus_zscore, here("all_twosidedPP_gtsB_effect_size.csv"), append = FALSE)


###
fold_structure<-read.csv("to_merge_gtsB_fold_codon_info_structure_gtsB.csv") ##this file contains the topology prediction at each codon position


merge_total_gtsb_fold_structure<-comparison_summary_plus_zscore%>%
  full_join(fold_structure, by = c("nucleotide_position",	"order_codon"	,"codon_number"
))%>%
  select(c(11:20,4:9))

####write_csv(join_predicti_obser_4_fold_completephylo, here("join_predicti_obser_4_fold_completephylo.csv"), append = FALSE)


###write_csv(pairs_less_1100, here("pairs_less_1100.csv"), append = FALSE)

#pertype (substitution profile) _separate first (observed) and second rows (posterior predictive) into their own dataframes, maintaining order and including row index
observed_first_rows <- df %>%
  filter(original_row_index %% 2 != 0) %>% #observed or posterior mappings
  select(original_row_index, everything())%>%
  mutate(nucleotide_position= as.numeric(str_extract(name, "(?<=_)[0-9]+"))
  )%>%
   group_by(name,nucleotide_position)%>%
  summarise(across(
  .cols = A_T:G_C, 
  .fns = list(avg = mean), na.rm = TRUE, 
  .names = "{col}_{fn}"
   ))
observed_second_rows <- df %>%  ##these are the posterior predictive mappings
  filter(original_row_index %% 2 == 0)%>% ##predic
  select(original_row_index, everything())%>%
  mutate(nucleotide_position= as.numeric(str_extract(name, "(?<=_)[0-9]+"))
  )%>%
  group_by(name,nucleotide_position)%>%
  summarise(across(
    .cols = A_T:G_C, 
    .fns = list(avg = mean), na.rm = TRUE, 
    .names = "{col}_{fn}"
  ))


observed_first_rows$transitions_avg <- rowSums(observed_first_rows[, c("A_G_avg","G_A_avg","C_T_avg","T_C_avg")])


observed_first_rows$transversions_avg <- rowSums(observed_first_rows[, c("A_C_avg","C_A_avg","A_T_avg","T_A_avg","T_G_avg","G_T_avg","C_G_avg","G_C_avg" )])  


observed_first_rows$total_type_avg<- rowSums(observed_first_rows[, c("A_G_avg","G_A_avg","C_T_avg","T_C_avg","A_C_avg","C_A_avg","A_T_avg","T_A_avg","T_G_avg","G_T_avg","C_G_avg","G_C_avg" )])  

observed_second_rows$transitions_avg <- rowSums(observed_second_rows[, c("A_G_avg","G_A_avg","C_T_avg","T_C_avg")])


observed_second_rows$transversions_avg <- rowSums(observed_second_rows[, c("A_C_avg","C_A_avg","A_T_avg","T_A_avg","T_G_avg","G_T_avg","C_G_avg","G_C_avg" )])  


observed_second_rows$total_type_avg<- rowSums(observed_second_rows[, c("A_G_avg","G_A_avg","C_T_avg","T_C_avg","A_C_avg","C_A_avg","A_T_avg","T_A_avg","T_G_avg","G_T_avg","C_G_avg","G_C_avg" )])  


merged_first_second<-observed_first_rows%>%
  full_join(observed_second_rows, by = c("nucleotide_position")
  )
write_csv(merged_first_second, here("gtsB_allfold_totalandpertype_obser_predicti.csv"), append = FALSE)


write_csv(observed_first_rows, here("gtsB_allfold_posterio_predicitve_totalandpertype.csv"), append = FALSE)


##Merging total counts, per type_fold

merge_total_gtsb_fold_structure_per_type<-merge_total_gtsb_fold_structure%>%
  full_join(observed_first_rows, by = c("nucleotide_position")
  )%>%
  select(c(-17,-32))


write_csv(merge_total_gtsb_fold_structure_per_type, here("final_merge_total_gtsb_fold_structure_per_type.csv"), append = FALSE)

predicted_second_rows <- df %>%
  filter(original_row_index %% 2 == 0) %>%
  select(original_row_index, everything()) # Retain original row index


###Now here I have the codon combination output from the python script file, converting it from wide to long this is only for 4-fold sites
my_data_sin_ancestro_observ_gtsB<-read.csv('ready_codon_gtsB_codones_sin_ancestros.csv')


my_data_sin_ancestro_observ_gtsB_long<-pivot_longer(my_data_sin_ancestro_observ_gtsB, cols=2:65, names_to = "codon_aa", values_to = "probability")%>%
  filter(!str_detect(codon_observ,regex("Stop", ignore_case = TRUE)))

##write_csv(my_data_sin_ancestro_observ_gtsB_long, here("my_data_sin_ancestro_observ_gtsB_long.csv"), append = FALSE)

my_data_sin_ancestro_observ_gtsB_long2<-my_data_sin_ancestro_observ_gtsB_long%>%
 mutate(aa_sbw25 = str_extract(codon_observ, "(?<=\\.\\.\\.)(.)"))

##Merging with the big dataframe

complete_total_pertype_structure_codon_gtsB<-merge_total_gtsb_fold_structure_per_type%>%
  full_join(my_data_sin_ancestro_observ_gtsB_long2, by = c("codon_number","aa_sbw25"))%>%
  filter(!nucleotide_position=='NA')%>%
  filter(!is.na(mean_freq))%>%
  pivot_wider(
    names_from = codon_aa,
    values_from = probability)


  write_csv(complete_total_pertype_structure_codon_gtsB, here("final_complete_total_pertype_structure_codon_gtsB.csv"), append = FALSE)


#Filtering by fold type and aa to just have 4-fold sites present in the reference genome P.fluorescens SBW25


complete_total_pertype_structure_codon_gtsB_4_v <- complete_total_pertype_structure_codon_gtsB %>%
  filter(fold==4)%>%
  filter(aa_sbw25=="V" | aa_sbw25=="A" | aa_sbw25=="G" |aa_sbw25=="T" |aa_sbw25=="P")%>%
  select(where(~ !all(is.na(.))))

df<-complete_total_pertype_structure_codon_gtsB_4_v        

complete_total_pertype_structure_codon_gtsB_4_vf<-complete_total_pertype_structure_codon_gtsB_4_v%>%
  rename_with(~ c("T-ending_codon_probability",	"C-ending codon probability",	"A-ending codon probability",	"G-ending codon probability"),
              .cols = tail(names(df), 4))



##Now merging all 4-fold dataframes

complete_4_fold_gtsB_codon<- bind_rows(complete_total_pertype_structure_codon_gtsB_4_vf, complete_total_pertype_structure_codon_gtsB_4_af, complete_total_pertype_structure_codon_gtsB_4_gf, complete_total_pertype_structure_codon_gtsB_4_pf, complete_total_pertype_structure_codon_gtsB_4_tf)
  
write_csv(complete_4_fold_gtsB_codon, here("final_complete_4_fold_gtsB_codon_ready.csv"), append = FALSE)


###Now codon Relative probability for each codon ending probability of 5 aminoacids with codons with a 4-fold site (20 codons)


relative_frq_4fold<-complete_4_fold_gtsB_codon%>%
    group_by(nucleotide_position)%>%
  mutate(
    total_sum_codon_fre= T-ending_codon_probability +	C-ending_codon_probability +	A-ending_codon_probability +	G-ending_codon_probability,
    T_ending_codon_relative_probability = T-ending_codon_probability/total_sum_codon_fre,
    A_ending_codon_relative_probability = A-ending_codon_probability/total_sum_codon_fre,
    G_ending_codon_relative_probability = G-ending_codon_probability/total_sum_codon_fre,
    C_ending_codon_relative_probability = C-ending_codon_probability/total_sum_codon_fre,
    total_relati= relativef_T_ending + relativef_A_ending + relativef_G_ending + relativef_C_ending
    ) %>%
  ungroup()

###Merging with CAI-w and RSCU, plus region 4-fold

gtsB_codonusage_region_4fold<-read.csv("gtsD_gtsB_to_merge_region_rscu_wi_4fold.csv")%>%
  filter(gene=="gtsB")

gtsB_4fold_completo<-relative_frq_4fold%>%
  full_join(gtsB_codonusage_region_4fold, by = c("nucleotide_position","codon_number"))

write_csv(gtsB_4fold_completo, here("gtsB_4fold_completo.csv"), append = FALSE)


##2-4 neighboring codon position and 4-fold site comparison


fold_2_4_gtsB<-merge_total_gtsb_fold_structure_per_type%>%
  filter(aa_sbw25=="A" | aa_sbw25=="G" | aa_sbw25=="T" | aa_sbw25=="P" | aa_sbw25=="V")%>%
  filter(order_codon==2 | order_codon==3)

write_csv(fold_2_4_gtsB, here("fold_2_4_gtsB.csv"), append = FALSE)

