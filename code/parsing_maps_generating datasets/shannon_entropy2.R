library(here)
library(tidyverse)
library(data.table)
library(stringi)
library(ggstatsplot)
library(ggplot2)
library(sandwich)
library(lmtest)
library(emmeans)
library(dplyr)
library(purrr)


#Data

gtsb_gtsD_2_4<-read.csv("gtsb_gtsD_all_2_4_fixed_entropy_final.csv")


gtsD1<-gtsb_gtsD_2_4%>%
  filter(gene=="gtsD")%>%
  mutate(
  site_type = case_when(
    order_codon == 3 & Cluster == 1 ~ "syn_constrained",
    order_codon == 3 & !is.na(Cluster) & Cluster != 1 ~ "syn_other",
    order_codon == 2 ~ "nonsyn_2nd",
##    order_codon == 2 & Cluster == 1 ~ "nonsyn_2nd",
    TRUE ~ NA_character_
  ),
  site_type = factor(site_type, levels = c("syn_other", "syn_constrained", "nonsyn_2nd"))
) %>%
  filter(!is.na(site_type), !is.na(Normalized_entropy))

m2 <- lm(Normalized_entropy ~ site_type, data = gtsD1)
summary(m2)


coeftest(m2, vcov. = vcovHC(m2, type = "HC3"))


emm2 <- emmeans(m2, ~ site_type)
pairs(emm2, adjust = "holm")

confint(pairs(emm2, adjust = "holm"))

gtsb1<-gtsb_gtsD_2_4%>%
  filter(gene=="gtsB")%>%
  mutate(
    site_type = case_when(
      order_codon == 3 & Cluster == 1 ~ "syn_constrained",
      order_codon == 3 & !is.na(Cluster) & Cluster != 1 ~ "syn_other",
##      order_codon == 2 ~ "nonsyn_2nd",
          order_codon == 2 & Cluster == 1 ~ "nonsyn_2nd",
      TRUE ~ NA_character_
    ),
    site_type = factor(site_type, levels = c("syn_other", "syn_constrained", "nonsyn_2nd"))
  ) %>%
  filter(!is.na(site_type), !is.na(Normalized_entropy))

m1 <- lm(Normalized_entropy ~ site_type, data = gtsb1)
summary(m1)


coeftest(m1, vcov. = vcovHC(m1, type = "HC3"))

emm1 <- emmeans(m1, ~ site_type)
pairs(emm1, adjust = "holm")


confint(pairs(emm1, adjust = "holm"))

