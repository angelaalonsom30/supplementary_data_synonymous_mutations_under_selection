##https://aosmith.rbind.io/2019/03/25/getting-started-with-emmeans/


library(here)
library(tidyverse)
library(data.table)
library(stringi)
library(stringr)
library(ggstatsplot)
library(ggplot2)
library(dplyr)
library(MASS) 
library(emmeans)
library(readxl)
library(openxlsx)

#Here loading all fold info gtsB and gtsD

allfold_gtsB_gtsD_model<-read_xlsx("all_fold_gtsB_gtsD.xlsx")

allfold_gtsB_gtsD_model$fold <- as.factor(allfold_gtsB_gtsD_model$fold)

allfold_gtsB_gtsD_model$gene<- as.factor(allfold_gtsB_gtsD_model$gene)


#rounding data, as nb needs integer
allfold_gtsB_gtsD_model$susbtitution_rate <- round(allfold_gtsB_gtsD_model$susbtitution_rate)

reduced_model <- glm.nb(susbtitution_rate ~gene+fold , data =allfold_gtsB_gtsD_model)
full_model <- glm.nb(susbtitution_rate ~gene*fold , data =allfold_gtsB_gtsD_model)


#Perform Likelihood Ratio Test
lrt <- anova(reduced_model, full_model, test = "Chisq")
print(lrt)

##Estimated Marginal Means
mymodelgtsB_gtsD_nbwe<-glm.nb(susbtitution_rate ~gene*fold , data =allfold_gtsB_gtsD_model)

summary(mymodelgtsB_gtsD_nbwe)

em_interaction <- emmeans(mymodelgtsB_gtsD_nbwe, ~ gene * fold)

pairwise_comparisons_fdr <- contrast(em_interaction, method = "pairwise", adjust = "fdr")

print(pairwise_comparisons_fdr)
































allfold_gtsB_gtsD_model$obser_total_mut_Mean_celi<-  ifelse(allfold_gtsB_gtsD_model$obser_total_mut_Mean %% 1 >= 0.5, ceiling(allfold_gtsB_gtsD_model$obser_total_mut_Mean), floor(allfold_gtsB_gtsD_model$obser_total_mut_Mean))
  
allfold_gtsB_gtsD_model<-allfold_gtsB_gtsD_model%>%
  filter(!fold=='ZERO' & !fold=='TWO'& !fold=='SIX' )
  

allfold_gtsB_gtsD_modelwe<-read.csv("gtsb_gtsd_allfold_try2_we.csv")


##histogram

mygraph <- ggplot(allfold_gtsB_gtsD_model, aes(x = obser_total_mut_Mean, xlab="substitutions"))

mygraph <- mygraph +
  # add data density smooth
  geom_density() +
  # add rug (bars at the bottom of the plot)
  geom_rug() +
  # add black semitransparent histogram
  geom_histogram(aes(y = ..density..),
                 color = "black",
                 alpha = 0.3) +
  # add normal curve in red, with mean and sd from fklength
  stat_function(fun = dnorm,
                args = list(
                  mean = mean(allfold_gtsB_gtsD_model$obser_total_mut_Mean),
                  sd = sd(allfold_gtsB_gtsD_model$obser_total_mut_Mean)
                ),
                color = "red")
#display graph, by location
mygraph+facet_grid(.~gene+fold+pos_codon)

### data table form (contingency table)

allfold_gtsB_gtsD_modelfwe<-allfold_gtsB_gtsD_modelwe%>%
  select(c(1,3,6))

allfold_gtsB_gtsD_modelff<-xtabs(obser_total_mut_Mean~fold+gene, data=allfold_gtsB_gtsD_modelf)

chisq.test(allfold_gtsB_gtsD_modelff) # runs chi square test of independence of gene and  fold
###model to test for differences between genes 
##https://www.biostars.org/p/9504330/

library(tidyr)
library(dplyr)
library(broom)

allfold_gtsB_gtsD_model$fold <- as.factor(allfold_gtsB_gtsD_model$fold)

allfold_gtsB_gtsD_model$gene<- as.factor(allfold_gtsB_gtsD_model$gene)

allfold_gtsB_gtsD_model$pos_codon<- as.factor(allfold_gtsB_gtsD_model$pos_codon)

glimpse(allfold_gtsB_gtsD_model)


"susbtitution_rate" %in% names(allfold_gtsB_gtsD_model)
mymodelgtsB_gtsD_poisson<-glm(susbtitution_rate~ fold ,family = poisson(link = "log"),control = list(maxit = 1000),data=allfold_gtsB_gtsD_model)

#mymodelgtsB_gtsD_poissonwe<-glm(weighted~fold*gene,family = poisson(link = "log"),control = list(maxit = 100),data=allfold_gtsB_gtsD_model)


summary(mymodelgtsB_gtsD_poisson)

mymodelgtsB_gtsD_quasipoisson<-glm(susbtitution_rate~gene+fold,family = quasipoisson(link = "log"),control = list(maxit = 1000),data=allfold_gtsB_gtsD_model)


summary(mymodelgtsB_gtsD_quasipoisson)


library(tidyr)
library(dplyr)
library(broom)

##negative binomial

library(MASS)
##https://aosmith.rbind.io/2019/03/25/getting-started-with-emmeans/

##no weight
mymodelgtsB_gtsD_nbwe<-glm.nb(susbtitution_rate ~gene*fold , data =allfold_gtsB_gtsD_model)

summary(mymodelgtsB_gtsD_nbwe)

emmwe <- emmeans(mymodelgtsB_gtsD_nbwe, ~ gene*fold, data = allfold_gtsB_gtsD_model)
#emmwe2 <- emmeans(mymodelgtsB_gtsD_nbwe, specs = pairwise ~ gene*fold, data = allfold_gtsB_gtsD_model,  type = "response")
#emmwe2$emmeans

x<-pairs(emmwe , adjust="fdr")

#emmwe$contrasts %>%
#  summary(infer = TRUE)


#str(regrid(x))

summary(x, infer = TRUE)

library(car)

Anova(mymodelgtsB_gtsD_nbwe, type=3, test = "LR")
## random effects gene size fix fold



pairs(emm1 , adjust="tukey")


##regression just to see if codon usage predicts subs rate in 4-fold

codon_usage_4_gtsb_gtsd<-read.csv('data_cluster_all4codons_gtsd_gtsb.csv')

codon_usage_4_gtsb_gtsd$gene<- as.factor(codon_usage_4_gtsb_gtsd$gene)







##permutation

# Load required libraries
library(MASS)  # For the glm.nb() function
library(dplyr) # For data manipulation

# Simulated data
set.seed(123)
n <- 100  # Sample size
x <- rnorm(n)  # Predictor variable
mu <- exp(0.5 + 0.8 * x)  # True mean parameter
size <- 1.5  # Dispersion parameter
y <- rnbinom(n, mu = mu, size = size)  # Response variable (Negative Binomial distributed)

# Create a dataframe
data <- data.frame(x = x, y = y)

# Fit Negative Binomial regression model
model <- glm.nb(y ~ x, data = data)

# Function to compute test statistic (e.g., coefficient of predictor variable)
compute_test_statistic <- function(mymodelgtsB_gtsD_nbwe) {
  coef(mymodelgtsB_gtsD_nbwe)["x"]  # Coefficient of predictor variable
}

# Observed test statistic
observed_statistic <- compute_test_statistic(mymodelgtsB_gtsD_nbwe)

# Number of permutations
n_permutations <- 10

# Permutation test
permuted_statistics <- replicate(n_permutations, {
  # Permute the response variable without replacement
  permuted_data <- mutate(allfold_gtsB_gtsD_model, y_permuted = sample(obser_total_mut_Mean, replace = FALSE))
  
  # Refit the model with permuted data
  permuted_model <- glm.nb(y_permuted ~ gene, data = permuted_data)
  
  # Compute test statistic for permuted data
  compute_test_statistic(permuted_model)
})

# Compute empirical p-value
empirical_p_value <- mean(permuted_statistics >= observed_statistic)

# Print results
cat("Observed Test Statistic:", observed_statistic, "\n")
cat("Empirical p-value:", empirical_p_value, "\n")




##install.packages("glmmTMB")
library(glmmTMB)
ggplot(allfold_gtsB_gtsD_modelfwe, aes(obser_total_mut_Mean)) + geom_histogram() + scale_x_log10()

# Assuming your data frame is called 'data' and the response variable is 'y'
model <- glmmTMB(obser_total_mut_Mean~fold*gene+(1+fold | gene), data = allfold_gtsB_gtsD_modelfwe, family = nbinom2)
summary(model)
+ (1 + pizza + time |subject)
allfold_gtsB_gtsD_modelfwe$obser_total_mut_Mean<-floor(allfold_gtsB_gtsD_modelfwe$obser_total_mut_Mean)

library(pscl)

summary(m1 <- zeroinfl(obser_total_mut_Mean~fold+ gene, data = allfold_gtsB_gtsD_modelfwe))
##weigted

allfold_gtsB_gtsD_modelfwe$weight <- 1 / ave(allfold_gtsB_gtsD_modelfwe$obser_total_mut_Mean, allfold_gtsB_gtsD_modelfwe$fold, allfold_gtsB_gtsD_modelfwe$gene, FUN = length)



allfold_gtsb_gtsdtry$fold <- as.factor(allfold_gtsb_gtsdtry$fold)

allfold_gtsb_gtsdtry$gene<- as.factor(allfold_gtsb_gtsdtry$gene)
allfold_gtsb_gtsdtry$gene_length<- as.factor(allfold_gtsb_gtsdtry$gene_length)


model_poisson<-glm(obser_total_mut_Mean ~gene*fold, data =allfold_gtsb_gtsdtry, family = poisson(link = "log"),control = list(maxit = 100))
summary(model_poisson)
emm <- emmeans(model_poisson, ~ gene, data = allfold_gtsb_gtsdtry)

pairs(emm , adjust="holm")


mymodel1<-glm.nb(obser_total_mut_Mean ~gene*fold , data =allfold_gtsb_gtsdtry )
summary(mymodel1)

emm1 <- emmeans(mymodel1, ~gene*fold, data = allfold_gtsb_gtsdtry)

pairs(emm1 , adjust="tukey")

coe<-exp(-0.14758)


####weighted model

allfold_gtsb_gtsdtry$weight <- 1 / ave(allfold_gtsb_gtsdtry$obser_total_mut_Mean, allfold_gtsb_gtsdtry$fold, FUN = length)

# Fit weighted negative binomial regression model
model_weighted_nb <- glm.nb( obser_total_mut_Mean ~gene*fold, 
                            data = allfold_gtsb_gtsdtry, weights = weight)

# Check overall model significance
summary(model_weighted_nb)

emm1 <- emmeans(model_weighted_nb , ~gene*fold, data = allfold_gtsb_gtsdtry)

pairs(emm1 , adjust="tukey")







# Modify your code to use glm.nb for negative binomial regression
per_gene_models <- allfold_gtsb_gtsd %>%
  group_by(gene) %>%
  do(tidy(glm.nb(obser_total_mut_Mean ~ fold, data = .)))

# Collect the results into a single data frame
per_gene_models <- collect(per_gene_models)

# Convert "fold" into an unordered factor


# Now you can use relevel to change the reference category
allfold_gtsb_gtsd$fold <- relevel(allfold_gtsb_gtsd$fold, ref = "ZERO")

# Fit the negative binomial regression model with the new reference
library(MASS)
per_gene_models <- allfold_gtsb_gtsd %>%
  group_by(gene) %>%
  do(tidy(glm.nb(obser_total_mut_Mean ~ fold, data = .)))

# Collect the results into a single data frame
per_gene_models <- collect(per_gene_models)

type_diff <- per_gene_models %>%
  mutate(padj = p.adjust(p.value, method="BH"))


##accounting for gene length

library(dplyr)

allfold_gtsb_gtsd2 <- allfold_gtsb_gtsd %>%
  group_by(gene) %>%
  mutate(gene_length = n())

allfold_gtsb_gtsd2f<-allfold_gtsb_gtsd2%>%
  sel


##now model by gene

# Convert "gene" into an unordered factor
allfold_gtsb_gtsdtry$gene <- as.factor(allfold_gtsb_gtsdtry$gene)
allfold_gtsb_gtsdtry$nuc_position <- as.factor(allfold_gtsb_gtsdtry$nuc_position)
allfold_gtsb_gtsdtry$gene_length <- as.factor(allfold_gtsb_gtsdtry$gene_length)
class(allfold_gtsb_gtsdtry$gene_length)


sapply(lapply(allfold_gtsb_gtsd2, unique), length)

lapply(allfold_gtsb_gtsd2[c('gene', 'gene_length')], unique)

glimpse(allfold_gtsb_gtsdtry)


allfold_gtsb_gtsdtryf<-allfold_gtsb_gtsdtry%>%
  filter(gene=="gtsB")%>%
  filter(fold=="FOUR")
##finally I added the gene length and the conflict stopped and I run the model with gene and gene lenght
library(MASS)

##position
per_position_models4 <- allfold_gtsb_gtsdtryf %>%
  do(tidy(glm.nb(obser_total_mut_Mean ~nuc_position , data = .)))
per_gene_models2pos <- collect(per_position_models4)

per_gene_models4 <- allfold_gtsb_gtsdtry %>%
  do(tidy(glm.nb(obser_total_mut_Mean ~gene+gene_length+gene*gene_length , data = .)))

# Collect the results into a single data frame
per_gene_models2 <- collect(per_gene_models2)

type_diff_gene <- per_gene_models2 %>%
  mutate(padj = p.adjust(p.value, method="BH"))

###model dsicussion

fold_model_gtsb_4ffff2<-read.csv("fitness_four.csv")

# Convert factors
fold_model_gtsb_4ffff2$position <- as.factor(fold_model_gtsb_4ffff$position )
fold_model_gtsb_4ffff$mean_subs  <- as.factor(fold_model_gtsb_4ffff$mean_subs)
fold_model_gtsb_4ffff$rscu<-as.factor(fold_model_gtsb_4ffff$rscu)
fold_model_gtsb_4ffff$W_CAI<-as.factor(fold_model_gtsb_4ffff$W_CAI)

modelunof<-glm(mean_subs~position, data =fold_model_gtsb_4ffff2, family = "poisson")
summary(modelunof)


modeluno<-glm(mean_w~position*rscu, data =fold_model_gtsb_4ffff, family = "gaussian")
summary(modeluno)

emm <- emmeans(modeluno, ~ gene, data = allfold_gtsb_gtsdtry)

pairs(emm , adjust="holm")
