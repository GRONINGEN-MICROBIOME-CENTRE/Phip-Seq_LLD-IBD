library(tidyverse)
library(readxl)

#1. Pull metadata
read_tsv("~/Resilio Sync/Antibodies_WIS (1)/Plate records/Covariates.tsv") %>% filter(grepl("LL", ID)) %>% 
  mutate( TimePoint= ifelse( grepl("_F", ID), "Followup", "Baseline" ) ) -> Covariates
#2. Get Matrix of interest
readRDS("~/Desktop/Immuno_matrix_postSelection.rds") -> Peptides
#3. Get Oeotudes if interest
read_excel("~/Downloads/Peptide_selection_May2022_SZ.xls") %>% filter(For_LLD == "Yes") -> Peptides_selection

#4. Quantify prevalence 
Check_prevalence = function( Peptides, Samples, Selection  ){
    Peptides %>% filter(ID %in% Samples) %>% select(c(ID, Selection)) -> Plate_selection
    #Plate_selection %>% select(-ID) %>% apply(2, function(x){ sum(x)/length(x) } ) %>% print()
    Plate_selection %>% select(-ID) %>% apply(2, function(x){ sum(x) } ) %>% print()
}


Samples = Covariates %>% filter(TimePoint == "Baseline") %>% filter(plate_id == "SF29990795")   %>% select(sample_id) %>% as_vector() %>% as.vector()
Samples = paste0("32_", Samples)
Selection = Peptides_selection$Peptide
Check_prevalence(Peptides, Samples, Selection)

#5. Quantify Prevalence in a random selection enriched for the less prevalent peptide
#Do a selection
set.seed(20)
Peptides %>% filter(ID %in% Samples)  %>% filter(twist_15595 == 1) %>% select(ID) -> For_experiment_enrich
Choices = tibble()
for ( i in c(16, 21, 31, 41)){
  sample(Samples, i) -> For_experiment
  For_experiment = c(For_experiment, For_experiment_enrich$ID)
  length(For_experiment)
  Check_prevalence(Peptides, For_experiment, Selection)
  rbind(Choices, tibble(ID = For_experiment, N=length(For_experiment) )) -> Choices
}
Choices %>% filter(N==40) -> Choice
Covariates %>% mutate( ID_LLD = ID, ID = paste0("32_", sample_id, join="") ) %>% select(ID, ID_LLD) %>%
  left_join(Choice, .) %>% View()





################################################################################################################################################
##### Check Null Distributions of correlation / Similarity to assess how many matches we need to find statistical significance##################
################################################################################################################################################
set.seed(89)

jaccard <- function(a, b) {
  (a+b)/2 -> Common
  Common = Common[Common == 1]
  intersection = length( Common )
  union = sum(a) + sum(b) - intersection
  return (intersection/union)
}

#########Parameters simulation
N_samples = 30 #40
N_1s = 4 #8 
Positions_with_1 =  c(5, 10, 15, 20) #c(5, 10, 15, 20, 25, 30, 35, 40)
N_random = 10000
Frq = N_1s/N_samples
Weights = c(1-Frq, Frq )
########Generate a simulated sample, and all possible matches
#Real Sample
Real_sample =  sample( c(0), N_samples, replace=T )
Real_sample[ Positions_with_1] = 1
#Possible matches
Pseudo_real = sample( c(0), N_samples, replace=T )
#Get a correlation/similarity per possible match
Similarities = tibble()
N = 1
for (i in Positions_with_1 ){
  Pseudo_real[i] = 1
  cor(Real_sample, Pseudo_real) -> C
  jaccard(Real_sample, Pseudo_real) -> D
  Similarities = rbind(Similarities, tibble(Correlation = C, Similarity=D, N=N ))
  N = N + 1
}
##########################################################################
#Generate Null distributions for correlation coeffients/similarities######
##########################################################################
Rho_dist = c()
Distances = c()
for (i in seq(N_random)){
  Random_sample = sample( c(0,1), N_samples, replace=T, prob=Weights )
  Rho = cor(Real_sample, Random_sample)
  Rho_dist = c(Rho_dist,  Rho)
  D = jaccard(Real_sample, Random_sample)
  Distances = c(Distances,  D)
}
Rho_dist[! is.na(Rho_dist)] -> Rho_dist
ggplot() + geom_density(aes(x=Rho_dist)) + theme_bw()
ggplot() + geom_density(aes(x=Distances)) + theme_bw()

#Compute the P of the possible matches, accounting for the null distribution
P_corr = c()
P_dist = c()
for (M in Similarities$N){
  Matches = Similarities %>% filter(N == M) %>% select(Correlation) %>% as_vector() %>% as.vector()
  Matches2 = Similarities %>% filter(N == M) %>% select(Similarity) %>% as_vector() %>% as.vector()
  P_corr = c(P_corr, sum(Rho_dist > Matches) / length(Rho_dist) )
  P_dist = c(P_dist, sum(Distances > Matches2) / length(Distances) )
}
Similarities %>% mutate(P_distance = P_dist) -> Similarities
Similarities %>% mutate(P_correlation = P_corr) -> Similarities
print(Similarities)


### Results 8 / 40

#With sampling proababilities 0.5/0.5 -> 100000 samples
#Correlation Similarity     N P_correlation P_distance
#<dbl>      <dbl> <dbl>         <dbl>      <dbl>
# 1       0.320      0.125     1       0.0233     0.730  
#2       0.459      0.25      2       0.00101    0.0904 
#3       0.569      0.375     3       0.00002    0.00182
#4       0.667      0.5       4       0          0      
#5       0.756      0.625     5       0          0      
#6       0.840      0.75      6       0          0      
#7       0.921      0.875     7       0          0      
#8       1          1         8       0          0      


#If Prob 1 8/40 --> 10000 samples
#Correlation Similarity     N P_correlation P_distance
#<dbl>      <dbl> <dbl>         <dbl>      <dbl>
#  1       0.320      0.125     1      0.0298      0.409  
#2       0.459      0.25      2      0.00460     0.0549 
#3       0.569      0.375     3      0.000700    0.00630
#4       0.667      0.5       4      0           0      
#5       0.756      0.625     5      0           0      
#6       0.840      0.75      6      0           0      
#7       0.921      0.875     7      0           0      
#8       1          1         8      0           0      




