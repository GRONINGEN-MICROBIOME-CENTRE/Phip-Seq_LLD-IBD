library(tidyverse)

process =  function(Data){
        #cleanup of the data before other functions can be applied
        Data %>% select(-ID) -> Data
        Data[ , apply(Data, 2, function(x) !any(is.na(x)))] -> Data
        return(Data)
}


Data = readRDS("Data/Immuno_matrix_postSelection.rds")
EBV = read_tsv("Data/EBV_peptides.txt", col_names = FALSE)
Data %>% select(c("ID", EBV$X1) ) -> Data_EBV
Data_EBV %>% select(-ID) %>% apply(1, function(x){ y = x[!is.na(x)] ; sum(y) } ) -> EBV_breadth
Data_EBV %>% mutate(EBV_breadth = EBV_breadth) %>% select(ID, EBV_breadth) -> Data_EBV
process(Data) -> Data_all

#Find local minima of the kernel density
density(Data_EBV$EBV_breadth) -> D
tibble(X = D$x, y= D$y) -> De
De %>% filter(X > 0 & X<60) %>% arrange(y)
Cluster = Data_EBV$EBV_breadth > 10.6

mutate(Data_EBV,C =as.factor(Cluster))  %>% ggplot(aes(x=EBV_breadth, fill =C)) + geom_histogram() + theme_bw() +geom_vline(xintercept =  11)

write_tsv(mutate(Data_EBV, EBV_cluster = Cluster),  "Data/EBV_clusters.txt", col_names = T)


