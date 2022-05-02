library(tidyverse)
library(readxl)
library(ggforce)
library(patchwork)
read_excel("Data/Covariates_CMV.xlsx") -> D
D %>% filter( grepl("LLD",Sample_name)) -> D

#Breadth vs age
lm(cmv_peptides_positives ~ Age, D) %>% summary()

#Breadth vs age, if we remove samples with few
lm(cmv_peptides_positives ~ Age, filter(D, cmv_peptides_positives>0) ) %>% summary()
lm(cmv_peptides_positives ~ Age, filter(D, cmv_peptides_positives>10) ) %>% summary()
#There is an increase in breadth with age (> 1 and >0), however this increase is minor with people who already have quite some
ggplot( D, aes(x=Age, y= cmv_peptides_positives) ) + geom_point() + theme_bw() + geom_smooth(method='lm') -> Plot_CMV
ggsave("Results/CMV_breadth_gain.pdf", Plot_CMV)

print("Binned analysis")

D %>% mutate(`Bins CMV peptide number` = cut(cmv_peptides_positives, breaks = c(-Inf,0,16,Inf))) -> D #16 the mean, so 0, below mean, above mean

aov(formula = D$Age ~ D$`Bins CMV peptide number`) -> M
print(summary(M))
TukeyHSD(M) -> M
M %>% print(digits=20)
D  %>% ggplot(aes(x=`Bins CMV peptide number`, y=Age)) + geom_boxplot(outlier.shape =  NA) + ggforce::geom_sina(alpha=0.4) +  theme_bw() -> Plot_1
D %>% ggplot(aes(x=cmv_peptides_positives)) + geom_histogram(bins=80) + geom_vline(xintercept = 16) + scale_y_log10() + theme_bw() + xlab("CMV peptide number") -> Plot2

Plot_1 + Plot2 -> Plot

ggsave("Results/CMV_breadth_distribution.pdf", height= 10, width=15, units="cm", Plot)

