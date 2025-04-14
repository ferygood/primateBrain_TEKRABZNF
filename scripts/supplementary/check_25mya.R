#1. find the 9 KRAB-ZNFs and 92 TEs that are evolutionary young. These could be looked back to the figure 2C
#2. How to define 25 mya from the young group?
#3. How many of them are negative correlated? Is it different than old group?
#4. Does the 25 mya show different pattern?

library(dplyr)

#1
hsc1 <- read.csv("tables/hsC1_corr_sig.csv")
hsc2 <- read.csv("tables/hsC2_corr_sig.csv")

znf_age <- read.csv("data/hmKZNFs337_ageInfer.csv")
znf_25mya <- znf_age %>%
    filter(branch >= 10) %>%
    filter(external_gene_name %in% kznf_infer[kznf_infer$age=="young",]$external_gene_name)
znf_young <- kznf_infer %>% filter(age=="young")

#2 get Homo sapiens, Hominoidea, Hominidae, Homininae
te_age_25mya <- te_infer %>% filter(OS != "Simiiformes")
te_young <- te_infer

#3 let us use cluster 1 as example:
hsc1_young <- hsc1 %>%
    filter(geneName %in% znf_young$external_gene_name & teName %in% te_young$NM)
# 3533 positive correlation, 1324 negative correlation

hsc1_25mya <- hsc1 %>%
    filter(geneName %in% znf_25mya$external_gene_name & teName %in% te_age_25mya$NM)
# 2 positive correlation and 7 negative correlation

#4
length(unique(hsc1_young$geneName)) #87
length(unique(hsc1_young$teName)) #150

length(unique(hsc1_25mya$geneName)) #3, ZNF439, ZNF763, ZNF98
length(unique(hsc1_25mya$teName)) #9, AluYb9, AluYk11, LTR14B, LTR13, LTR13A, LTR5_Hs, LTR7C, LTR7Y, MER11C


