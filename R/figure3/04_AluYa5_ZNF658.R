library(twice)
library(cowplot)
library(ggpubr)
kznf_age_infer <- read.csv("/path/to/hmKZNFs337_ageInfer.csv")

exp <- readRDS("/path/to/primateBrain_kznfTEexp.rds")
HmC1_id <- colnames(HmPtC1$HmPtC1_DE$geneCorrInputRef)

df_exp <- exp[HmC1_id, c("ZNF658", "AluYa5")]
g <- ggplot(df_exp, aes(x=log(ZNF658), y=log(AluYa5))) +
    geom_point(color="black", fill="#E4D4B3", shape=21, size=3) +
    geom_smooth(method="lm", se = FALSE) +
    xlab("ZNF658") +
    ylab("AluYa5") +
    theme_bw()
