library(ggvenn)
library(ggpubr)
library(ggplot2)

load("/path/to/data/mayo_cntTable/mayoTEKRABber_balance.RData")

kznf_gene_infer <- read.csv("/path/to/data/hmKZNFs337_ageInfer.csv")
young_kznf <- kznf_gene_infer %>%
    filter(branch>=8)
te_infer <- read.csv("/path/to/data/Dfam_TE_simiiformes.csv")

cbeControlCor <- mayoTEKRABber$cbeControlCorr %>%
    filter(pvalue<0.05 & coef<0) %>%
    mutate(pair=paste0(geneName, "-", teName))

cbeADCor <- mayoTEKRABber$cbeADCorr %>%
    filter(pvalue<0.05 & coef<0) %>%
    mutate(pair=paste0(geneName, "-", teName))

tcxControlCor <- mayoTEKRABber$tcxControlCorr %>%
    filter(pvalue<0.05 & coef<0) %>%
    mutate(pair=paste0(geneName, "-", teName))

tcxADCor <- mayoTEKRABber$tcxADCorr %>%
    filter(pvalue<0.05 & coef<0) %>%
    mutate(pair=paste0(geneName, "-", teName))

pair_list <- list(
    cbe_control = cbeControlCor$pair,
    cbe_AD = cbeADCor$pair,
    tcx_control = tcxControlCor$pair,
    tcx_AD = tcxADCor$pair
)

vc <- ggvenn(list(CBE_control=cbeControlCor$pair, CBE_AD=cbeADCor$pair), fill_color = c("#6699FF", "#FFCC66"))

vt <- ggvenn(list(TCX_control=tcxControlCor$pair, TCX_AD=tcxADCor$pair), fill_color = c("#CC99FF", "#99CC99"))

vall <- ggvenn(pair_list, fill_color = c("#6699FF", "#FFCC66", "#CC99FF", "#99CC99"))

gall <- ggarrange(vc, vt, vall, ncol=3)
