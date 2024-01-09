library(ggplot2)
library(ggvenn)
library(ggpubr)
library(RCy3)
library(twice)

cor <- readRDS("~/github/pBrain/data/hmCorr_tpm.rds")
rmsk <- read.csv("~/github/pBrain/data/hg38rmsk_info.csv")

kznf_infer_young <- kznf_infer %>%
    filter(age=="young")

data("hg38rmsk_info")
hg38TE <- hg19rmsk_info

# select human specific TE:KRAB-ZNF from primate brain data
hmc1 <- cor$c1 %>%
    filter(pvalue<0.05 & coef<0) %>%
    mutate(pair=paste0(geneName, "-", teName))

ptc1 <- HmPtC1$corrCompare %>%
    filter(pvalue<0.05 & coef<0) %>%
    mutate(pair=paste0(geneName, "-", teName))

ppc1 <- HmPpC1$corrCompare %>%
    filter(pvalue<0.05 & coef<0) %>%
    mutate(pair=paste0(geneName, "-", teName))

mmc1 <- HmMmC1$HmMmC1_MmCorr %>%
    filter(pvalue<0.05 & coef<0) %>%
    mutate(pair=paste0(geneName, "-", teName))

nhp_union_c1 <- union(union(ptc1$pair, ppc1$pair), mmc1$pair)

# select human specific link
hs_link_c1 <- setdiff(hmc1$pair, nhp_union_c1)

# 2. control sample
tcx_control <- mayoTEKRABber$tcxControlCorr %>%
    filter(pvalue<0.05 & coef<0) %>%
    mutate(pair=paste0(geneName, "-", teName))

tcx_control_hslink <- tcx_control %>%
    filter(pair %in% hs_link_c1)

tcx_AD <- mayoTEKRABber$tcxADCorr %>%
    filter(pvalue<0.05 & coef<0) %>%
    mutate(pair=paste0(geneName, "-", teName))

lost_c1 <- setdiff(tcx_control_hslink$pair, tcx_AD$pair)

df_result <- tcx_control %>%
    filter(pair %in% lost_c1) %>%
    mutate(kznf_age = ifelse(geneName %in% kznf_infer_young$external_gene_name,
                             "young", "old")) %>%
    mutate(te_age = ifelse(teName %in% te_infer$NM, "young", "old")) %>%
    mutate(count = paste0(kznf_age, "-", te_age)) %>%
    left_join(hg38TE[,c(1,2)], join_by(teName==gene_id))

write.csv(df_result, file="tables/tcx_lostAD_link.csv", row.names = F)

# see 2492 TE:KRAB-ZNF in AD sample
lost_list <- df_result$pair

tcxCorr_ad <- mayoTEKRABber$tcxADCorr %>%
    mutate(pair = paste0(geneName, "-", teName)) %>%
    filter(pair %in% lost_list)

# create network
node.tcx <- data.frame(
    id = c(c(unique(hs_tcx$geneName), unique(hs_tcx$teName)))
)

yy_kznfs_TE <- c(hs_tcx_yy$geneName, hs_tcx_yy$teName)
node.tcx <- node.tcx %>%
    mutate(age = ifelse(id %in% yy_kznfs_TE, "young", "old"))


edge.tcx <- hs_tcx
edge.tcx <- edge.tcx %>%
    mutate(up = ifelse(pair %in% hs_tcx_select$pair, 1, 0)) %>%
    select(c(1,2,9,11))
colnames(edge.tcx) <- c("source", "target", "interaction", "DE")

createNetworkFromDataFrames(node.tcx, edge.tcx)
