library(twice)
library(cowplot)
library(ggpubr)
library(RCy3)
young_kznf <- kznf_infer %>% filter(age == "young")

#HmPtC1 <- readRDS("~/github/pBrain/posts/brain/CorrDEanalysis/results_rdata/HmPtC1.rds")

c1_gene_res <- HmPtC1$DEobject$gene_res %>%
    data.frame() %>%
    filter(abs(log2FoldChange) >= -1.5 & pvalue<0.05) %>%
    filter(rownames(.) %in% kznf_infer$external_gene_name)

c1_te_res <- HmPtC1$DEobject$te_res %>%
    data.frame() %>%
    filter(abs(log2FoldChange) >= 1.5 & pvalue<0.05)

hm_corr <- HmPtC1$corrRef %>%
    filter(pvalue<0.05 & coef<0) %>%
    filter(geneName %in% rownames(c1_gene_res) &
               teName %in% rownames(c1_te_res)) %>%
    mutate(pair = paste0(geneName, "-", teName))

pt_corr <- HmPtC1$corrCompare %>%
    filter(pvalue<0.05 & coef<0) %>%
    filter(geneName %in% rownames(c1_gene_res) &
               teName %in% rownames(c1_te_res)) %>%
    mutate(pair = paste0(geneName, "-", teName))

# in human
hm_corr_age <- hm_corr %>%
    mutate(kznf_age = ifelse(geneName %in% young_kznf$external_gene_name, "y", "o")) %>%
    mutate(te_age = ifelse(teName %in% te_infer$NM, "y", "o")) %>%
    mutate(age_c = paste0(kznf_age, "-", te_age))

node.hm <- data.frame(
    id = c(c(unique(hm_corr$geneName), unique(hm_corr$teName)))
)

edge.hm <- hm_corr_age[,c(1,2,9)]
colnames(edge.hm) <- c("source", "target", "interaction")

createNetworkFromDataFrames(node.hm, edge.hm)

# in chimpanzee
pt_corr_age <- pt_corr %>%
    mutate(kznf_age = ifelse(geneName %in% young_kznf$external_gene_name, "y", "o")) %>%
    mutate(te_age = ifelse(teName %in% te_infer$NM, "y", "o")) %>%
    mutate(age_c = paste0(kznf_age, "-", te_age))

node.pt <- data.frame(
    id = c(c(unique(pt_corr$geneName), unique(pt_corr$teName)))
)

edge.pt <- pt_corr_age[,c(1,2,9)]
colnames(edge.pt) <- c("source", "target", "interaction")

createNetworkFromDataFrames(node.pt, edge.pt)

