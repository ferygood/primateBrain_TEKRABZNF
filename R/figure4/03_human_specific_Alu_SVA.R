library(RCy3)

hmCorr_tpm <- readRDS("/path/to/hmCorr_tpm.rds")

hm_c1 <- hmCorr_tpm$c1 %>%
    filter(coef < 0 & pvalue < 0.05) %>%
    mutate(pair = paste0(geneName, "-", teName))

# read NHPs
pt_c1 <- HmPtC1$corrCompare %>%
    filter(coef < 0 & pvalue < 0.05) %>%
    mutate(pair = paste0(geneName, "-", teName))

pp_c1 <- HmPpC1$corrCompare %>%
    filter(coef < 0 & pvalue < 0.05) %>%
    mutate(pair = paste0(geneName, "-", teName))

mm_c1 <- HmMmC1$HmMmC1_MmCorr %>%
    filter(coef < 0 & pvalue < 0.05) %>%
    mutate(pair = paste0(geneName, "-", teName))

uni_pair <- union(union(pt_c1$pair, pp_c1$pair), mm_c1$pair)

# filter only human
hm_c1_hmonly <- hm_c1 %>%
    filter(!pair %in% uni_pair)

hg38TE <- read.csv("/path/to/hg38rmsk_info.csv")
youngTE <- read.csv("/path/to/Dfam_TE_simiiformes.csv")
kznf_age <- read.csv("/path/to/hmKZNFs337_ageInfer.csv")
young_kznf <- kznf_infer %>% filter(age == "young")

alu_id <- hg38TE %>%
    filter(repFamily == "Alu")

hm_c1_hmonly_yy <- hm_c1_hmonly %>%
    filter(teName %in% youngTE$NM) %>%
    filter(teName %in% alu_id$repName) %>%
    filter(geneName %in% young_kznf$external_gene_name) %>%
    left_join(kznf_infer[,c(2,6)], join_by(geneName==external_gene_name)) %>%
    left_join(youngTE[,c(3,4)], join_by(teName==NM))

n1.node <- data.frame(
    id=c(c(unique(hm_c1_hmonly_yy$geneName), unique(hm_c1_hmonly_yy$teName))))
n1.edge <- hm_c1_hmonly_yy[,c(1,2,4)]
colnames(n1.edge) <- c("source", "target", 'interaction')

createNetworkFromDataFrames(n1.node, n1.edge)

sva_id <- c("SVA_A", "SVA_B", "SVA_C", "SVA_D", "SVA_E", "SVA_F")

hm_c1_hmonly_yy_SVA <- hm_c1_hmonly %>%
    filter(teName %in% youngTE$NM) %>%
    filter(teName %in% sva_id) %>%
    filter(geneName %in% young_kznf$external_gene_name) %>%
    left_join(kznf_infer[,c(2,6)], join_by(geneName==external_gene_name)) %>%
    left_join(youngTE[,c(3,4)], join_by(teName==NM))

n2.node <- data.frame(
    id=c(c(unique(hm_c1_hmonly_yy_SVA$geneName),
           unique(hm_c1_hmonly_yy_SVA$teName))))
n2.edge <- hm_c1_hmonly_yy_SVA[,c(1,2,4)]
colnames(n2.edge) <- c("source", "target", 'interaction')

createNetworkFromDataFrames(n2.node, n2.edge)
