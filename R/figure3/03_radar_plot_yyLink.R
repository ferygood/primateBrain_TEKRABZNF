library(TEKRABber)
library(twice)

data("hmKZNFs337")

meta <- metadata %>%
    left_join(brain_meta, join_by(brain_region==region))

kznf_age_infer <- read.csv("/path/to/hmKZNFs337_ageInfer.csv")
kznf_age_infer$age <- ifelse(kznf_age_infer$branch >= 8, "young", "old")

te_age_infer <- read.csv("/path/to/Dfam_TE_simiiformes.csv")

read_rds_files("/path/to/results_rdata/")

convert_tpm <- function(sample_id){

    # convert to tpm: gene and TE
    # gene
    gene_counts <- colSums(hmGene[, -1])
    scaling_gene <- gene_counts / 1e6
    hmGene_tpm <- hmGene
    hmGene_tpm[, -1] <- t(t(hmGene[,-1])/ scaling_gene + 1) * 1e6
    hmGene_tpm <- hmGene_tpm %>%
        inner_join(hmKZNFs337, join_by(geneID == ensembl_gene_id))
    rownames(hmGene_tpm) <- hmGene_tpm$external_gene_name

    hmGene_tpm <- hmGene_tpm[,c(sample_id)]

    # TE
    te_counts <- colSums(hmTE[,-c(1,2,3)])
    scaling_te <- te_counts/ 1e6
    hmTE_tpm <- hmTE[,-c(1,2,3)]
    rownames(hmTE_tpm) <- hmTE$name
    hmTE_tpm[,-c(1,2,3)] <- t(t(hmTE[,-c(1,2,3)] ) / scaling_te + 1) * 1e6

    hmTE_tpm <- hmTE_tpm[,c(sample_id)]

    tpm_result <- list("df_gene"=hmGene_tpm, "df_te"=hmTE_tpm)

    tpm_result
}


test_result <- convert_tpm(sampleID$Run)

sample_counts <- colSums(hmGene[, -1])

scaling_factor <- sample_counts / 1e6

hmGene_tpm <- hmGene
hmGene_tpm[, -1] <- t(t(hmGene[, -1] )/ scaling_factor + 1) * 1e6

hmGene_tpm <- hmGene_tpm %>%
    inner_join(hmKZNFs337, join_by(geneID==ensembl_gene_id))

rownames(hmGene_tpm) <- hmGene_tpm$external_gene_name
hmGene_tpm <- hmGene_tpm[,c(2:133)]

# convert te count to TPM
te_counts <- colSums(hmTE[, -c(1,2,3)])
scaling_te <- te_counts / 1e6
hmTE_tpm <- hmTE
hmTE_tpm[,-c(1,2,3)] <- t(t(hmTE[,-c(1,2,3)] )/ scaling_te + 1) * 1e6
rownames(hmTE_tpm) <- hmTE_tpm$name
hmTE_tpm <- hmTE_tpm[,-c(1,2,3)]

yy_kznfsTEs <- function(cluster_name, nhp_list){
    # 1. use brain region2 for testing
    sample_id <- meta %>%
        filter(cluster==cluster_name & Organism=="Homo sapiens") %>%
        select(Run)
    print(length(sample_id$Run))

    # 2. select expression from gene and TE
    hmGene_tpm_select <- hmGene_tpm %>% select(c(sample_id$Run))
    hmTE_tpm_select <- hmTE_tpm %>% select(c(sample_id$Run))

    # 3. Use TEKRABber to calculate correlation
    hmCorr <- corrOrthologTE(hmGene_tpm_select, hmTE_tpm_select)

    # 4. filter result and add age inference
    hmCorr_sig <- hmCorr %>%
        filter(pvalue<0.05 & coef<0)

    hmCorr_sig_annot <- hmCorr_sig %>%
        mutate(pair = paste0(geneName, '-', teName)) %>%
        left_join(kznf_age_infer[,c(2,4)], join_by(geneName==external_gene_name)) %>%
        mutate(te_age = ifelse(teName %in% te_age_infer$NM, "young", "old")) %>%
        mutate(group = case_when(
            age=="old" & te_age=="old" ~ "o-o",
            age=="young" & te_age=="old" ~ "y-o",
            age=="old" & te_age=="young" ~ "o-y",
            age=="young" & te_age=="young" ~ "y-y"
        ))

    # 5. calculate stats
    kznf_num <- length(unique(hmCorr_sig_annot$geneName))
    te_num <- length(unique(hmCorr_sig_annot$teName))
    young_kznfs <- hmCorr_sig_annot %>%
        filter(age=="young") %>%
        pull(geneName) %>%
        unique() %>%
        length()
    young_TEs <- hmCorr_sig_annot %>%
        filter(te_age=="young") %>%
        pull(teName) %>%
        unique() %>%
        length()
    yy_count <- hmCorr_sig_annot %>%
        filter(group=="y-y") %>%
        nrow()
    yo_count <- hmCorr_sig_annot %>%
        filter(group=="y-o") %>%
        nrow()
    oy_count <- hmCorr_sig_annot %>%
        filter(group=="o-y") %>%
        nrow()
    oo_count <- hmCorr_sig_annot %>%
        filter(group=="o-o") %>%
        nrow()

    total_count <- yy_count + yo_count + oy_count + oo_count

    #6. See if there are human specific yy counts
    #nhp_list <- list("pt" = HmPtC2$HmPtC2_PtCorr,
    #                 "pp" = HmPpC2$HmPpC2_PpCorr,
    #                 "mm" = HmMmC2$HmMmC2_MmCorr)

    hm_yy_pair <- hmCorr_sig_annot %>% filter(group=="y-y") %>% pull(pair)

    pt_pair <- nhp_list$pt %>%
        filter(pvalue<0.05 & coef<0) %>%
        mutate(pair=paste0(geneName, "-", teName)) %>%
        pull(pair)

    pp_pair <- nhp_list$pp %>%
        filter(pvalue<0.05 & coef<0) %>%
        mutate(pair=paste0(geneName, "-", teName)) %>%
        pull(pair)

    mm_pair <- nhp_list$mm %>%
        filter(pvalue<0.05 & coef<0) %>%
        mutate(pair=paste0(geneName, "-", teName)) %>%
        pull(pair)

    u_list <- union(union(pt_pair, pp_pair), mm_pair)

    hm_specific_yypair <- setdiff(hm_yy_pair, u_list)
    hm_specific_yypair_count <- length(hm_specific_yypair)
    hm_specific_yypair_percentage <- round(hm_specific_yypair_count * 100 / yy_count, 2)

    output_list <- c(kznf_num, te_num, young_kznfs, young_TEs,
                     total_count, yy_count, yo_count, oy_count, oo_count,
                     hm_specific_yypair_count, hm_specific_yypair_percentage)

    output_list
}

# cluster 1
nhp_list_c1 <- list("pt"=HmPtC1$HmPtC1_PtCorr,
                    "pp"=HmPpC1$HmPpC1_PpCorr,
                    "mm"=HmMmC1$HmMmC1_MmCorr)
c1 <- yy_kznfsTEs("cluster1", nhp_list_c1)

# cluster 2
nhp_list_c2 <- list("pt"=HmPtC2$HmPtC2_PtCorr,
                    "pp"=HmPpC2$HmPpC2_PpCorr,
                    "mm"=HmMmC2$HmMmC2_MmCorr)
c2 <- yy_kznfsTEs("cluster2", nhp_list_c2)

# cluster 3
nhp_list_c3 <- list("pt"=HmPtC3$HmPtC3_PtCorr,
                    "pp"=HmPpC3$HmPpC3_PpCorr,
                    "mm"=HmMmC3$HmMmC3_MmCorr)
c3 <- yy_kznfsTEs("cluster3", nhp_list_c3)

# cluster 4
nhp_list_c4 <- list("pt"=HmPtC4$HmPtC4_PtCorr,
                    "pp"=HmPpC4$HmPpC4_PpCorr,
                    "mm"=HmMmC4$HmMmC4_MmCorr)
c4 <- yy_kznfsTEs("cluster4", nhp_list_c4)

# cluster 5
nhp_list_c5 <- list("pt"=HmPtC5$HmPtC5_PtCorr,
                    "pp"=HmPpC5$HmPpC5_PpCorr,
                    "mm"=HmMmC5$HmMmC5_MmCorr)
c5 <- yy_kznfsTEs("cluster5", nhp_list_c5)

# cluster 6
nhp_list_c6 <- list("pt"=HmPtC6$HmPtC6_PtCorr,
                    "pp"=HmPpC6$HmPpC6_PpCorr,
                    "mm"=HmMmC6$HmMmC6_MmCorr)
c6 <- yy_kznfsTEs("cluster6", nhp_list_c6)

# cluster 1
nhp_list_c7 <- list("pt"=HmPtC7$HmPtC7_PtCorr,
                    "pp"=HmPpC7$HmPpC7_PpCorr,
                    "mm"=HmMmC7$HmMmC7_MmCorr)
c7 <- yy_kznfsTEs("cluster7", nhp_list_c7)

rows <- matrix(c(c1, c2, c3, c4, c5, c6, c7), ncol=11, byrow = TRUE)
columnNames <- c("total_kznfs", "total_TEs", "young_kznfs","young_TEs",
                 "total_correlation", "yy_count",
                 "yo_count", "oy_count", "oo_count",
                 "hm_specific_yy", "hm_specific_yy_percentage")
df <- data.frame(rows)
colnames(df) <- columnNames

brain_value <- c("primary & secondary\n cortices", "limbic & association\n cortices", "archicortex", "hypothalamus & thalamus", "cerebellar \n white matter", "cerebellar \n grey matter", "striatum")

df$brain_region <- brain_value
df <- df[,c(12, 1:11)]

write.csv(df, file="/path/to/yy_kznfTEs_table_hm.csv",row.names = FALSE)

library(ggplot2)
library(geomtextpath)
library(dplyr)
library(tidyr)

yydf <- read.csv("/path/to/yy_kznfTEs_table_hm.csv")

yydf_select <- yydf[,c(1,7,11)]
yydf_select <- yydf_select %>%
    pivot_longer(!brain_region, names_to = "group", values_to = "count")

g <- ggplot(yydf_select, aes(brain_region, count)) +
    geom_col(aes(fill=group), position = position_dodge(width=1)) +
    geom_vline(xintercept = 1:8 - 0.5, color = "darkgray") +
    geom_hline(yintercept = 0:5 * 5, color = "darkgray") +
    scale_fill_manual(values = c("#E4D4B3", "#E59E00")) +
    ggtitle("Percentage of young Human specific TE:KRAB-ZNF in all Young TE:KRAB-ZNF") +
    theme_bw() +
    coord_curvedpolar() +
    theme(axis.text.x = element_text(color = "black", size = 14),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size=24))
