library(ggplot2)
library(hrbrthemes)
library(twice)
library(stringr)
library(ggpubr)

data("hg38rmsk_info")

# read table
tcx_lost <- read.csv("tables/tcx_lostAD_link.csv")

# let see the distribution of these TEs in both sample
te_list <- unique(tcx_lost$teName) #551

tcx_de <- mayoTEKRABber$tcxDE$te_res %>% data.frame()

tcx_de_filter <- tcx_de %>%
    filter(rownames(.) %in% te_list)

tcx_de_filter <- tcx_de_filter %>%
    mutate(group = ifelse(pvalue<0.05, "sig", "no-sig")) %>%
    mutate(exp = ifelse(log2FoldChange < 0, "down-regulated in AD", "up-regulated in AD"))

df_tcx_selectTE <- tcx_de_filter %>%
    select(c(6,7)) %>%
    mutate(teName = rownames(.)) %>%
    select(c(3,1,2))
rownames(df_tcx_selectTE) <- 1:nrow(df_tcx_selectTE)

df_tcx_selectTE <- df_tcx_selectTE %>%
    left_join(hg19rmsk_info, join_by(teName==gene_id))

# we then pick Alu, hAT elements, ERVs, L1, SVAs
df_tcx_selectTE_subset <- df_tcx_selectTE %>%
    mutate(TEgroup = case_when(
        family_id == "Alu" ~ "Alu elements",
        str_detect(family_id, "^hAT") ~ "hAT elements",
        family_id == "L1" ~ "L1",
        str_detect(family_id, "^ERV") ~ "ERVs",
        str_detect(family_id, "^SVA") ~ "SVAs"
    ))

df_tcx_selectTE_subset <- na.omit(df_tcx_selectTE_subset)

contigency_table <- table(df_tcx_selectTE_subset$exp, df_tcx_selectTE_subset$TEgroup)

fisher_test_result <- fisher.test(contigency_table)

df_tcx_selectTE_subset <- df_tcx_selectTE_subset %>%
    mutate(exp_numeric = ifelse(exp=="up-regulated in AD", 1, -1))

df_summary <- df_tcx_selectTE_subset %>%
    select(c(6, 7)) %>%
    group_by(TEgroup, exp_numeric) %>%
    summarise(total_exp = sum(exp_numeric))

df_summary$exp_numeric <- ifelse(df_summary$exp_numeric==1, "up-regulated in AD", "down-regulated in AD")
df_summary$total_exp <- abs(df_summary$total_exp)

df_summary$exp_numeric <- factor(df_summary$exp_numeric,
                                 levels=c("up-regulated in AD",
                                          "down-regulated in AD"))

p <- ggplot(df_summary, aes(x=TEgroup, y=total_exp, fill=exp_numeric)) +
    geom_bar(stat = "identity", position=position_dodge(width=0.9), color="black") +
    scale_fill_manual(values=c("#d74c4f", "#0a75ad")) +
    labs(x="", y="number of TEs", fill="") +
    theme_bw()

