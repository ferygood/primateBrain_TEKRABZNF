# Figure 2B Volcanoplot

library(twice)
library(ggpubr)
library(ggplot2)
library(tidyverse)

data("hmKZNFs337")
data("hg19rmsk_info")

# import brain data
read_rds_files("~/github/pBrain/posts/brain/CorrDEanalysis/results_rdata/")

HmPtC1_resKZNFs <- HmPtC1$HmPtC1_DE$gene_res %>%
    data.frame() %>%
    filter(rownames(.) %in% hmKZNFs337$external_gene_name) %>%
    mutate(log2FoldChange = log2FoldChange * -1)
HmPtC1_resTE <- HmPtC1$HmPtC1_DE$te_res %>%
    data.frame() %>%
    mutate(log2FoldChange = log2FoldChange * -1)


volcano_plot <- function(df){

    v <- EnhancedVolcano::EnhancedVolcano(
        df,
        lab = rownames(df),
        x = 'log2FoldChange',
        y = 'pvalue',
        labSize = 3.0,
        FCcutoff = 1.5,
        pCutoff = 0.05,
        legendPosition = 'none',
        legendLabSize = 10,
        title = "",
        subtitle = "",
        caption = ""
    )

    v

}

geneVol <- volcano_plot(HmPtC1_resKZNFs)
teVol <- volcano_plot(HmPtC1_resTE)
gvol <- ggarrange(geneVol, teVol, ncol=1)
ggsave(filename="figures/volcano_c1.jpg", dpi=500, width=4.5, height=10)

