library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)

hmExp_tpm <- readRDS("/path/to/hmExp_tpm.rds")
hmCorr_tpm <- readRDS("/path/to/hmCorr_tpm.rds")
read_rds_files("/path/to/results_rdata/")
te_infer <- read.csv("/path/to/Dfam_TE_simiiformes.csv")

detect_TEexp_level <- function(hm_te_res, hmCorr_tpm_c, hmExp_tpm_c, xlab_species_cluster){

    # get data
    de_te <- hm_te_res %>%
        data.frame() %>%
        mutate(log2FoldChange = log2FoldChange * (-1)) %>%
        mutate(teName = rownames(.))

    # correlation results
    corr <- hmCorr_tpm_c %>%
        filter(coef<0 & pvalue<0.05)
    corr_TEs <- unique(corr$teName)

    # expression data
    te_exp <- hmExp_tpm_c
    te_exp_mean <- te_exp %>%
        mutate(mean_exp = rowMeans(.)) %>%
        select(ncol(.)) %>%
        mutate(teName = rownames(.)) %>%
        left_join(de_te[,c(2,4,6)], join_by(teName==teName))

    te_exp_mean_nona <- na.omit(te_exp_mean)
    rownames(te_exp_mean_nona) <- 1:nrow(te_exp_mean_nona)

    te_exp_mean_nona <- te_exp_mean_nona %>%
        mutate(group = ifelse(
            (abs(log2FoldChange) >= 1.5) &
                (pvalue<0.05) &
                (teName %in% corr_TEs), "DE-TE:KRAB-ZNF", "others"))

    gg <- ggplot(te_exp_mean_nona, aes(x=log2FoldChange, y=mean_exp)) +
        geom_point(size=3, shape=21, colour="black", aes(fill=group)) +
        scale_fill_manual(values = c("#D55E00", "#999999")) +
        xlab(paste0(xlab_species_cluster)) +
        ylab("") +
        geom_vline(xintercept = c(1.5, -1.5), linetype="dashed", color="#56B4E9") +
        theme_bw() +
        theme(legend.position = "bottom")

    nrow_count <- te_exp_mean_nona %>% filter(group=="DE-TE:KRAB-ZNF") %>% nrow()
    nrow_base <- te_exp_mean_nona %>% filter(abs(log2FoldChange) >= 1.5 & pvalue < 0.05) %>% nrow()
    percentage <- round( 100*nrow_count / nrow_base, 2)

    result <- list(
        plot=gg,
        percentage=percentage
    )

    result

}

# comparing human to chimpanzee
hmptc1_TElowexp <- detect_TEexp_level(HmPtC1$HmPtC1_DE$te_res, hmCorr_tpm$c1, hmExp_tpm$c1_tpm$te, "primary & secondary cortices")
hmptc2_TElowexp <- detect_TEexp_level(HmPtC2$HmPtC2_DE$te_res, hmCorr_tpm$c2, hmExp_tpm$c2_tpm$te, "limbic & association cortices")
hmptc3_TElowexp <- detect_TEexp_level(HmPtC3$HmPtC3_DE$te_res, hmCorr_tpm$c3, hmExp_tpm$c3_tpm$te, "archicortex")
hmptc4_TElowexp <- detect_TEexp_level(HmPtC4$HmPtC4_DE$te_res, hmCorr_tpm$c4, hmExp_tpm$c4_tpm$te, "thalamus & hypothalamus")
hmptc5_TElowexp <- detect_TEexp_level(HmPtC5$HmPtC5_DE$te_res, hmCorr_tpm$c5, hmExp_tpm$c5_tpm$te, "cerebellar white matter")
hmptc6_TElowexp <- detect_TEexp_level(HmPtC6$HmPtC6_DE$te_res, hmCorr_tpm$c6, hmExp_tpm$c6_tpm$te, "cerebellar grey matter")
hmptc7_TElowexp <- detect_TEexp_level(HmPtC7$HmPtC7_DE$te_res, hmCorr_tpm$c7, hmExp_tpm$c7_tpm$te, "striatum")


# combine plot
g_hmpt <- ggarrange(hmptc1_TElowexp$plot, hmptc2_TElowexp$plot, hmptc3_TElowexp$plot,
                    hmptc4_TElowexp$plot, hmptc5_TElowexp$plot, hmptc6_TElowexp$plot,
                    hmptc7_TElowexp$plot, ncol=3, nrow=3, common.legend = TRUE)

# comparing human to bonobo
hmppc1_TElowexp <- detect_TEexp_level(HmPpC1$HmPpC1_DE$te_res, hmCorr_tpm$c1, hmExp_tpm$c1_tpm$te, "primary & secondary cortices")
hmppc2_TElowexp <- detect_TEexp_level(HmPpC2$HmPpC2_DE$te_res, hmCorr_tpm$c2, hmExp_tpm$c2_tpm$te, "limbic & association cortices")
hmppc3_TElowexp <- detect_TEexp_level(HmPpC3$HmPpC3_DE$te_res, hmCorr_tpm$c3, hmExp_tpm$c3_tpm$te, "archicortex")
hmppc4_TElowexp <- detect_TEexp_level(HmPpC4$HmPpC4_DE$te_res, hmCorr_tpm$c4, hmExp_tpm$c4_tpm$te, "thalamus & hypothalamus")
hmppc5_TElowexp <- detect_TEexp_level(HmPpC5$HmPpC5_DE$te_res, hmCorr_tpm$c5, hmExp_tpm$c5_tpm$te, "cerebellar white matter")
hmppc6_TElowexp <- detect_TEexp_level(HmPpC6$HmPpC6_DE$te_res, hmCorr_tpm$c6, hmExp_tpm$c6_tpm$te, "cerebellar grey matter")
hmppc7_TElowexp <- detect_TEexp_level(HmPpC7$HmPpC7_DE$te_res, hmCorr_tpm$c7, hmExp_tpm$c7_tpm$te, "striatum")


# combine plot
g_hmpp <- ggarrange(hmppc1_TElowexp$plot, hmppc2_TElowexp$plot, hmppc3_TElowexp$plot,
                    hmppc4_TElowexp$plot, hmppc5_TElowexp$plot, hmppc6_TElowexp$plot,
                    hmppc7_TElowexp$plot, ncol=3, nrow=3, common.legend = TRUE)

# comparing human to macaque
hmmmc1_TElowexp <- detect_TEexp_level(HmMmC1$HmMmC1_DE$te_res, hmCorr_tpm$c1, hmExp_tpm$c1_tpm$te, "primary & secondary cortices")
hmmmc2_TElowexp <- detect_TEexp_level(HmMmC2$HmMmC2_DE$te_res, hmCorr_tpm$c2, hmExp_tpm$c2_tpm$te, "limbic & association cortices")
hmmmc3_TElowexp <- detect_TEexp_level(HmMmC3$HmMmC3_DE$te_res, hmCorr_tpm$c3, hmExp_tpm$c3_tpm$te, "archicortex")
hmmmc4_TElowexp <- detect_TEexp_level(HmMmC4$HmMmC4_DE$te_res, hmCorr_tpm$c4, hmExp_tpm$c4_tpm$te, "thalamus & hypothalamus")
hmmmc5_TElowexp <- detect_TEexp_level(HmMmC5$HmMmC5_DE$te_res, hmCorr_tpm$c5, hmExp_tpm$c5_tpm$te, "cerebellar white matter")
hmmmc6_TElowexp <- detect_TEexp_level(HmMmC6$HmMmC6_DE$te_res, hmCorr_tpm$c6, hmExp_tpm$c6_tpm$te, "cerebellar grey matter")
hmmmc7_TElowexp <- detect_TEexp_level(HmMmC7$HmMmC7_DE$te_res, hmCorr_tpm$c7, hmExp_tpm$c7_tpm$te, "striatum")


# combine plot
g_hmmm <- ggarrange(hmmmc1_TElowexp$plot, hmmmc2_TElowexp$plot, hmmmc3_TElowexp$plot,
                    hmmmc4_TElowexp$plot, hmmmc5_TElowexp$plot, hmmmc6_TElowexp$plot,
                    hmmmc7_TElowexp$plot, ncol=3, nrow=3, common.legend = TRUE)

# supplementary figure S5
df_heatmap <- data.frame(
    hs_pt <- c(hmptc1_TElowexp$percentage, hmptc2_TElowexp$percentage, hmptc3_TElowexp$percentage,
               hmptc4_TElowexp$percentage, hmptc5_TElowexp$percentage, hmptc6_TElowexp$percentage,
               hmptc7_TElowexp$percentage),
    hs_pp <- c(hmppc1_TElowexp$percentage, hmppc2_TElowexp$percentage, hmppc3_TElowexp$percentage,
               hmppc4_TElowexp$percentage, hmppc5_TElowexp$percentage, hmppc6_TElowexp$percentage,
               hmppc7_TElowexp$percentage),
    hs_mm <- c(hmmmc1_TElowexp$percentage, hmmmc2_TElowexp$percentage, hmmmc3_TElowexp$percentage,
               hmmmc4_TElowexp$percentage, hmmmc5_TElowexp$percentage, hmmmc6_TElowexp$percentage,
               hmmmc7_TElowexp$percentage)
)

colnames(df_heatmap) <- c("vs. Chimpanzee", "vs. Bonobo", "vs. Macaque")
rownames(df_heatmap) <- c("primary & secondary cortices", "limbic & association cortices",
                          "archicortex", "thalamus & hypothalamus", "cerebellar white matter",
                          "cerebellar grey matter","striatum")

col_fun <- colorRamp2(c(50, 75, 100), c("blue", "white", "brown"))

h_percentage <- Heatmap(
    df_heatmap,
    cluster_rows = F,
    cluster_columns = F,
    column_names_rot = 70,
    width = 3 * unit(12, "mm"),
    height = 7 * unit(12, "mm"),
    rect_gp = gpar(col="black", lwd=0.8),
    col = col_fun,
    cell_fun = function(j, i, x, y, width, height, fill){
        grid.text(sprintf("%.2f", df_heatmap[i, j]), x, y, gp=gpar(fontsize=10))
    },
    heatmap_legend_param = list(
        title="DE-TE:KRAB-ZNF percentage",
        title_position = "lefttop-rot"
    )
)
