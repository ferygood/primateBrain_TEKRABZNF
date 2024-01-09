library(ComplexHeatmap)
library(ggplot2)

young_kznf <- kznf_infer %>% filter(age=="young")

# preprocess tcx
tcx_res_kznfs <- mayoTEKRABber$tcxDE$gene_res %>%
    data.frame() %>%
    filter(abs(log2FoldChange) >= 0.5 & pvalue < 0.05) %>%
    filter(rownames(.) %in% kznf_infer$external_gene_name)

tcx_exp_kznfs <- mayoTEKRABber$tcxDE$normalized_gene_counts %>%
    data.frame() %>%
    filter(rownames(.) %in% rownames(tcx_res_kznfs)) %>%
    mutate(tcx_control = rowMeans(.[1:23])) %>%
    mutate(tcx_AD = rowMeans(.[24:47])) %>%
    select(c(48, 49)) %>%
    mutate(branch = ifelse(rownames(.) %in% young_kznf$external_gene_name,
                           "young", "old")) %>%
    arrange(desc(branch))

# preprocess cbe
cbe_res_kznfs <- mayoTEKRABber$cbeDE$gene_res %>%
    data.frame() %>%
    filter(abs(log2FoldChange) >= 0.5 & pvalue < 0.05) %>%
    filter(rownames(.) %in% kznf_infer$external_gene_name)

cbe_exp_kznfs <- mayoTEKRABber$cbeDE$normalized_gene_counts %>%
    data.frame() %>%
    filter(rownames(.) %in% rownames(cbe_res_kznfs)) %>%
    mutate(cbe_control = rowMeans(.[1:23])) %>%
    mutate(cbe_AD = rowMeans(.[24:45])) %>%
    select(c(46, 47)) %>%
    mutate(branch = ifelse(rownames(.) %in% young_kznf$external_gene_name,
                           "young", "old")) %>%
    arrange(desc(branch))

png("figures/tcx_kznfs_heatmap.jpg", width=6, height=10, units="in", res=1200)
age_ha <- rowAnnotation(
    age = c(rep("young", 2), rep("old", 4)),
    col = list(age = c("young" = "#E59E00", "old" = "#5A8FBB")),
    show_annotation_name = FALSE
)


tcx_kznfs_heatmap <- ComplexHeatmap::Heatmap(
    log(tcx_exp_kznfs[,c(1,2)] + 1),
    right_annotation = age_ha,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    name = "Log Expr.",
    width = 2 * unit(5, "mm"),
    height = 6 * unit(5, "mm"),
    column_names_rot = 45,
    row_names_gp = gpar(fontsize=8),
    column_names_gp = gpar(fontsize=8),
    rect_gp = gpar(col = "black", lwd = 0.8),
    na_col="#dbd9d9"
)

tcx_kznfs_heatmap
dev.off()

png("figures/cbe_kznfs_heatmap.jpg", width=6, height=10, units="in", res=1200)
age_ha <- rowAnnotation(
    age = c(rep("young", 8), rep("old", 6)),
    col = list(age = c("young" = "#E59E00", "old" = "#5A8FBB")),
    show_annotation_name = FALSE
)

cbe_kznfs_heatmap <- ComplexHeatmap::Heatmap(
    log(cbe_exp_kznfs[,c(1,2)] + 1),
    right_annotation = age_ha,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    name = "Log Expr.",
    width = 2 * unit(5, "mm"),
    height = 12 * unit(5, "mm"),
    column_names_rot = 45,
    row_names_gp = gpar(fontsize=8),
    column_names_gp = gpar(fontsize=8),
    rect_gp = gpar(col = "black", lwd = 0.8),
    na_col="#dbd9d9"
)

cbe_kznfs_heatmap
dev.off()

# TEs
# preprocess tcx
tcx_res_TEs <- mayoTEKRABber$tcxDE$te_res %>%
    data.frame() %>%
    filter(abs(log2FoldChange) >= 0.5 & pvalue < 0.05)

tcx_exp_TEs <- mayoTEKRABber$tcxDE$normalized_te_counts %>%
    data.frame() %>%
    filter(rownames(.) %in% rownames(tcx_res_TEs)) %>%
    mutate(tcx_control = rowMeans(.[1:23])) %>%
    mutate(tcx_AD = rowMeans(.[24:47])) %>%
    select(c(48, 49)) %>%
    mutate(branch = ifelse(rownames(.) %in% te_infer$NM,
                           "young", "old")) %>%
    arrange(desc(branch))

# preprocess cbe
cbe_res_TEs <- mayoTEKRABber$cbeDE$te_res %>%
    data.frame() %>%
    filter(abs(log2FoldChange) >= 0.5 & pvalue < 0.05)

cbe_exp_TEs <- mayoTEKRABber$cbeDE$normalized_te_counts %>%
    data.frame() %>%
    filter(rownames(.) %in% rownames(cbe_res_TEs)) %>%
    mutate(cbe_control = rowMeans(.[1:23])) %>%
    mutate(cbe_AD = rowMeans(.[24:45])) %>%
    select(c(46, 47)) %>%
    mutate(branch = ifelse(rownames(.) %in% te_infer$NM,
                           "young", "old")) %>%
    arrange(desc(branch))

png("figures/tcx_TEs_heatmap.jpg", width=6, height=10, units="in", res=1200)
age_ha <- rowAnnotation(
    age = c(rep("young", 5), rep("old", 5)),
    col = list(age = c("young" = "#E59E00", "old" = "#5A8FBB")),
    show_annotation_name = FALSE
)

tcx_TEs_heatmap <- ComplexHeatmap::Heatmap(
    log(tcx_exp_TEs[,c(1,2)] + 1),
    right_annotation = age_ha,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    name = "Log Expr.",
    width = 2 * unit(5, "mm"),
    height = 10 * unit(5, "mm"),
    column_names_rot = 45,
    row_names_gp = gpar(fontsize=8),
    column_names_gp = gpar(fontsize=8),
    rect_gp = gpar(col = "black", lwd = 0.8),
    na_col="#dbd9d9"
)
tcx_TEs_heatmap
dev.off()

png("figures/cbe_TEs_heatmap.jpg", width=6, height=10, units="in", res=1200)

age_ha <- rowAnnotation(
    age = c(rep("young", 2), rep("old", 5)),
    col = list(age = c("young" = "#E59E00", "old" = "#5A8FBB")),
    show_annotation_name = FALSE
)

cbe_TEs_heatmap <- ComplexHeatmap::Heatmap(
    log(cbe_exp_TEs[,c(1,2)] + 1),
    right_annotation = age_ha,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    name = "Log Expr.",
    width = 2 * unit(5, "mm"),
    height = 7 * unit(5, "mm"),
    column_names_rot = 45,
    row_names_gp = gpar(fontsize=8),
    column_names_gp = gpar(fontsize=8),
    rect_gp = gpar(col = "black", lwd = 0.8),
    na_col="#dbd9d9"
)

cbe_TEs_heatmap
