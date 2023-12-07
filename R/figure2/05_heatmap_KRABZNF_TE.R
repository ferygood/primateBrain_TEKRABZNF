# Figure 2E Heatmap

library(ComplexHeatmap)
library(ggplot2)
library(ggpubr)

####################
#### KRAB-ZNFs #####
####################

#1. check human
c1_kznfs <- unique(HmPtC1$corrRef$geneName) #285
c1_kznfs_res <- HmPtC1$DEobject$gene_res %>%
    data.frame() %>%
    filter(abs(log2FoldChange) >= 1.5 & pvalue < 0.05)

select_kznfs <- c1_kznfs_res %>% filter(rownames(.) %in% c1_kznfs) #53

#2. check bonobo
c1_res_HmPp <- HmPpC1$DEobject$gene_res %>%
    data.frame() %>%
    filter(abs(log2FoldChange) >= 1.5 & pvalue < 0.05) %>%
    filter(rownames(.) %in% rownames(select_kznfs)) #24

#3. check macaque
c1_res_HmMm <- HmMmC1$HmMmC1_DE$gene_res %>%
    data.frame() %>%
    filter(abs(log2FoldChange) >= 1.5 & pvalue < 0.05) %>%
    filter(rownames(.) %in% rownames(select_kznfs)) #20

#4. check intersect
union_result <- union(union(rownames(select_kznfs), rownames(c1_res_HmPp)), rownames(c1_res_HmMm)) #53

#5. load expression files
df <- readRDS("~/github/pBrain/data/primateBrain_kznfTEexp.rds")

df_select <- df %>%
    filter(region=="Primary & Secondary Cortices")
df_select_kznfs <- df_select %>%
    select(c(all_of(union_result), "species"))
df_select_kznfs$species <- factor(
    df_select_kznfs$species,
    levels=c("Human", "Chimpanzee", "Bonobo", "Macaque"))

mean_values <- df_select_kznfs %>%
    group_by(species) %>%
    summarize(across(everything(), mean))

df_heatmap_input <- data.frame(t(mean_values[,-1]))
colnames(df_heatmap_input) <- c("Human", "Chimpanzee", "Bonobo", "Macaque")

df_heatmap_input[df_heatmap_input==0] <- NA

df_heatmap_input_branch <- df_heatmap_input %>%
    mutate(name = rownames(.)) %>%
    left_join(kznf_infer[,c(2,6)], join_by(name==external_gene_name)) %>%
    arrange(desc(age))
rownames(df_heatmap_input_branch) <- df_heatmap_input_branch$name

df_check <- df_heatmap_input_branch[,c(1,2,3,4)]
df_check <- df_check[!is.na(df_check$Human),]
result <- df_check[!(rowSums(!is.na(df_check[-1]) & (df_check[-1] > df_check$Human) == 0)) | !(rowSums(!is.na(df_check[-1]) & (df_check[-1] < df_check$Human) == 0)) ,] #36

hm_specific_kznf <- rownames(result)

df_heatmap_input_branch_hm <- df_heatmap_input_branch[hm_specific_kznf,]

#6. plot kznfs heatmap
png("figures/DE_ZNF_heatmap_all.jpg", width=6, height=13, units="in", res=1200)
age_ha <- rowAnnotation(
    age = c(rep("young", 13), rep("old", 23)),
    col = list(age = c("young"="#E59E00", "old"="#5A8FBB")),
    show_annotation_name = FALSE
)

hc1_kznfs <- ComplexHeatmap::Heatmap(
    log(df_heatmap_input_branch_hm[,c(1:4)] + 1),
    right_annotation = age_ha,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    name = "Log Expr.",
    width = 2 * unit(8, "mm"),
    height = 12 * unit(12, "mm"),
    column_names_rot = 45,
    row_names_gp = gpar(fontsize=8),
    column_names_gp = gpar(fontsize=8),
    rect_gp = gpar(col = "black", lwd = 0.8),
    na_col="#706f6f"
)
hc1_kznfs
dev.off()


#############
#### TE #####
#############

# check intersect and te age annotation

HmPt_te_res <- HmPtC1$DEobject$te_res %>%
    data.frame() %>%
    filter(abs(log2FoldChange) >= 1.5 & pvalue < 0.05)

select_TE <- rownames(HmPt_te_res) #53

HmPp_te <- HmPpC1$DEobject$te_res %>%
    data.frame() %>%
    filter(abs(log2FoldChange) >= 1.5 & pvalue < 0.05)  %>%
    filter(rownames(.) %in% select_TE)

HmMm_te <- HmMmC1$HmMmC1_DE$te_res %>%
    data.frame() %>%
    filter(abs(log2FoldChange) >= 1.5 & pvalue < 0.05)  %>%
    filter(rownames(.) %in% select_TE)

intersect_te <- intersect(intersect(select_TE, rownames(HmPp_te)), rownames(HmMm_te))


dfhm <- HmPtC1$DEobject$teCorrInputRef %>%
    mutate(Human = rowMeans(.)) %>%
    select(41) %>%
    filter(rownames(.) %in% intersect_te)

dfchimp <- HmPtC1$DEobject$teCorrInputCompare %>%
    mutate(Chimpanzee = rowMeans(.)) %>%
    select(31) %>%
    filter(rownames(.) %in% intersect_te)

dfbonobo <- HmPpC1$DEobject$teCorrInputCompare %>%
    mutate(Bonobo = rowMeans(.)) %>%
    select(31) %>%
    filter(rownames(.) %in% intersect_te)

dfmacaque <- HmMmC1$HmMmC1_DE$teCorrInputCompare %>%
    mutate(Macaque = rowMeans(.)) %>%
    select(30) %>%
    filter(rownames(.) %in% intersect_te)

dfTE <- cbind(dfhm, dfchimp, dfbonobo, dfmacaque)


# add TE infer age
df_heatmap_TE <- dfTE %>%
    mutate(name = rownames(.)) %>%
    mutate(age = ifelse(name %in% te_infer$NM, "young", "old")) %>%
    arrange(desc(age))

# plot TE heatmap
png("figures/DE_TE_heatmap_all.jpg", width=6, height=13, units="in", res=1200)
age_ha <- rowAnnotation(
    age = c(rep("young", 6), rep("old", 12)),
    col = list(age = c("young"="#E59E00", "old"="#5A8FBB")),
    show_annotation_name = FALSE
)

hc1_TEs <- ComplexHeatmap::Heatmap(
    log(df_heatmap_TE[,c(1:4)]),
    right_annotation = age_ha,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    name = "Log Expr.",
    width = 2 * unit(8, "mm"),
    height = 12 * unit(6, "mm"),
    column_names_rot = 45,
    row_names_gp = gpar(fontsize=8),
    column_names_gp = gpar(fontsize=8),
    rect_gp = gpar(col = "black", lwd = 0.8),
    na_col="#706f6f"
)
hc1_TEs
dev.off()
