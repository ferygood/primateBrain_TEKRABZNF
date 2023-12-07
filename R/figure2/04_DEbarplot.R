# Figure 2D

library(ComplexHeatmap)
library(ggplot2)
library(ggpubr)

c1_kznfs <- unique(HmPtC1$corrRef$geneName) #285
c1_kznfs_res <- HmPtC1$DEobject$gene_res %>%
    data.frame() %>%
    filter(log2FoldChange <= -1.5 & pvalue < 0.05)

select_kznfs <- c1_kznfs_res %>% filter(rownames(.) %in% c1_kznfs) #12

#2. check bonobo
c1_res_HmPp <- HmPpC1$DEobject$gene_res %>%
    data.frame() %>%
    filter(log2FoldChange <= -1.5 & pvalue < 0.05) %>%
    filter(rownames(.) %in% rownames(select_kznfs)) #2

#3. check macaque
c1_res_HmMm <- HmMmC1$HmMmC1_DE$gene_res %>%
    data.frame() %>%
    filter(log2FoldChange <= -1.5 & pvalue < 0.05) %>%
    filter(rownames(.) %in% rownames(select_kznfs)) #3

#4. check intersect
union_result <- union(union(rownames(select_kznfs), rownames(c1_res_HmPp)), rownames(c1_res_HmMm))

# krab-znfs
c1_de_kznfs <- HmPtC1$DEobject$gene_res %>%
    data.frame() %>%
    mutate(geneName = rownames(.)) %>%
    inner_join(kznf_infer[,c(2,6)], join_by(geneName==external_gene_name)) %>%
    mutate(DE = ifelse((abs(log2FoldChange) >= 1.5 & pvalue < 0.05), 1, 0))

de_young_kznfs <- c1_de_kznfs %>%
    filter(DE==1 & age=="young") #21
de_old_kznfs <- c1_de_kznfs %>%
    filter(DE==1 & age=="old") #32


# TEs
c1_de_TEs <- HmPtC1$DEobject$te_res %>%
    data.frame() %>%
    mutate(teName = rownames(.)) %>%
    mutate(age = ifelse(teName %in% te_infer$NM, "young", "old")) %>%
    mutate(DE = ifelse((abs(log2FoldChange) >= 1.5 & pvalue < 0.05), 1, 0))

de_young_TEs <- c1_de_TEs %>%
    filter(DE==1 & age=="young") #23
de_old_TEs <- c1_de_TEs %>%
    filter(DE==1 & age=="old") #30


# create dataframe
df <- data.frame(
    raw_count = c(32, 21, 30, 23),
    scale = c(237, 100, 842, 309),
    age = c("old", "young",
            "old","young"),
    name = c("DE KRAB-ZNFs", "DE KRAB-ZNFs", "DE TEs", "DE TEs")
)

df <- df %>%
    mutate(normalized_count = raw_count*100/scale)

gc1 <- ggplot(df, aes(x=name, y=normalized_count, fill=age)) +
    geom_bar(stat='identity', position=position_dodge()) +
    scale_fill_manual(values = c("#5A8FBB", "#E59E00")) +
    ylab("Percentage counts (%)") +
    xlab("") +
    ggtitle("Primary & Secondary cortices") +
    theme_bw()

# in cerebellum white matter
# krab-znfs
c5_de_kznfs <- HmPtC5$DEobject$gene_res %>%
    data.frame() %>%
    mutate(geneName = rownames(.)) %>%
    inner_join(kznf_infer[,c(2,6)], join_by(geneName==external_gene_name)) %>%
    mutate(DE = ifelse((abs(log2FoldChange) >= 1.5 & pvalue < 0.05), 1, 0))

de_young_kznfs_c5 <- c5_de_kznfs %>%
    filter(DE==1 & age=="young") #22
de_old_kznfs_c5 <- c5_de_kznfs %>%
    filter(DE==1 & age=="old") #32

# TEs
c5_de_TEs <- HmPtC5$DEobject$te_res %>%
    data.frame() %>%
    mutate(teName = rownames(.)) %>%
    mutate(age = ifelse(teName %in% te_infer$NM, "young", "old")) %>%
    mutate(DE = ifelse((abs(log2FoldChange) >= 1.5 & pvalue < 0.05), 1, 0))

de_young_TEs_c5 <- c5_de_TEs %>%
    filter(DE==1 & age=="young") #19
de_old_TEs_c5 <- c5_de_TEs %>%
    filter(DE==1 & age=="old") #25

# create dataframe
df_c5 <- data.frame(
    raw_count = c(32, 22, 25, 19),
    scale = c(205, 78, 817, 290),
    age = c("old", "young",
            "old","young"),
    name = c("DE KRAB-ZNFs", "DE KRAB-ZNFs", "DE TEs", "DE TEs")
)

df_c5 <- df_c5 %>%
    mutate(normalized_count = raw_count*100/scale)

gc2 <- ggplot(df_c5, aes(x=name, y=normalized_count, fill=age)) +
    geom_bar(stat='identity', position=position_dodge()) +
    scale_fill_manual(values = c("#5A8FBB", "#E59E00")) +
    ylab("Percentage counts (%)") +
    xlab("") +
    ggtitle("Cerebellar White Matter") +
    theme_bw()

gc <- ggarrange(gc1, gc2, ncol=1, common.legend = TRUE)
ggsave(filename="figures/DE_KZNFS_TEs_percentage.jpg", gc, dpi=400, width=3, height=8)
