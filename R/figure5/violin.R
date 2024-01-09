library(twice)
library(dplyr)
library(ggplot2)
library(ggpubr)

load("/path/to/mayoTEKRABber_balance.RData")
kznf_age_infer <- read.csv("/path/to/hmKZNFs337_ageInfer.csv")
te_infer <- read.csv("/path/to/Dfam_TE_simiiformes.csv")


# cerebellum expression data load
cbe_gene <- mayoTEKRABber$cbeDE$normalized_gene_counts %>%
    data.frame()
cbe_kznfs <- cbe_gene[kznf_age_infer$external_gene_name, ]

cbe_TE <- mayoTEKRABber$cbeDE$normalized_te_counts %>%
    data.frame()


cbe_kznfs <- data.frame(t(cbe_kznfs))

df_cbe_kznfs <- cbe_kznfs %>%
    mutate(group=c(rep("control", 23), rep("AD", 22)))

df_cbe_kznfs <- df_cbe_kznfs[, colSums(is.na(df_cbe_kznfs)) < nrow(df_cbe_kznfs)]

df_cbe_kznfs <- df_cbe_kznfs %>%
    pivot_longer(cols = c(1:330), names_to = "KRAB_ZNFs", values_to = "expression") %>%
    left_join(kznf_infer[,c(2,6)], join_by(KRAB_ZNFs==external_gene_name)) %>%
    mutate(expression=log(expression))

cbe_kznfs <- ggviolin(df_cbe_kznfs, x="age", y="expression", fill="age",
                      palette = c("#5A8FBB", "#E59E00"), add="boxplot",
                      add.params = list(fill="white")) +
    stat_compare_means(label="p.signif", label.x=1.5, method="wilcox.test") +
    theme(legend.position = "right",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ylab("Log expression") +
    ggtitle("KRAB-ZNFs in cerebellumn (AD vs. Control)") +
    facet_wrap(~group)

cbe_TE <- data.frame(t(cbe_TE))

df_cbe_TE <- cbe_TE %>%
    mutate(group=c(rep("control", 23), rep("AD", 22)))

df_cbe_TE <- df_cbe_TE[, colSums(is.na(df_cbe_TE)) < nrow(df_cbe_TE)]

df_cbe_TE <- df_cbe_TE %>%
    pivot_longer(cols = c(1:1048), names_to = "TEs", values_to = "expression") %>%
    mutate(age = ifelse(TEs %in% te_infer$NM, "young", "old")) %>%
    mutate(expression = log(expression))

cbe_TE <- ggviolin(df_cbe_TE, x="age", y="expression", fill="age",
                   palette = c("#5A8FBB", "#E59E00"), add="boxplot",
                   add.params = list(fill="white")) +
    stat_compare_means(label="p.signif", label.x=1.5, method="wilcox.test") +
    theme(legend.position = "right",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ylab("Log expression") +
    ggtitle("TEs in cerebellumn (AD vs. Control)") +
    facet_wrap(~group)

# load expression files
tcx_gene <- mayoTEKRABber$tcxDE$normalized_gene_counts %>%
    data.frame()
tcx_kznfs <- tcx_gene[kznf_age_infer$external_gene_name, ]

tcx_TE <- mayoTEKRABber$tcxDE$normalized_te_counts %>%
    data.frame()

# select kznfs
tcx_kznfs <- data.frame(t(tcx_kznfs))

df_tcx_kznfs <- tcx_kznfs %>%
    mutate(group=c(rep("control", 23), rep("AD", 24)))

df_tcx_kznfs <- df_tcx_kznfs[, colSums(is.na(df_tcx_kznfs)) < nrow(df_tcx_kznfs)]

df_tcx_kznfs <- df_tcx_kznfs %>%
    pivot_longer(cols = c(1:330), names_to = "KRAB_ZNFs", values_to = "expression") %>%
    left_join(kznf_infer[,c(2,6)], join_by(KRAB_ZNFs==external_gene_name)) %>%
    mutate(expression=log(expression))

# select TEs
tcx_TE <- data.frame(t(tcx_TE))

df_tcx_TE <- tcx_TE %>%
    mutate(group=c(rep("control", 23), rep("AD", 24)))

df_tcx_TE <- df_tcx_TE[, colSums(is.na(df_tcx_TE)) < nrow(df_tcx_TE)]

df_tcx_TE <- df_tcx_TE %>%
    pivot_longer(cols = c(1:1048), names_to = "TEs", values_to = "expression") %>%
    mutate(age = ifelse(TEs %in% te_infer$NM, "young", "old")) %>%
    mutate(expression = log(expression))

tcx_kznfs <- ggviolin(df_tcx_kznfs, x="age", y="expression", fill="age",
                      palette = c("#5A8FBB", "#E59E00"), add="boxplot",
                      add.params = list(fill="white")) +
    stat_compare_means(label="p.signif", label.x=1.5, method="wilcox.test") +
    theme(legend.position = "right",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ylab("Log expression") +
    ggtitle("KRAB-ZNFs in temporal cortex (AD vs. Control)") +
    facet_wrap(~group)

tcx_TE <- ggviolin(df_tcx_TE, x="age", y="expression", fill="age",
                   palette = c("#5A8FBB", "#E59E00"), add="boxplot",
                   add.params = list(fill="white")) +
    stat_compare_means(label="p.signif", label.x=1.5, method="wilcox.test") +
    theme(legend.position = "right",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ylab("Log expression") +
    ggtitle("TEs in temporal cortex (AD vs. Control)") +
    facet_wrap(~group)

g <- ggarrange(cbe_kznfs, cbe_TE, tcx_kznfs, tcx_TE, ncol=2, nrow=2, common.legend = TRUE)

# temporal cortex
df_tcx_kznfs_select <- df_tcx_kznfs %>%
    select(c(1,3,4)) %>%
    mutate(type="KRAB-ZNFs")

df_tcx_TEs_select <- df_tcx_TE %>%
    select(c(1,3,4)) %>%
    mutate(type="TEs")

df_tcx <- rbind(df_tcx_kznfs_select, df_tcx_TEs_select)

df_tcx$group <- factor(df_tcx$group, levels = c("control", "AD"))

df_tcx <- filter(df_tcx, !is.infinite(expression))

tcx_v <- ggplot(df_tcx, aes(x=group, y=expression, fill=age)) +
    introdataviz::geom_split_violin(alpha=.8)+
    geom_boxplot(width = .2, alpha = .6, show.legend = FALSE) +
    facet_wrap(vars(type)) +
    stat_compare_means(label = "p.signif", label.x=1.5, method="wilcox.test") +
    scale_fill_manual(values = c("#5A8FBB", "#E59E00")) +
    labs(fill="evolutionary age") +
    theme_bw() +
    rotate_x_text(angle = 20) +
    xlab("") +
    ylab("logExp.") +
    ggtitle("Expression in temporal cortex") +
    theme(text = element_text(size = 16))

df_cbe_kznfs_select <- df_cbe_kznfs %>%
    select(c(1,3,4)) %>%
    mutate(type="KRAB-ZNFs")

df_cbe_TEs_select <- df_cbe_TE %>%
    select(c(1,3,4)) %>%
    mutate(type="TEs")

df_cbe <- rbind(df_cbe_kznfs_select, df_cbe_TEs_select)

df_cbe$group <- factor(df_cbe$group, levels = c("control", "AD"))

df_cbe <- filter(df_cbe, !is.infinite(expression))

cbe_v <- ggplot(df_cbe, aes(x=group, y=expression, fill=age)) +
    introdataviz::geom_split_violin(alpha=.8)+
    geom_boxplot(width = .2, alpha = .6, show.legend = FALSE) +
    facet_wrap(vars(type)) +
    stat_compare_means(label = "p.signif", label.x=1.5, method="wilcox.test") +
    scale_fill_manual(values = c("#5A8FBB", "#E59E00")) +
    labs(fill="evolutionary age") +
    theme_bw() +
    rotate_x_text(angle = 20) +
    xlab("") +
    ylab("logExp.") +
    ggtitle("Expression in cerebellum") +
    theme(text = element_text(size = 16))

g_split <- ggarrange(tcx_v, cbe_v, ncol=1, nrow=2, common.legend = TRUE, legend="bottom")
