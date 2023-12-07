# Figure 2C

library(twice)
library(ggpubr)
library(ggplot2)
library(tidyverse)

data("hmKZNFs337")
data("hg19rmsk_info")

## see young/old KZNFs, young/old TEs expression in different brain regions comparing species
exp <- readRDS("~/github/pBrain/data/primateBrain_kznfTEexp.rds")

# split gene and TE dataframe
kznf_exp <- exp[,c(1:285, 1191, 1192)]
kznf_exp_longer <- kznf_exp %>%
    mutate(sample = rownames(.)) %>%
    pivot_longer(cols = c(1:285), names_to = "KRAB-ZNFs",values_to = "value")

te_exp <- exp[,c(286:1192)]
te_exp_longer <- te_exp %>%
    mutate(sample = rownames(.)) %>%
    pivot_longer(cols = c(1:905), names_to = "TEs", values_to = "value")

# add age inference
kznf_exp_longer <- kznf_exp_longer %>%
    left_join(kznf_infer[,c(2,6)], by=c("KRAB-ZNFs"="external_gene_name")) %>%
    mutate(value = log(value + 1))

te_infer <- read.csv("~/github/pBrain/data/Dfam_TE_simiiformes.csv")
te_exp_longer <- te_exp_longer %>%
    mutate(age = ifelse(TEs %in% te_infer$NM, "young", "old")) %>%
    mutate(value = log(value + 1))

compare_violin <- function(kznfs, TEs, brain_region){
    c_kznfs <- kznfs %>%
        filter(region == brain_region)

    c_kznfs$species <- factor(c_kznfs$species,
                              levels=c("Human", "Chimpanzee", "Bonobo", "Macaque"))

    gv <- ggviolin(c_kznfs, x="age", y="value", fill="age",
                   palette=c("#5A8FBB", "#E59E00"),
                   add = "boxplot", add.params = list(fill="white")) +
        stat_compare_means(label="p.signif", label.x=1.5, method="wilcox.test") +
        theme(legend.position = "bottom", axis.text.x=element_blank()) +
        facet_wrap(~species, ncol=4) +
        ylab("Log expression") +
        xlab("") +
        ggtitle(paste0("KRAB-ZNFs in ", brain_region))

    #ggsave(filename=gfile, gv, dpi=500, width=10, height=4)

    c_te <- TEs %>%
        filter(region == brain_region)

    c_te$species <- factor(c_te$species,
                           levels=c("Human", "Chimpanzee", "Bonobo", "Macaque"))
    tv <- ggviolin(c_te, x="age", y="value", fill="age",
                   palette=c("#5A8FBB", "#E59E00"),
                   add = "boxplot", add.params = list(fill="white")) +
        stat_compare_means(label="p.signif", label.x=1.5, method="wilcox.test") +
        theme(legend.position = "bottom", axis.text.x=element_blank()) +
        facet_wrap(~species, ncol=4) +
        ylab("Log expression") +
        xlab("") +
        ggtitle(paste0("TEs in ", brain_region))

    #ggsave(filename = tfile, tv, dpi=500, width=10, height=4)

    result <- list("gv"=gv, "tv"=tv)
    result
}

c1_v <- compare_violin(kznf_exp_longer, te_exp_longer, "Primary & Secondary Cortices")
c2_v <- compare_violin(kznf_exp_longer, te_exp_longer, "Limbic & Association Cortices")
c3_v <- compare_violin(kznf_exp_longer, te_exp_longer, "Archicortex")
c4_v <- compare_violin(kznf_exp_longer, te_exp_longer, "Thalamus & Hypothalamus")
c5_v <- compare_violin(kznf_exp_longer, te_exp_longer, "Cerebellar White Matter")
c6_v <- compare_violin(kznf_exp_longer, te_exp_longer, "Cerebellar Grey Matter")
c7_v <- compare_violin(kznf_exp_longer, te_exp_longer, "Striatum")

# combine figures
g_all <- ggarrange(c1_v$gv, c2_v$gv, c3_v$gv, c4_v$gv, c5_v$gv, c6_v$gv, c7_v$gv, nrow = 7, common.legend = TRUE, legend = "bottom")
ggsave(filename="figures/kznfs_TEs_violin_brain/kznfs.jpg", dpi=500, width=10, height=28)

t_all <- ggarrange(c1_v$tv, c2_v$tv, c3_v$tv, c4_v$tv, c5_v$tv, c6_v$tv, c7_v$tv, nrow = 7, common.legend = TRUE, legend = "bottom")
ggsave(filename="figures/kznfs_TEs_violin_brain/TEs.jpg", dpi=500, width=10, height=28)

