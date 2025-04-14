library(introdataviz)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)

exp <- readRDS("~/Desktop/phd/pBrain/data/primateBrain_kznfTEexp.rds")

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

te_exp_longer <- te_exp_longer %>%
    mutate(age = ifelse(TEs %in% te_infer$NM, "young", "old")) %>%
    mutate(value = log(value + 1))

# plot the other region
gv1 <- split_violin_pbd("Primary & Secondary Cortices")
gv2 <- split_violin_pbd("Limbic & Association Cortices")
gv3 <- split_violin_pbd("Archicortex")
gv4 <- split_violin_pbd("Thalamus & Hypothalamus")
gv5 <- split_violin_pbd("Cerebellar White Matter")
gv6 <- split_violin_pbd("Cerebellar Grey Matter")
gv7 <- split_violin_pbd("Striatum")

g_arr <- ggarrange(gv2, gv3, gv4, gv5, gv6, gv7,
                   ncol=3, nrow=2, common.legend = TRUE, legend="bottom")

ggsave(gv1, filename="figures/JPG/2C_kznfs_TEs_violin_cortex.jpg", dpi=500, width=8, height=10)
ggsave(gv1, filename="figures/SVG/2C_kznfs_TEs_violin_cortex.svg", dpi=500, width=8, height=10)

ggsave(g_arr, filename="figures/JPG/S3_split_violin_pbd.jpg", dpi=300,
       width=20, height=20)
