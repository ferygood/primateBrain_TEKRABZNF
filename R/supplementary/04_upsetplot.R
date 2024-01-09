library(dplyr)
library(ggvenn)
library(ggpubr)
library(ComplexHeatmap)

# in primary and secondary cortices
HmC1 <- HmPtC1$corrRef %>% filter(pvalue < 0.05 & coef < 0) %>%
    mutate(link = paste0(geneName, "-", teName))

PtC1 <- HmPtC1$corrCompare %>% filter(pvalue < 0.05 & coef < 0) %>%
    mutate(link = paste0(geneName, "-", teName))

PpC1 <- HmPpC1$corrCompare %>% filter(pvalue < 0.05 & coef < 0) %>%
    mutate(link = paste0(geneName, "-", teName))

MmC1 <- HmMmC1$HmMmC1_MmCorr %>% filter(pvalue < 0.05 & coef < 0) %>%
    mutate(link = paste0(geneName, "-", teName))

a <- list(
    "human" = HmC1$link,
    "chimpanzee" = PtC1$link,
    "bonobo" = PpC1$link,
    "macaque" = MmC1$link
)

g1 <- ggvenn(a, stroke_linetype = 2, stroke_size = 0.5,
             fill_color = c("red", "blue", "purple", "green"),
             set_name_size = 4) +
    ggtitle("All Negative TE:KRAB-ZNF")

young_kznf <- kznf_infer %>%
    filter(age=="young") %>%
    select(2) %>%
    unlist()
young_TEs <- te_infer$NM

# old and young
HmC1_old <- HmC1 %>%
    filter(!(geneName %in% young_kznf & teName %in% young_TEs))
HmC1_young <- HmC1 %>%
    filter(geneName %in% young_kznf & teName %in% young_TEs)

PtC1_old <- PtC1 %>%
    filter(!(geneName %in% young_kznf & teName %in% young_TEs))
PtC1_young <- PtC1 %>%
    filter(geneName %in% young_kznf & teName %in% young_TEs)

PpC1_old <- PpC1 %>%
    filter(!(geneName %in% young_kznf & teName %in% young_TEs))
PpC1_young <- PpC1 %>%
    filter(geneName %in% young_kznf & teName %in% young_TEs)

MmC1_old <- MmC1 %>%
    filter(!(geneName %in% young_kznf & teName %in% young_TEs))
MmC1_young <- MmC1 %>%
    filter(geneName %in% young_kznf & teName %in% young_TEs)

b <- list(
    "human" = HmC1_young$link,
    "chimpanzee" = PtC1_young$link,
    "bonobo" = PpC1_young$link,
    "macaque" = MmC1_young$link
)

g2 <- ggvenn(c, stroke_linetype = 2, stroke_size = 0.5,
             fill_color = c("red", "blue", "purple", "green"),
             set_name_size = 4) +
    ggtitle("Young Negative TE:KRAB-ZNF")

jpeg(filename = "figures/upset_all.jpg", width=4, height=3, units="in", res=1200, pointsize = 8)
# all
a_all <- make_comb_mat(a)
UpSet(a_all, set_order = c("human", "chimpanzee", "bonobo", "macaque"), comb_order = order(desc(comb_size(a_all))))

dev.off()

# young
jpeg(filename = "figures/upset_young.jpg", width=4, height=3, units="in", res=1200, pointsize = 8)
b_all <- make_comb_mat(b)
UpSet(b_all, set_order = c("human", "chimpanzee", "bonobo", "macaque"), comb_order = order(desc(comb_size(b_all))))
dev.off()

