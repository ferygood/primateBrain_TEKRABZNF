library(twice)
library(cowplot)
library(ggpubr)
library(RCy3)
young_kznf <- kznf_infer %>% filter(age == "young")

HmPtC1 <- readRDS("/path/to/results_rdata/HmPtC1.rds")

c1_gene_res <- HmPtC1$DEobject$gene_res %>%
    data.frame() %>%
    filter(log2FoldChange <= -1.5 & pvalue<0.05) %>%
    filter(rownames(.) %in% kznf_infer$external_gene_name)

c1_te_res <- HmPtC1$DEobject$te_res %>%
    data.frame() %>%
    filter(log2FoldChange >= 1.5 & pvalue<0.05)

hm_corr <- HmPtC1$corrRef %>%
    filter(pvalue<0.05 & coef<0) %>%
    filter(geneName %in% rownames(c1_gene_res) &
               teName %in% rownames(c1_te_res)) %>%
    mutate(pair = paste0(geneName, "-", teName))

pt_corr <- HmPtC1$corrCompare %>%
    filter(pvalue<0.05 & coef<0) %>%
    filter(geneName %in% rownames(c1_gene_res) &
               teName %in% rownames(c1_te_res)) %>%
    mutate(pair = paste0(geneName, "-", teName))

intersect(hm_corr$pair, pt_corr$pair) #12

hmGene <- HmPtC1$HmPtC1_DE$geneCorrInputRef["ZNF224", ]
hmTE <- HmPtC1$HmPtC1_DE$teCorrInputRef["LTR3B", ]
ptGene <- HmPtC1$HmPtC1_DE$geneCorrInputCompare["ZNF224", ]
ptTE <- HmPtC1$HmPtC1_DE$teCorrInputCompare["LTR3B", ]

df_input <- data.frame(
    ZNF224 = c(t(hmGene), t(ptGene)),
    LTR3B = c(t(hmTE), t(ptTE)),
    species = factor(c(rep("Human", 40), rep("Chimpanzee", 30)),
                     levels=c("Human", "Chimpanzee"))
)

znf224 <- ggviolin(df_input, x="species", y="ZNF224", fill="species",
                   palette = c("red", "blue"), add="boxplot", add.params = list(fill="white", width=0.1)) +
    stat_compare_means(label="p.signif", label.x=1.4) +
    xlab("") +
    ylab("") +
    theme(legend.position = "right", axis.text.y = element_blank(), axis.title.y=element_blank()) +
    rotate()

ltr3b <-  ggviolin(df_input, x="species", y="LTR3B", fill="species",
                   palette = c("red", "blue"), add="boxplot", add.params = list(fill="white", width=0.1))+
    stat_compare_means(label="p.signif", label.x=1.4) +
    ylab("") +
    theme(legend.position = "right", axis.text.x = element_blank(), axis.title.x=element_blank())

gs <- ggscatter(df_input, x="ZNF224", y="LTR3B",
                color="species", palette = "jco") +
    stat_cor(aes(color=species), label.x=100) +
    scale_color_manual(values = c("red", "blue")) + border() +
    expand_limits(x = c(0, 300), y = c(0, 250))

gs <- gs + rremove("legend")
xplot <- znf224 + rremove("legend")
yplot <- ltr3b

gf <- plot_grid(xplot, NULL, gs, yplot, ncol=2, align = "hv", rel_widths = c(2,2), rel_heights = c(1.5,2))

