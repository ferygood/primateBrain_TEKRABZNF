library(TEKRABber)
library(dplyr)
library(twice)
library(RCy3)

load("/path/to/tables/primateBrainRMSK.rda")
c1 <- brain_meta %>% filter(cluster=="cluster1")
hmid <- metadata %>%
    filter(Organism=="Homo sapiens" & brain_region %in% c1$region)
ptid <- metadata %>%
    filter(Organism=="Pan troglodytes" & brain_region %in% c1$region)
ppid <- metadata %>%
    filter(Organism=="Pan paniscus" & brain_region %in% c1$region)
mmid <- metadata %>%
    filter(Organism=="Macaca mulatta" & brain_region %in% c1$region)

ptGene_c1 <- ptGene[,c("geneID", ptid$Run)]
ptTE_c1 <- ptTE[,c("name", ptid$Run)]
hmGene_c1 <- hmGene[,c("geneID", hmid$Run)]
hmTE_c1 <- hmTE[,c("name", hmid$Run)]
ppGene_c1 <- ppGene[,c("geneID",ppid$Run)]
ppTE_c1 <- ppTE[,c("name", ppid$Run)]
mmGene_c1 <- mmGene[,c("geneID", mmid$Run)]
mmTE_c1 <- mmTE[,c("name", mmid$Run)]

chimpHm <- orthologScale(
    speciesRef = "ptroglodytes",
    speciesCompare = "hsapiens",
    geneCountRef = ptGene_c1,
    geneCountCompare = hmGene_c1,
    teCountRef = ptTE_c1,
    teCountCompare = hmTE_c1,
    rmsk = primateBrainRMSK[,c(1,2,4,3)]
)

chimpPp <- orthologScale(
    speciesRef = "ptroglodytes",
    speciesCompare = "ppaniscus",
    geneCountRef = ptGene_c1,
    geneCountCompare = ppGene_c1,
    teCountRef = ptTE_c1,
    teCountCompare = ppTE_c1,
    rmsk = primateBrainRMSK[,c(1,2,4,5)]
)

chimpMm <- orthologScale(
    speciesRef = "ptroglodytes",
    speciesCompare = "mmulatta",
    geneCountRef = ptGene_c1,
    geneCountCompare = mmGene_c1,
    teCountRef = ptTE_c1,
    teCountCompare = mmTE_c1,
    rmsk = primateBrainRMSK[,c(1,2,4,6)]
)

inputPtHm <- DECorrInputs(chimpHm)
inputPtPp <- DECorrInputs(chimpPp)
inputPtMm <- DECorrInputs(chimpMm)

# chimpanzee vs human
meta_PtHm <- data.frame(
    species = c(rep("chimpanzee", ncol(chimpHm$geneRef) -1),
                rep("human", ncol(chimpHm$geneCompare) -1))
)

rownames(meta_PtHm) <- colnames(inputPtHm$geneInputDESeq2)
meta_PtHm$species <- factor(meta_PtHm$species,
                            levels=c("chimpanzee", "human"))

PtHmDE <- DEgeneTE(
    geneTable = inputPtHm$geneInputDESeq2,
    teTable = inputPtHm$teInputDESeq2,
    metadata = meta_PtHm,
    expDesign = TRUE
)

# chimpanzee vs bonobo
meta_PtPp <- data.frame(
    species = c(rep("chimpanzee", ncol(chimpPp$geneRef) -1),
                rep("bonobo", ncol(chimpPp$geneCompare) -1))
)

rownames(meta_PtPp) <- colnames(inputPtPp$geneInputDESeq2)
meta_PtPp$species <- factor(meta_PtPp$species,
                            levels=c("chimpanzee", "bonobo"))

PtPpDE <- DEgeneTE(
    geneTable = inputPtPp$geneInputDESeq2,
    teTable = inputPtPp$teInputDESeq2,
    metadata = meta_PtPp,
    expDesign = TRUE
)

# chimpanzee vs macaque
meta_PtMm <- data.frame(
    species = c(rep("chimpanzee", ncol(chimpMm$geneRef) -1),
                rep("macaque", ncol(chimpMm$geneCompare) -1))
)

rownames(meta_PtMm) <- colnames(inputPtMm$geneInputDESeq2)
meta_PtMm$species <- factor(meta_PtMm$species,
                            levels=c("chimpanzee", "macaque"))

PtMmDE <- DEgeneTE(
    geneTable = inputPtMm$geneInputDESeq2,
    teTable = inputPtMm$teInputDESeq2,
    metadata = meta_PtMm,
    expDesign = TRUE
)


# KRAB-ZNFs intersect
de_PtHm <- data.frame(PtHmDE$gene_res) %>%
    mutate(DE = case_when(
        log2FoldChange < -1.2 & pvalue < 0.05 ~ "up",
        log2FoldChange > 1.2 & pvalue < 0.05 ~ "down",
        TRUE ~ NA_character_
    )) %>%
    filter(!is.na(DE)) %>%
    mutate(gene = rownames(.)) %>%
    select(c(7,6))

de_PtPp <- data.frame(PtPpDE$gene_res) %>%
    mutate(DE = case_when(
        log2FoldChange < -1.2 & pvalue < 0.05 ~ "up",
        log2FoldChange > 1.2 & pvalue < 0.05 ~ "down",
        TRUE ~ NA_character_
    )) %>%
    filter(!is.na(DE)) %>%
    mutate(gene = rownames(.)) %>%
    select(c(7,6))


de_PtMm <- data.frame(PtMmDE$gene_res) %>%
    mutate(DE = case_when(
        log2FoldChange < -1.2 & pvalue < 0.05 ~ "up",
        log2FoldChange > 1.2 & pvalue < 0.05 ~ "down",
        TRUE ~ NA_character_
    )) %>%
    filter(!is.na(DE)) %>%
    mutate(gene = rownames(.)) %>%
    select(c(7,6))

de_gene <- inner_join(de_PtHm, de_PtPp, by=c("gene", "DE")) %>%
    inner_join(de_PtMm, by=c("gene", "DE"))

rownames(de_gene) <- de_gene$gene
de_geneName <- ensIDtoGeneName(de_gene, species = "ptroglodytes")
de_geneName <- de_geneName %>%
    mutate(geneName = rownames(.)) %>%
    filter(geneName %in% kznf_infer$external_gene_name) #only five KRAB-ZNFs
# ZNF331, ZNF230, ZNF165, ZNF551, ZNF224

# TEs intersect
deTE_PtHm <- data.frame(PtHmDE$te_res) %>%
    mutate(DE = case_when(
        log2FoldChange < -1.2 & pvalue < 0.05 ~ "up",
        log2FoldChange > 1.2 & pvalue < 0.05 ~ "down",
        TRUE ~ NA_character_
    )) %>%
    filter(!is.na(DE)) %>%
    mutate(te = rownames(.)) %>%
    select(c(7,6))

deTE_PtPp <- data.frame(PtPpDE$te_res) %>%
    mutate(DE = case_when(
        log2FoldChange < -1.2 & pvalue < 0.05 ~ "up",
        log2FoldChange > 1.2 & pvalue < 0.05 ~ "down",
        TRUE ~ NA_character_
    )) %>%
    filter(!is.na(DE)) %>%
    mutate(te = rownames(.)) %>%
    select(c(7,6))


deTE_PtMm <- data.frame(PtMmDE$te_res) %>%
    mutate(DE = case_when(
        log2FoldChange < -1.2 & pvalue < 0.05 ~ "up",
        log2FoldChange > 1.2 & pvalue < 0.05 ~ "down",
        TRUE ~ NA_character_
    )) %>%
    filter(!is.na(DE)) %>%
    mutate(te = rownames(.)) %>%
    select(c(7,6))

de_te <- inner_join(deTE_PtHm, deTE_PtPp, by=c("te", "DE")) %>%
    inner_join(deTE_PtMm, by=c("te", "DE"))

rownames(de_te) <- de_gene$te #only six TEs
# Arthur1C, HERVL32-int, LTR15, LTR39, LTR6A, MER88

kznf_list <- c("ZNF331", "ZNF230", "ZNF165", "ZNF551", "ZNF224")
te_list <- c("Arthur1C", "HERVL32-int", "LTR15", "LTR39", "LTR6A", "MER88")

kznf_id <- de_geneName %>%
    filter(geneName %in% kznf_list)

# calculate correlation
PtCorr <- corrOrthologTE(
    geneInput = PtHmDE$geneCorrInputRef[kznf_id$gene,],
    teInput = PtHmDE$teCorrInputRef[te_list,],
    corrMethod = "pearson",
    padjMethod = "fdr"
)

HmCorr <- corrOrthologTE(
    geneInput = PtHmDE$geneCorrInputCompare[kznf_id$gene,],
    teInput = PtHmDE$teCorrInputCompare[te_list,],
    corrMethod = "pearson",
    padjMethod = "fdr"
)

# select negative significant results
PtCorr_sig <- PtCorr %>%
    filter(pvalue<0.05 & coef<0)

HmCorr_sig <- HmCorr %>%
    filter(pvalue<0.05 & coef<0)

PtCorr_sig_DE <- PtCorr_sig %>%
    left_join(de_geneName, join_by(geneName==gene)) %>%
    left_join(de_te, join_by(teName==te))

node.Pt <- data.frame(id=c(unique(PtCorr_sig_DE$geneName.y),
                           unique(PtCorr_sig_DE$teName)))
node.Pt$DE <- c("up", "down", "up", "up", "up", "up")

link.Pt <- data.frame(source = PtCorr_sig_DE$geneName.y,
                      target = PtCorr_sig_DE$teName)

createNetworkFromDataFrames(node.Pt, link.Pt)

HmCorr_sig_DE <- HmCorr_sig %>%
    left_join(de_geneName, join_by(geneName==gene)) %>%
    left_join(de_te, join_by(teName==te))

node.Hm <- data.frame(
    id=c(unique(HmCorr_sig_DE$geneName.y),
         unique(HmCorr_sig_DE$teName))
)

node.Hm$DE <- c("up", "down", "up", "up", "up")

link.Hm <- data.frame(source = HmCorr_sig_DE$geneName.y,
                      target = HmCorr_sig_DE$teName)

createNetworkFromDataFrames(node.Hm, link.Hm)
