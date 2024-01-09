library(TEKRABber)
library(twice)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(corrplot)
library(ggpmisc)

data("hmKZNFs337")

meta_combine <- metadata %>%
    inner_join(brain_meta, join_by(brain_region==region))

group_noHa <- meta_combine %>%
    filter(individual %in% c("hb", "hc", "hd")) %>%
    filter(cluster=="cluster1")
group_noHb <- meta_combine %>%
    filter(individual %in% c("ha", "hc", "hd")) %>%
    filter(cluster=="cluster1")
group_noHc <- meta_combine %>%
    filter(individual %in% c("ha", "hb", "hd")) %>%
    filter(cluster=="cluster1")
group_noHd <- meta_combine %>%
    filter(individual %in% c("ha", "hb", "hc")) %>%
    filter(cluster=="cluster1")

data("hg38_panTro6_rmsk")

check_human <- function(human_id){
    fetch <- orthologScale(
        speciesRef = "hsapiens",
        speciesCompare = "ptroglodytes",
        geneCountRef = hmGene[,c("geneID", human_id)],
        geneCountCompare = ptGene,
        teCountRef = hmTE[,c("name", human_id)],
        teCountCompare = ptTE[,-c(2,3)],
        rmsk = hg38_panTro6_rmsk
    )

    inputBundle <- DECorrInputs(fetch)

    meta <- data.frame(species=c(rep("human", ncol(fetch$geneRef) - 1),
                                 rep("chimpanzee", ncol(fetch$geneCompare) - 1))
    )

    rownames(meta) <- colnames(inputBundle$geneInputDESeq2)

    meta$species <- factor(meta$species, levels = c("human", "chimpanzee"))

    hmchimpDE <- DEgeneTE(
        geneTable = inputBundle$geneInputDESeq2,
        teTable = inputBundle$teInputDESeq2,
        expDesign = TRUE,
        metadata = meta
    )

    hmkznf_exp <- hmchimpDE$geneCorrInputRef %>%
        filter(rownames(.) %in% kznf_infer$ensembl_gene_id)

    hmCorr <- corrOrthologTE(
        geneInput = hmkznf_exp,
        teInput = hmchimpDE$teCorrInputRef,
        corrMethod = "pearson",
        padjMethod = "fdr"
    )

    hmCorr
}

hmCorr_noHa <- check_human(group_noHa$Run)
hmCorr_noHb <- check_human(group_noHb$Run)
hmCorr_noHc <- check_human(group_noHc$Run)
hmCOrr_noHd <- check_human(group_noHd$Run)

# all human
HmAll <- HmPtC1$corrRef %>% filter(pvalue<0.05 & coef<0)
PtC1 <- HmPtC1$corrCompare
PpC1 <- HmPpC1$corrCompare
MmC1 <- HmMmC1$HmMmC1_MmCorr

HmAll <- HmAll %>%
    left_join(kznf_infer[,c(1,2)],
              join_by(geneName==external_gene_name)) %>%
    select(c(6,2,3,4,5))
colnames(HmAll)[1] <- "geneName"

get_pair_list <- function(df){

    df_sig <- df %>%
        filter(coef<0 & pvalue<0.05) %>%
        mutate(pair = paste0(geneName, "-", teName))
    mylist <- df_sig$pair
    mylist
}

allhuman <- get_pair_list(HmAll)
noHa <- get_pair_list(hmCorr_noHa)
noHb <- get_pair_list(hmCorr_noHb)
noHc <- get_pair_list(hmCorr_noHc)
noHd <- get_pair_list(hmCOrr_noHd)

# Define the five character lists

# Create a list of the character lists
lists <- list(allhuman, noHa, noHb, noHc, noHd)

# Function to calculate the intersection percentage between two character lists
intersection_percentage <- function(list1, list2) {
    intersect_count <- length(intersect(list1, list2))
    union_count <- length(union(list1, list2))
    percentage <- (intersect_count / union_count) * 100
    return(percentage)
}

# Calculate the pairwise intersection percentages
result <- matrix(nrow = length(lists), ncol = length(lists))
colnames(result) <- names(lists)
rownames(result) <- names(lists)

for (i in 1:length(lists)) {
    for (j in 1:length(lists)) {
        result[i, j] <- intersection_percentage(lists[[i]], lists[[j]])
    }
}

col_names <- c("all human", "noHa", "noHb", "noHc", "noHd")
rownames(result) <- col_names
colnames(result) <- col_names

result_upper <- result
result_upper[lower.tri(result_upper)] <- NA


jpeg("figures/hm_corrplot_deseq2.jpg", units="in", width=5, height=5, res=300)
corrplot::corrplot(
    result/100, addCoef.col = 'black', col=corrplot::COL2("PiYG"), cl.pos = 'n',
    type='upper'
)
dev.off()
