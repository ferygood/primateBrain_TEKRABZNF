library(tidyverse)

#1 read all rds file and save as a global variable

read_rds_files <- function(folder_path) {
    # Get a list of all .rds files in the folder
    file_list <- list.files(path = folder_path, pattern = "\\.rds$", full.names = TRUE)

    # Read each .rds file and save it as a global variable with the file prefix as the variable name
    for (file_path in file_list) {
        file_name <- basename(file_path)
        var_name <- gsub("\\.rds$", "", file_name)
        assign(var_name, readRDS(file_path), envir = .GlobalEnv)
    }
}

#read_rds_files("~/path/to/results_rdata")

#2 distributed negative correlation collection
# human
hm_df_list <- list(HmPtC1$HmPtC1_HmCorr, HmPtC2$HmPtC2_HmCorr, HmPtC3$HmPtC3_HmCorr, HmPtC4$HmPtC4_HmCorr, HmPtC5$HmPtC5_HmCorr, HmPtC6$HmPtC6_HmCorr, HmPtC7$HmPtC7_HmCorr)

# chimpanzee
chimp_df_lsit <- list(
    HmPtC1$HmPtC1_PtCorr, HmPtC2$HmPtC2_PtCorr, HmPtC3$HmPtC3_PtCorr,
    HmPtC4$HmPtC4_PtCorr, HmPtC5$HmPtC5_PtCorr, HmPtC6$HmPtC6_PtCorr,
    HmPtC7$HmPtC7_PtCorr
)

# bonobo
bonobo_df_lsit <- list(
    HmPpC1$HmPpC1_PpCorr, HmPpC2$HmPpC2_PpCorr, HmPpC3$HmPpC3_PpCorr,
    HmPpC4$HmPpC4_PpCorr, HmPpC5$HmPpC5_PpCorr, HmPpC6$HmPpC6_PpCorr,
    HmPpC7$HmPpC7_PpCorr
)

# macaque
macaque_df_lsit <- list(
    HmMmC1$HmMmC1_MmCorr, HmMmC2$HmMmC2_MmCorr, HmMmC3$HmMmC3_MmCorr,
    HmMmC4$HmMmC4_MmCorr, HmMmC5$HmMmC5_MmCorr, HmMmC6$HmMmC6_MmCorr,
    HmMmC7$HmMmC7_MmCorr
)


get_neg_count <- function(df_list, init_var){
    for (df in df_list){
        df_sig <- df %>%
            dplyr::filter(pvalue < 0.05 & coef < 0)
        scaling <- length(unique(df_sig$geneName)) * length(unique(df_sig$teName))
        #print(scaling)
        # scale
        row_count <- df_sig %>% nrow()
        row_count <- (row_count * 1000) / scaling
        init_var <- c(init_var, row_count)
    }
    round(init_var,2)
}



human <- c()
human_count <- get_neg_count(hm_df_list, human)

chimpanzee <- c()
chimp_count <- get_neg_count(chimp_df_lsit, chimpanzee)

bonobo <- c()
bonobo_count <- get_neg_count(bonobo_df_lsit, bonobo)

maca <- c()
maca_count <- get_neg_count(macaque_df_lsit, maca)


df <- data.frame(
    human = human_count,
    chimpanzee = chimp_count,
    bonobo = bonobo_count,
    macaque = maca_count)
rownames(df) <- c("primary & secondary cortices", "limbic & association cortices", "archicortex", "hypothalamus & thalamus", "cerebellar white matter", "cerebellar grey matter", "striatum")

df$brain <- rownames(df)

df_long <- df %>%
    pivot_longer(cols=c("human", "chimpanzee", "bonobo", "macaque"), names_to = "species", values_to = "counts")
df_long$brain <- factor(df_long$brain, levels=c(
    "primary & secondary cortices", "limbic & association cortices", "archicortex", "hypothalamus & thalamus", "cerebellar white matter", "cerebellar grey matter", "striatum"
))
df_long$species <- factor(df_long$species, levels=c("human", "chimpanzee", "bonobo","macaque"))

g <- ggplot(df_long, aes(x=forcats::fct_rev(brain), y=counts, fill=species)) +
    geom_bar(stat='identity') +
    scale_fill_manual(values=c("red", "blue", "purple", "darkgreen")) +
    facet_wrap(~species, nrow = 1) + coord_flip() +
    ylab("Normalized negatively correlated TE:KRAB-ZNF counts")+
    xlab("") +
    theme_bw()
