library(dplyr)
library(tidyr)
library(readr)

human <- read.table("~/Downloads/hg38_rmsk.tsv", sep="/t")
chimp <- read.table("~/Downloads/panTro6_rmsk.tsv", sep="/t")
bonobo <- read.table("~/Downloads/panPan3_rmsk.tsv", sep="/t")
macaque <- read.table("~/Downloads/rheMac10_rmsk.tsv", sep="/t")

human_copy <- human %>%
    filter(V12 %in% c("LINE", "SINE", "LTR", "DNA", "Retroposon")) %>%
    count(V11, name = "human")

chimp_copy <- chimp %>%
    filter(V12 %in% c("LINE", "SINE", "LTR", "DNA", "Retroposon")) %>%
    count(V11, name = "chimp")

bonobo_copy <- bonobo %>%
    filter(V12 %in% c("LINE", "SINE", "LTR", "DNA", "Retroposon")) %>%
    count(V11, name = "bonobo")

macaque_copy <- macaque %>%
    filter(V12 %in% c("LINE", "SINE", "LTR", "DNA", "Retroposon")) %>%
    count(V11, name = "macaque")


merged_df <- human_copy %>%
    full_join(chimp_copy, by = "V11") %>%
    full_join(bonobo_copy, by = "V11") %>%
    full_join(macaque_copy, by = "V11") %>%
    mutate(across(where(is.numeric), ~ replace_na(.x, 0))) %>%  # Replace NA with 0 only in numeric columns
    rename(TE_subfamily = V11)  # Rename V11 to TE_subfamily

# Save to CSV
write_csv(merged_df, "TE_counts.csv")
