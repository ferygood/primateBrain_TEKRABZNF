library(dplyr)
library(twice)

# loading tables
data("primateGRFs")
kznf_infer_new <- kznf_infer
gentree <- readxl::read_xlsx("~/Downloads/hg38_ver95_age(2).xlsx")

kznf_infer_new <- kznf_infer_new %>%
    left_join(gentree[,c(1,3)], join_by(ensembl_gene_id==gene))

# 127 needs to be check
kznf_infer_check <- subset(kznf_infer_new, branch.y >= 8 | is.na(branch.y))

# select ensembl id and gene name
kznf_infer_select <- kznf_infer_check[,c(1,2)]


# read Vladi's three bucket file
table_v <- read.table("~/Downloads/yaos_KRAB_ZNFs_ages.txt", header = TRUE)

# combine
kznf_infer_select <- kznf_infer_select %>%
    left_join(table_v, join_by(ensembl_gene_id==ensembl_gene_id))

# fill in the forms for empty data: znf747, znf91, znf723
kznf_infer_select[8,c(3,4,5)] <- c(3, "no", "no") #znf747
kznf_infer_select[51,c(3,4,5)] <- c(1, "no", "no") #znf91
kznf_infer_select[98,c(3,4,5)] <- c(2, "yes", "no") #znf723

# combine with all old data
kznf_all <- kznf_infer[,c(1,2)] %>%
    left_join(
        kznf_infer_select[,c(1,3,4,5)], join_by(ensembl_gene_id==ensembl_gene_id))

#fill na with 1 bucket
kznf_all$bucket <- ifelse(is.na(kznf_all$bucket), 1, kznf_all$bucket)
kznf_all$great_apes <- ifelse(is.na(kznf_all$great_apes), "no", kznf_all$great_apes)
kznf_all$chimp_human <- ifelse(
    is.na(kznf_all$chimp_human), "no", kznf_all$chimp_human)

kznf_all$age <- ifelse(kznf_all$bucket %in% c(1,3), "old", "young")

write.csv(kznf_all, file="../../posts/manuscript/tables/kznf_bucket.csv", row.names=FALSE)
