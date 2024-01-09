library(RCy3)
library(twice)
library(condor)
library(ggpubr)

data("hg19rmsk_info")
hs_tcx <- read.csv("tables/tcx_lostAD_link.csv")

#1.  We first use condor to calculate the modularity:

df <- hs_tcx[,c(2,1,4)] # select TE, KRAB-ZNFs, coefficient
df$coef <- abs(df$coef)
condor.object <- create.condor.object(df)
condor.object <- condor.cluster(condor.object, project = F, cs.method = "LEC")

#2.  Then we explore the distribution of nodes from each clusters (we have 10)

#TE nodes
te_com <- condor.object$red.memb
te_com_df <- te_com %>%
    left_join(hg19rmsk_info[,c(1,2)], join_by(red.names==gene_id))
te_com_df$com <- as.character(te_com_df$com)

df <- te_com_df %>%
    group_by(com, family_id) %>%
    summarise(Count = n())

df$com <- factor(df$com, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))

other_list <- c("hAT-Tip100", "CR1", "DNA", "TcMar", "TcMar-Mariner", "TcMar-Tc2",
                "RTE","Gypsy", "Dong-R4", "LTR", "hAT", "hAT-Blackjack", "PiggyBac",
                "MIR", "MuDR")
te_com_merge <- df %>%
    mutate(family_id = ifelse(family_id %in% other_list, "others", family_id)) %>%
    mutate(family_id =
               ifelse(family_id %in% c("SVA_A", "SVA_C", "SVA_D", "SVA_E", "SVA_F"), "SVAs", family_id)) %>%
    mutate(family_id=
               ifelse(family_id %in% c("ERV1", "ERVK", "ERVL", "ERVL-MaLR"), "ERVs", family_id))

te_com_merge$com <- as.character(te_com_merge$com)
te_com_merge$com <- factor(te_com_merge$com, levels=c(1:10))
te_com_merge$family_id <- factor(
    te_com_merge$family_id,
    levels=c("Alu", "ERVs", "hAT-Charlie", "L1", "SVAs", "TcMar-Tigger", "others"))

# Define color palette
custom_colors <- c("#ffa6c9", "#ffb81c", "#003399", "#a31f34", "#00a4b4", "#2d4a64", "#942ccc")

# Use ggplot to create the visualization
condor_TE <- ggplot(te_com_merge, aes(fill = family_id, y = Count, x = com)) +
    geom_bar(position = "stack", stat = "identity") +
    scale_fill_manual(values = custom_colors) +
    xlab("Clusters") +
    theme_bw()

# KRAB-ZNF nodes
znf_com <- condor.object$blue.memb
colnames(znf_com) <- c("id", "com")
znf_com_process <- znf_com %>%
    left_join(kznf_infer[c(2,6)], join_by(id==external_gene_name))
znf_com_process$com <- factor(znf_com_process$com, levels=c(seq(1:10)))
znf_com_process_count <- znf_com_process %>%
    group_by(com, age) %>%
    summarise(count=n())

condor_kznf <- ggplot(znf_com_process_count, aes(fill=age, y=count, x=com)) +
    geom_bar(position="stack", stat="identity") +
    scale_fill_manual(values = c("#5A8FBB", "#E59E00")) +
    xlab("Clusters") +
    ylab("Count") +
    theme_bw()

# combine figure
g_znf_te <- ggarrange(condor_kznf, condor_TE, nrow=1)
ggsave(g_znf_te, file="figures/condor_node_count.jpg", dpi=300, width=8, height=3)
