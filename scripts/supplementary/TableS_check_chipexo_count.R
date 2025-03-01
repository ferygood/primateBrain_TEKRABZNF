# Table S

# We download the ChIP-exo data (GSE78099) to calculate if there are overlap
# between KAP1 with (1) Zinc finger proteins and (2) RepeatMasker Annotation

# the first part is using bedtools to do intersection, an example of the script:
#!/bin/bash
# for znf_file in *_ZNF*.bed; do
#     znf_name=$(echo "$znf_file" | awk -F'_' '{print $2}')
#     output_file="intersect_${znf_name}.bed"
#     bedtools intersect -a "$znf_file" -b GSM2067350_KAP1_exo_H1_peaks.bed ~/Downloads/hg38_rmsk.bed > "output/$output_file"
#     echo "Processed: $znf_file -> $output_file"
# done

# The second part is to extract the information and save it as counts.csv
#!/bin/bash
# echo "gene,count" > counts.csv
# for file in intersect_*.bed; do
#     gene=$(echo "$file" | sed 's/intersect_//; s/\.bed//')
#     count=$(wc -l < "$file")
#     echo "$gene, $count" >> counts.csv
# done

library(dplyr)
library(ggplot2)
library(ggpubr)

df_count <- read.csv("tables/counts.csv")

df_merge <- df_count %>%
    inner_join(kznf_infer[,c(2,6)], join_by(gene==external_gene_name)) %>%
    filter(count!=0) #ZNF182 has no overlap detected

g <- ggplot(df_merge, aes(x = age, y = log(count, base = 10))) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.5, color = "blue") +
    ylim(0, NA) +
    stat_compare_means(method = "wilcox.test", label = "p.signif") +
    labs(x = "Evolutionary age", y = "Count (log)", title = "Overlap binding events") +
    theme_minimal()

ggsave(filename="../../figures/S10_check_chipexo_overlap.jpg", g, dpi=400, width=4, height=4)

