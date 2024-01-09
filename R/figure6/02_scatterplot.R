library(twice)
library(ggplot2)

tcx_kznfs <- mayoTEKRABber$tcxDE$normalized_gene_counts %>%
    data.frame()
tcx_TE <- mayoTEKRABber$tcxDE$normalized_te_counts %>%
    data.frame()

# select znf736, AluYh9
df_znf736_AluYh9 <- data.frame(
    t(tcx_kznfs["znf736", ]),
    t(tcx_TE["AluYh9", ]),
    c(rep("control", 23), rep("AD", 24))
)

colnames(df_znf736_AluYh9) <- c("znf736", "AluYh9", "group")
df_znf736_AluYh9$group <- factor(df_znf736_AluYh9$group,
                                 levels=c("control", "AD"))

gznf736_AluYh9 <- ggplot(df_znf736_AluYh9,
                         aes(x=log(znf736), y=log(AluYh9), fill=group)) +
    geom_point(colour="black", size=3, shape=21) +
    scale_fill_manual(values = c("#CC99FF", "#99CC99")) +
    facet_wrap(~group) +
    geom_smooth(method="lm", se=FALSE, show.legend = F) +
    xlab("znf736 (Log exp.)") +
    ylab("AluYh9 (Log exp.)") +
    theme_bw()
