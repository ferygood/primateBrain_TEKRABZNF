library(ggplot2)
hg38TE <- read.csv("~/github/pBrain/data/hg38rmsk_info.csv")
young_TE <- read.csv("~/github/pBrain/data/Dfam_TE_simiiformes.csv")

hg38TE_age <- hg38TE %>%
    mutate(age = ifelse(repName %in% young_TE$NM, "young", "old")) %>%
    select(c(3,4))

hg38_old <- hg38TE_age %>%
    filter(age=="old") %>%
    count(repFamily) %>%
    mutate(group=">44.2 mya")

hg38_young <- hg38TE_age %>%
    filter(age=="young") %>%
    count(repFamily) %>%
    mutate(group="after Simiiformes \n(<44.2mya)")

dfTE <- rbind(hg38_old, hg38_young)
svas <- c("SVAs", 4, "after Simiiformes \n(<44.2mya)")
dfTE <- rbind(dfTE, svas)
dfTE$n <- as.numeric(dfTE$n)

gTE <- ggplot(dfTE, aes(x=repFamily, y=n, fill=group)) +
    geom_bar(stat="identity", position="dodge") +
    theme_bw() +
    scale_fill_manual(values = c("#5A8FBB", "#E59E00")) +
    coord_flip() +
    ylab("number of TEs") +
    xlab("TE family") +
    theme(legend.position = "bottom", legend.text = element_text(size=8))

ggsave(filename="figures/age_TE.jpg", gTE, dpi=300, height=6, width=4)
