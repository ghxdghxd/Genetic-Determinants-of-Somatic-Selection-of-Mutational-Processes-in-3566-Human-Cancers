# Figure S4

```R
IV_immu_ICI_hallmakr <- read.csv("79_immu_ICI_PUCH.tsv", header = T, sep="\t",stringsAsFactors = F, nrows = 10, skip = 9)
IV_immu_ICI_hallmakr_gene <- read.csv("79_immu_ICI_PUCH.tsv", header = T, sep="\t",skip=24, stringsAsFactors = F)
I_genes <- read.table("../TableS2/Table_S2C.csv", header=T,sep="\t",stringsAsFactors = F)

IV_immu_ICI_hallmakr_gene = IV_immu_ICI_hallmakr_gene[rowSums(IV_immu_ICI_hallmakr_gene[,-c(1:3)]!="") > 0, ]
IV_immu_ICI_hallmakr_gene[, -c(1:3)][IV_immu_ICI_hallmakr_gene[,-c(1:3)]!=""] = 1
a = merge(IV_immu_ICI_hallmakr_gene, I_genes %>% filter(cancer_type=="pan-cancer") %>% select(mutational_propensity, gene), by.x = "Gene.Symbol", by.y = "gene", all.x = T)
a = melt(a[,-c(2:3)], id.vars = c("Gene.Symbol", "mutational_propensity"))
a = a[which(a$value==1),]
a = merge(a, IV_immu_ICI_hallmakr, by.x = "variable", by.y = "Gene.Set.Name", all.x = T)
# a = unique(a[,-c(2,4)])
a = a %>% group_by(variable, mutational_propensity) %>% mutate(n = sum(as.numeric(value))) %>%select(FDR.q.value, n) %>% unique() %>% as.data.frame()

pdf("FigureS7.pdf", width = 8, height = 4)
ggplot(a %>% filter(FDR.q.value < 0.1), aes(variable, mutational_propensity, color = -log10(FDR.q.value), size = n)) +
    geom_point() + scale_colour_steps2() +
    labs(size = "Gene counts", color = "-log10(FDR)") +
    facet_grid(mutational_propensity ~ variable, scales = "free", space = "free") +
    scale_size(range = c(3,6)) +
    theme(axis.text.x = element_text(angle = -30, hjust = 0, vjust = 1, size = 12),
          axis.text.y = element_text(size = 10),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          strip.text = element_blank(),
          strip.background = element_blank(),
          legend.box = "horizontal",
          panel.spacing = unit(0.1, 'lines'),
          panel.grid.major = element_blank(),
          aspect.ratio = 1)
dev.off()

```
