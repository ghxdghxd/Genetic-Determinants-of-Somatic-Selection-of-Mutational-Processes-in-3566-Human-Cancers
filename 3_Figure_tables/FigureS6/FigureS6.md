# Figure S6

```R
file_Candidate_gene = "../TableS2/Table_S2A.csv"
file_E_gene = "../TableS2/Table_S2B.csv"
file_I_gene = "../TableS2/Table_S2C.csv"

eff_RData = "../Figure5/DepMap_19q2_dependency_score.RData"

E_genes <- read.table(file_E_gene, header = T, sep = "\t", stringsAsFactors = F)
I_genes <- read.table(file_I_gene, header = T, sep = "\t", stringsAsFactors = F)

load(eff_RData)
eff_raw = melt(depmap_eff, id.vars="cancer")
index_0.2 = lapply(seq(-4, -0.2, by = 0.2), function(x){return(c(x, x+0.2))})
eff_raw_median = eff_raw %>% group_by(variable) %>% summarise(value = median(value))

library(ggpubr)
library(egg)
library(ggsci)
pdf("FigureS6.pdf", width = 8, height = 6)
# gene number = 485/350/1427
rbind(eff_raw_median %>% filter(variable %in% E_genes$gene) %>% mutate(group = "E-genes"),
    eff_raw_median %>% filter(variable %in% intersect(E_genes$gene, I_genes$gene)) %>% mutate(group = "Overlapped-genes"),
    eff_raw_median %>% filter(variable %in% I_genes$gene) %>% mutate(group = "I-genes")) %>%
    ggboxplot(fill = "group", outlier.shape = NA, add.params = list(size=0.8), add = c("jitter"),
          order = c("E-genes", "Overlapped-genes","I-genes"), x= "group", xlab = "", y ="value",
          ylab = "Median CERES", 
          width = 0.98, ggtheme = theme_presentation()) + scale_fill_d3() +
    stat_compare_means(comparisons = list(c("E-genes","Overlapped-genes"), c("Overlapped-genes","I-genes"),c("E-genes","I-genes"))) +
    theme(aspect.ratio = 1, legend.position = "none")
dev.off()
```
