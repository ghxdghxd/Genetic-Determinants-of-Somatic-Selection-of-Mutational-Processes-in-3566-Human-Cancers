# Figure S1

```R
library(ggplot2)
library(RColorBrewer)
library(ggsci)
library(ggrepel)

table_pca = "Table_PCA.csv"

pca <- read.table(table_pca, header=T, sep="\t", stringsAsFactors = F)
pca$pop = pca$population
pca$pop[which(pca$source=='TCGA')]="TCGA"

p_pca = ggplot() + geom_point(data=pca[which(pca$source=="TCGA"), ],
                      aes(x = PC1, y = PC2, colour = pop), size = 2.5, alpha = 1) +
    geom_point(data = pca[which(pca$source!="TCGA"), ],
               aes(x = PC1, y = PC2, colour = pop), size = 3, alpha = 0.8) +
    labs(x="Principal component 1",y="Principal component 2", title = "Population Stratification") +
    scale_colour_manual(values = c("CHB" = '#E41A1C', "MXL" = "#377EB8", "CEU" = "#4DAF4A", "ASW" = "#984EA3","YRI"="#FF7F00", "TCGA" = "gray"),
                        breaks = c("CEU", "CHB", "ASW", "MXL", "YRI", "TCGA"),
                        labels = c("CEU", "CHB", "ASW", "MXL", "YRI", "TCGA")) +
    theme_classic() +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 15),
          legend.title = element_blank(),
          legend.position = c(0.8, 0.7),
          legend.text = element_text(size = 12),
          legend.background = element_blank(),
          panel.background = element_blank(),
          plot.background = element_blank(),
          plot.title = element_text(hjust=0.5, size = 15)) +
    coord_equal() +
    geom_rect(aes(xmin=-12, xmax=10, ymin=27, ymax=40), color="gray70", alpha=0.1) +
    # geom_text(aes(x=-5, y=35, label="CHB"), size=4) +
    geom_rect(aes(xmin=-12, xmax=5, ymin=2, ymax=27), color="gray70", alpha=0.1) +
    # geom_text(aes(x=-5, y=15, label="MXL"), size=4) +
    geom_rect(aes(xmin=-12, xmax=2, ymin=-10, ymax=2), color="gray70", alpha=0.1) +
    # geom_text(aes(x=-5, y=0, label="CEU"), size=4) +
    geom_rect(aes(xmin=15, xmax=38, ymin=-12, ymax=2), color="gray70", alpha=0.1) +
    # geom_text(aes(x=30, y=0, label="ASW"), size=4) +
    geom_rect(aes(xmin=38, xmax=45, ymin=-12, ymax=2), color="gray70", alpha=0.1)
    # geom_text(aes(x=41, y=0, label="YRI"), size=4)


dat = as.data.frame(table(pca[pca$population=="CEU" & pca$source=="TCGA", "cancer_type"])) %>% mutate(text_y = cumsum(Freq) - Freq/2)
dat$Var1 = factor(dat$Var1, levels = rev(dat$Var1))
dat$lab = paste(dat$Var1, dat$Freq, sep=",")
p_cancer = ggplot(dat, aes(x="", y=Freq, fill=Var1)) + geom_col() +
    geom_label_repel(aes(label = lab, x=1.5, y = text_y),color = "white", segment.size = 0.2, nudge_x = 0.1) +
    coord_polar(theta = 'y') + scale_fill_simpsons() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          legend.position="none",
          plot.background = element_blank(),
          panel.background = element_blank())
pdf(width=3, height=3)
p_cancer
dev.off()

pdf("Figure_S1.pdf", width = 8, height = 5)
ggarrange(p_pca, p_cancer, nrow=1, widths = c(1,1.2), labels = c("A","B"), newpage = F)
dev.off()
```
