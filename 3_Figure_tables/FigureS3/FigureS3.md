# Figure S3

```R
library(reshape2)
library(tidyverse)
library(ggpubr)
library(cowplot)

pmLogSigMat <- read.table("TCGA_mutational_propensities.txt", header=T,sep="\t",stringsAsFactors = F)
immu <- read.table("The_Immune_Landscape_of_Cancer_ImmuneSubtype.csv", header=T,sep="\t",stringsAsFactors=F)
rownames(immu) = immu$TCGA.Participant.Barcode

sample_int = intersect(immu$TCGA.Participant.Barcode, rownames(pmLogSigMat))

mat = cbind(pmLogSigMat[sample_int, ], immu[sample_int, ])
mat$Immune.Subtype = factor(mat$Immune.Subtype, levels = paste0("C", 1:6))

#### Immune.Subtype VS signature
library(ggsci)
plot_sig_ImmSub <- function(sig){
    res = melt(pairwise.wilcox.test(mat[, sig], mat[, "Immune.Subtype"])$p.value) %>% filter(value < 0.05)
    my_comparisons = lapply(1:nrow(res), function(x){return(as.character(unlist(res[x, 1:2])))})
    p = ggviolin(na.omit(mat[, c(sig, "Immune.Subtype")]), x = "Immune.Subtype", y = sig, color = "Immune.Subtype", add = "boxplot", 
        add.params = list(fill = "white"), xlab = "Immune Subtype", ylab = paste(sig, "relative activities")) +
        theme(legend.position = "none", axis.text = element_text(size = 12), axis.title = element_text(size = 14)) +
        scale_color_d3() +
        stat_compare_means(aes(label = ..p.signif.., size = 10), comparisons = my_comparisons, method = "wilcox.test", paried = T, p.adjust.method = "bonferroni",) +
        stat_compare_means(label.y = max(mat[, sig], na.rm = T) * 3)
    return(p)
}

library(ggsci)
immu_sig = function(mat, immu, ylab="Activity"){
    mat = cbind(mat[sample_int, ], immu[sample_int, ])
    mat = mat[mat$Immune.Subtype != "C5", ]
    mat$Immune.Subtype = factor(mat$Immune.Subtype, levels = paste0("C", c(1:4,6)))
    m = na.omit(melt(na.omit(mat[, c(paste0("MP", 1:6), "Immune.Subtype","TCGA.Study")], 
        id.vars=c("Immune.Subtype", "TCGA.Study"))))
    # m$variable = gsub("MMR","dMMR",m$variable)
    # m$variable = factor(m$variable, levels = c("NpCpG","APOBEC", "dMMR", "Tobacco", "EPP1", "EPP2"))
    p = ggplot(m) +
        # geom_boxplot(aes(x = variable, y = value, fill = variable)) +
        geom_violin(aes(x = variable, y = value, fill = variable), scale = "width", trim = FALSE, adjust = .5, na.rm = FALSE) +
        facet_grid(~Immune.Subtype, scales = "free") + 
        guides(fill = guide_legend(nrow = 1)) + 
        ylab(ylab) + scale_fill_d3() + 
        theme_bw() +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.y = element_text(size = 10),
              axis.line= element_blank(),
              axis.title.x = element_blank(),
              plot.title = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.position = "top",
              legend.direction = "horizontal",
              legend.title = element_blank())
    return(p)
}

p_list = lapply(colnames(pmLogSigMat), function(x){
    return(plot_sig_ImmSub(x))
})
pdf(file = "sig_vs_ImmSub.pdf", width = 12, height = 6)
plot_grid(plotlist= p_list, nrow = 2, ncol = 3, rel_heights = c(5, 5, 5), rel_widths = c(5, 5))
dev.off()
```
