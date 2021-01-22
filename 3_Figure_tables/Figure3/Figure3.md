# Figure 3

## Figure 3A

```R
library(tidyverse)
library(readxl)
mat_gs_split_FDR <- read.table(gzfile("mat_rank_gs_split_FDR.txt.gz"), header = T, sep = "\t", stringsAsFactors = F)

plot_feature_adjR2 = function(signature, feature_index){
    colors = pal_nejm("default")(7)
    names(colors) = c("snv","scna","methy","miss","trun","large", "fusion")
    features = c("Germline_missenses", "Germline_truncation", "Germline_structural", "Somatic_nsy-mutation", "SCNAs", "TSS-methylation", "Fusion")
    names(features) = c("miss", "trun", "large", "snv", "scna", "methy", "fusion")
    p = ggplot(mat_gs_split_FDR %>% filter(sig == signature) %>% mutate(feature=factor(feature, levels = rev(feature_index)))) +
        geom_histogram(aes(x = adjr2, y = ..count.., fill = feature), position = "identity", binwidth = 0.001) +
        labs(x="adjusted Rsquare", y = "", fill="Gene statuses", colour="") +
        facet_grid(sig~.) + theme_article() +
        scale_y_sqrt(expand = c(-0.01, 0), limits = c(0, 2300)) +
        scale_x_continuous(expand = c(-0.01, 0), limits = c(-0.004, 0.541)) +
        scale_fill_manual(breaks = feature_index, values = colors[feature_index], labels = features[feature_index]) +
        theme(legend.position = "right",
              legend.key.size = unit(0.5, 'cm'),
              legend.text = element_text(size=12),
              strip.text.y = element_text(angle = 0, size = 14),
              strip.background = element_blank(),
              axis.title = element_text(size=14),
              axis.text.x = element_text(size=12),
              panel.grid = element_blank(),
              panel.spacing = unit(0.1,"lines"),
              plot.margin = margin(0,0,0,0,"lines"))
    return(p)
}

t = theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())
t1 = theme(legend.position = c(0.8,1),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())
t2 = theme(legend.position = "none")

pdf("Figure3A.pdf", width = 10, height = 6)
plot_grid(plot_feature_adjR2("MP1", feature_index = c("fusion", "large", "trun", "miss", "methy","snv", "scna")) + t,
          plot_feature_adjR2("MP2", feature_index = c("fusion", "large", "trun", "miss", "snv", "scna", "methy")) + t,
          plot_feature_adjR2("MP3", feature_index = c("fusion", "large", "trun", "miss", "methy", "scna", "snv")) + t1,
          plot_feature_adjR2("MP4", feature_index = c("fusion", "large", "trun", "miss", "methy", "scna", "snv")) + t,
          plot_feature_adjR2("MP5", feature_index = c("fusion", "large", "trun", "methy","miss", "scna", "snv")) + t,
          plot_feature_adjR2("MP6", feature_index = c("fusion", "large", "trun", "miss", "methy", "scna", "snv")) + t2,
          ncol = 1, align = "hv", axis = 'rl')
dev.off()
```

## Figure 3B

```R
color_palette = pal_nejm("default")(8)
names(color_palette) = c("Somatic_nsy_mutation", "SCNAs", "TSS_methylation", "Germline_Missenses", "Germline_Truncation", "Germline_Structural", "Fusion")
plot_gene_status_overlap = function(){
    miss = unique(mat_gs_split_FDR$gene[mat_gs_split_FDR$feature == "miss"])
    trun = unique(mat_gs_split_FDR$gene[mat_gs_split_FDR$feature == "trun"])
    large = unique(mat_gs_split_FDR$gene[mat_gs_split_FDR$feature == "large"])
    snv = unique(mat_gs_split_FDR$gene[mat_gs_split_FDR$feature == "snv"])
    scna = unique(mat_gs_split_FDR$gene[mat_gs_split_FDR$feature == "scna"])
    methy = unique(mat_gs_split_FDR$gene[mat_gs_split_FDR$feature == "methy"])
    fusion = unique(mat_gs_split_FDR$gene[mat_gs_split_FDR$feature == "fusion"])
    a = list(Germline_Missenses = miss, Germline_Truncation = trun, Germline_Structural = large,
        Somatic_nsy_mutation = snv, SCNAs = scna, Fusion = fusion, TSS_methylation = methy)
    venn(a, opacity = 0.7, box = F, ilabels = TRUE, size = 25, cexil = 1.2, cexsn = 1.5, zcolor = color_palette[names(a)])
}

pdf("Figure3B.pdf", width = 5, height = 5)
plot_gene_status_overlap()
dev.off()
```

# Figure 3C and 3D

```R
library(ggsci)

plot_driverGene_num = function(mat){
    a = as.data.frame.matrix(table(mat[,c("gene","sig")]))
    a[a>0]=1
    a = as.data.frame(table(rowSums(a)))
    a$type = "Signature"
    b = as.data.frame.matrix(table(mat[,c("gene","cancer")]))
    b[b>0]=1
    bb = as.data.frame(table(rowSums(b[,-5])))
    bb$type = "Cancer"
    bb$Freq[bb$Var1==0] = as.data.frame(table(b$panCan))[2, 2]
    dat = rbind(data.frame(Var1=NA, Freq=2314,type="bg"), a, bb)
    dat$Var1 = gsub("0","pan-Cancer",dat$Var1)
    dat$Var1[dat$Var1!="pan-Cancer" & dat$type == "Cancer"] = "Cancer-specific"
    dat$Var1 = factor(dat$Var1, levels = c(1:4,"pan-Cancer","Cancer-specific"))
    d3 = pal_d3("category10")(10)
    p = ggplot(na.omit(dat), aes(x=Var1, y=Freq)) + geom_bar(stat = "identity", width = 0.5) +
        facet_wrap(type~., scale="free") + theme_bw() +
        labs(x = "", y = "Gene count") +
        theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 14),
            strip.text = element_text(size =12),
            legend.position = "none",
            aspect.ratio = 1)
    return(p)
}

p_driverGene_num = plot_driverGene_num(mat_sig)

pdf("Figure3C_D.pdf", width = 8, height = 4)
p_driverGene_num
dev.off()
```
