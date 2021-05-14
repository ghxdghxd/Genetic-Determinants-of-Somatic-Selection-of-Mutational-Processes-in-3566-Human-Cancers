# Figure 3

## Figure 3A

```R
library(tidyverse)
library(readxl)
library(ggrepel)
mat_gs_split_FDR <- read.table(gzfile("mat_rank_gs_split_FDR.txt.gz"), header = T, sep = "\t", stringsAsFactors = F)

plot_feature_adjR2 = function(signature, feature_index, lab_genes = NULL, lab_features = NULL){
    colors = pal_nejm("default")(7)
    names(colors) = c("snv","scna","methy","miss","trun","large", "fusion")
    features = c("Germline_missenses", "Germline_truncation", "Germline_structural", "Somatic_nsy-mutation", "SCNAs", "TSS-methylation", "Fusion")
    names(features) = c("miss", "trun", "large", "snv", "scna", "methy", "fusion")
    p = ggplot(mat_gs_split_FDR %>% filter(sig == signature) %>% mutate(feature=factor(feature, levels = rev(feature_index)))) +
        geom_histogram(aes(x = adjr2, y = ..count.., fill = feature), position = "identity", binwidth = 0.001)
    if(!is.null(lab_genes) & is.null(lab_features)){
        m = mat_gs_split_FDR %>% filter(sig == signature, gene %in% lab_genes) %>%
            mutate(feature=factor(feature, levels = rev(feature_index)),
                lab = ifelse(gene %in% lab_genes, paste(gene, cancer, sep = ","), ""))
        p = p + geom_text_repel(data = m, aes(x = adjr2, y = 0, color = feature, label = lab), segment.color = 'grey50', nudge_y = 1000)
    }else if(!is.null(lab_genes) & !is.null(lab_features)){
        names(lab_features) = lab_genes
        m = mat_gs_split_FDR %>% filter(sig == signature, gene %in% lab_genes, feature %in% lab_features) %>%
            mutate(feature=factor(feature, levels = rev(feature_index)),
                lab = ifelse(gene %in% lab_genes & feature %in% lab_features[gene], paste(gene, cancer, sep = ","), ""))
        p = p + geom_text_repel(data = m, aes(x = adjr2, y = 0, color = feature, label = lab), segment.color = 'grey50', nudge_y = 1000)
    }
    p = p + labs(x="adjusted Rsquare", y = "", fill="Gene statuses", colour="") +
        facet_grid(sig~.) + theme_article() +
        scale_y_sqrt(expand = c(-0.01, 0), limits = c(0, 2300)) +
        scale_x_continuous(expand = c(-0.01, 0), limits = c(-0.004, 0.541)) +
        scale_fill_manual(breaks = feature_index, values = colors[feature_index], labels = features[feature_index]) +
        scale_color_manual(breaks = feature_index, values = colors[feature_index], labels = features[feature_index]) +
        theme(legend.position = "right",
              legend.key.size = unit(0.5, 'cm'),
              legend.text = element_text(size=12),
              strip.text.y.left = element_text(angle = 0, hjust = 1, vjust = 0.5),
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
    axis.title.x = element_blank())
tt = t + theme(axis.ticks.x = element_blank())


library(ggrepel)
pdf("Figure3A.pdf", width = 10, height = 6)
plot_grid(plot_feature_adjR2("MP1", feature_index = c("fusion", "large", "trun", "miss", "methy","snv", "scna"),
            lab_genes = c("MBD4", "BRCA2"), lab_features = c("trun", "trun")) + tt,
          plot_feature_adjR2("MP2", feature_index = c("fusion", "large", "trun", "miss", "snv", "scna", "methy"),
            lab_genes = c("APOBEC3H","APOBEC3H","APOBEC1","APOBEC2","APOBEC3B"), lab_features = c("miss", "methy", "methy", "methy", "methy")) + tt,
          plot_feature_adjR2("MP3", feature_index = c("fusion", "large", "trun", "miss", "methy", "scna", "snv"),
            lab_genes = c("PMS2", "MSH6","PMS2","MLH1","MSH2","MSH6"), lab_features = c("miss", "miss","snv","snv","snv","snv")) + tt,
          plot_feature_adjR2("MP4", feature_index = c("fusion", "large", "trun", "miss", "methy", "scna", "snv")) + tt,
          plot_feature_adjR2("MP5", feature_index = c("fusion", "large", "trun", "methy","miss", "scna", "snv"),
          lab_genes = c("POLE", "POLE2","POLE3", "POLD1"), lab_features = c("snv", "snv","snv","snv")) + tt,
          plot_feature_adjR2("MP6", feature_index = c("fusion", "large", "trun", "miss", "methy", "scna", "snv"),
          lab_genes = c("POLE", "POLE2","POLE3", "POLD1"), lab_features = c("snv", "snv","snv","snv")) + t,
          ncol = 1, align = "hv", axis = 'rl')
dev.off()
```

## Figure 3B

```R
library(UpSetR)
upset_gene = as.data.frame.matrix(table(mat_gs_split_FDR[,c("gene","feature")]))
upset_gene[upset_gene > 0] = 1
upset_gene$gene = "B"
upset_gene$gene[which(rowSums(upset_gene[,c("fusion","snv","scna","methy")]) > 0 & rowSums(upset_gene[,c("miss","trun","large")]) > 0)]="A"

variant_list = c("Somatic_nsy_mutation", "SCNAs", "TSS_methylation", "Germline_Missenses", "Germline_Truncation", "Germline_Structural", "Fusion")
names(variant_list) = c("snv","scna", "methy", "miss", "trun", "large","fusion")
colnames(upset_gene)[1:7] = variant_list[colnames(upset_gene)[1:7]]

among <- function(row, list1, list2){
  newData <- (sum(as.numeric(row[list1])) > 0) & (sum(as.numeric(row[list2])) > 0)
}

pdf("Figure3B.pdf", width = 7, height = 5)
upset(upset_gene, sets = rev(c("Germline_Missenses", "Germline_Truncation", "Germline_Structural", "Somatic_nsy_mutation", "SCNAs", "Fusion", "TSS_methylation")),
      line.size = 0.3, point.size = 2, set_size.show = T, keep.order = T, shade.color = "grey50",
      queries = list(list(query = elements, params = list("gene","B"), color = "gray40", active = T)))
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
