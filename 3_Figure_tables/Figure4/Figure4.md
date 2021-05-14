# Figure 4

## Figure 4A

```R
file_Candidate_gene = "../TableS2/Table_S2A.csv"
file_E_gene = "../TableS2/Table_S2B.csv"
file_I_gene = "../TableS2/Table_S2C.csv"
file_cosmic_gene = "Cancer_Gene_Census_allThu_Mar_14_06_02_23_2019.csv"
file_susceptibility_gene = "Cancer-susceptibility-genes.txt"
file_Compendium_Cancer_Gene = "Compendium_Cancer_Genes.tsv"
file_Nucleotide_context_gene = "Nucleotide_context_gene.csv"


Candidate_genes <- read.table(file_Candidate_gene, header = T, sep = "\t", stringsAsFactors = F)
E_genes <- read.table(file_E_gene, header = T, sep = "\t", stringsAsFactors = F)
I_genes <- read.table(file_I_gene, header = T, sep = "\t", stringsAsFactors = F)

# Cancer Gene Census(CGCs) in Catalogue of Somatic Mutations in Cancer (COSMIC)
CGCs = read.table(file_cosmic_gene, header=T,sep=",",stringsAsFactors=F)

# Pathogenic Germline Variants in 10,389 Adult Cancers
Susceptibility_genes = read.csv(file_susceptibility_gene, header = T, sep="\t",stringsAsFactors = F)

# A compendium of mutational cancer driver genes
Compendium_Cancer_Genes = read.csv(file_Compendium_Cancer_Gene, sep = "\t",header = T, stringsAsFactors=F)

# Identification of cancer driver genes based on nucleotide context
Nucleotide_context_genes = read.csv(file_Nucleotide_context_gene, sep = "\t",header = T, stringsAsFactors=F)


upset_gene = data.frame(gene = unique(c(Candidate_genes$gene, E_genes$gene, I_genes$gene, CGCs$Gene.Symbol,
    Compendium_Cancer_Genes$SYMBOL, Nucleotide_context_genes$Gene, Susceptibility_genes$Gene..Symbol)))
upset_gene$Candidate_genes = ifelse(upset_gene$gene %in% Candidate_genes$gene, 1, 0)
upset_gene$Egenes = ifelse(upset_gene$gene %in% E_genes$gene, 1, 0)
upset_gene$Igenes = ifelse(upset_gene$gene %in% I_genes$gene, 1, 0)
upset_gene$COSMIC = ifelse(upset_gene$gene %in% CGCs$Gene.Symbol, 1, 0)
upset_gene$Compendium_Cancer_Genes = ifelse(upset_gene$gene %in% Compendium_Cancer_Genes$SYMBOL, 1, 0)
upset_gene$Nucleotide_context = ifelse(upset_gene$gene %in% Nucleotide_context_genes$Gene, 1, 0)
upset_gene$Susceptibility_genes = ifelse(upset_gene$gene %in% Susceptibility_genes$Gene..Symbol, 1, 0)

pdf("Figure4A.pdf", width = 7, height = 5)
upset(upset_gene, sets = rev(c("Candidate_genes","Egenes", "Igenes", "COSMIC", "Compendium_Cancer_Genes", "Nucleotide_context", "Susceptibility_genes")),
      sets.x.label = "Gene set size", keep.order = T, line.size = 0.3, point.size = 2, set_size.show = T,shade.color = "grey50",)
dev.off()
```

## Figure 4B

```R
# E-gene与I-gene overlap
C = cbind(unique(Candidate_genes[Candidate_genes$cancer_type=="pan-cancer", c("gene","mutational_propensity")]), type = "Candidate")
E = cbind(unique(E_genes[E_genes$cancer_type=="pan-cancer", c("gene","mutational_propensity")]), type = "E-gene")
I = cbind(unique(I_genes[I_genes$cancer_type=="pan-cancer", c("gene","mutational_propensity")]), type = "I-gene")
# D = cbind(merge(E,I, by=c("gene","mutational_propensity"))[,1:2], type = "Dual-gene")
# D = merge(E_genes[,c("gene","sig","cancer","feature")], I_genes[,c("gene","sig","cancer","feature")], by = "gene") %>% unique()
# write.table(D,"E_I_overlap_350.txt", c = T, r = F, sep = "\t", quote = F)
# # E = merge(E, D, by=c("gene","sig"), all.x = T) %>% filter(is.na(type.y)) %>% select(gene, sig) %>% mutate(type = "E-gene")
# # I = merge(I, D, by=c("gene","sig"), all.x = T) %>% filter(is.na(type.y)) %>% select(gene, sig) %>% mutate(type = "I-gene")
m = rbind(E, I)

library(GeneOverlap)
res = as.data.frame(do.call(rbind, lapply(c("I", "E"), function(x){
    res = t(sapply(sort(unique(C$mutational_propensity)), function(y){
        res = testGeneOverlap(newGeneOverlap(get(x)$gene, C$gene[C$mutational_propensity == y], genome.size = 2314))
        return(c(p = res@pval, OR = res@odds.ratio, mutational_propensity = y))
    }))
    return(cbind(res, type=x))
})))

res$type[res$type=="I"] = "I-gene"
res$type[res$type=="E"] = "E-gene"

a = merge(melt(table(m[,2:3])), res, by=c('mutational_propensity',"type")) %>% group_by(type) %>%
    mutate(rate = value/sum(value), p = as.numeric(p), OR = as.numeric(OR)) %>%
    mutate(type = factor(type, levels = c("E-gene","I-gene")), P = -log10(p)),
           mutational_propensity = factor(mutational_propensity, levels = c("MP1", "MP2", "MP3", "MP4", "MP5", "MP6")))
a$P[a$P > 10] = 10
a$OR[a$OR > 10] = 10

p_E_I_overlap = ggplot(a, aes(mutational_propensity, type, size = OR, color = P)) + geom_point() +
    facet_grid(type ~ mutational_propensity, scales = "free", space = "free") +
    labs(x="",y="", size = "Ratio of enrichment", color = "-log10(P)") +
    scale_color_steps(low = "#BDD7E7", high = "#08519C", n.breaks = 5) +
    scale_size(range = c(3, 8)) +
    theme(axis.text.x = element_text(hjust=0, vjust=1, size = 12, angle = -60),
          axis.text.y = element_text(size = 10),
          #axis.text = element_blank(),
          axis.title = element_text(size = 14),
          axis.ticks = element_blank(),
          legend.text = element_text(size = 12),
          legend.position = "bottom",
          legend.direction = "vertical",
          legend.box = "horizontal",
          legend.background = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          panel.grid = element_blank(),
          plot.background = element_blank(),
          panel.spacing = unit(0.1, 'lines'),
          aspect.ratio = 1)

pdf("Figure4B.pdf", width = 3, height = 5)
p_E_I_overlap
dev.off()
```

# Figure 4C

```R
IV_reactome <- read.csv("E_genes_reactome.tsv", header = T, sep="\t",stringsAsFactors = F)
IV_reactome$Description = factor(IV_reactome$Description, levels = rev(IV_reactome$Description[order(IV_reactome$k.K)]))
IV_reactome$Gene.Set.Name = factor(IV_reactome$Gene.Set.Name, levels = IV_reactome$Gene.Set.Name[order(IV_reactome$k.K)])

IV_immu_reactome <- read.csv("I_genes_reactome.tsv", header = T, sep="\t",stringsAsFactors = F)
IV_immu_reactome$Description = factor(IV_immu_reactome$Description, levels = rev(IV_immu_reactome$Description[order(IV_immu_reactome$k.K)]))
IV_immu_reactome$Gene.Set.Name = factor(IV_immu_reactome$Gene.Set.Name, levels = IV_immu_reactome$Gene.Set.Name[order(IV_immu_reactome$k.K)])

IV_immu_GO <- read.csv("I_genes_GOBP.tsv", header = T, sep="\t",stringsAsFactors = F)
IV_immu_GO = IV_immu_GO[IV_immu_GO$Gene.Set.Name == "GO_BIOLOGICAL_ADHESION", ]
IV_immu_GO$Description = "Go biological adhesion"
IV_immu_GO$FDR.q.value = 1e-8

IV_immu_hallmakr <- read.csv("I_genes_hallmark.tsv", header = T, sep="\t",stringsAsFactors = F)
IV_immu_hallmakr = IV_immu_hallmakr[IV_immu_hallmakr$Gene.Set.Name == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", ]
IV_immu_hallmakr$Description = "Hallmark_epithelial_mesenchymal_transition"

IV_immu_hallmakr$genes = "I-genes"
IV_immu_GO$genes = "I-genes"
IV_reactome$genes = "E-genes"
IV_immu_reactome$genes = "I-genes"
reactome = rbind(IV_reactome, IV_immu_reactome, IV_immu_GO,IV_immu_hallmakr)
reactome$genes = factor(reactome$genes, levels = c("E-genes", "I-genes"))
reactome$Description = factor(reactome$Description, levels = unique(reactome$Description))

p_reactom = ggplot(reactome, aes(Description, genes, size = k.K, color = -log10(FDR.q.value))) + geom_point() +
    facet_grid(genes~Description, space = "free", scale="free") +
    labs(y = "", x="", color = "-log10(FDR)", size = "k/K") +
    scale_size(range = c(3, 8)) + scale_colour_steps2() +
    theme(axis.text.x = element_text(hjust=1, vjust=1, size = 12, angle = 60),
        axis.text.y = element_text(size = 10),
        #axis.text = element_blank(),
        axis.title = element_text(size = 14),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 12),
        legend.position = "right",
        legend.direction = "vertical",
        legend.box = "vertical",
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_blank(),
        panel.spacing = unit(0.1, 'lines'),
        aspect.ratio = 1)

pdf("Figure4C.pdf", width = 7, height = 8)
p_reactom
dev.off()
```

## Figure 4D

```R
library(networkD3)

fix_name = c("M1 macrophages", "T follicular helper cells", "CD8 T cells", "Neutrophils", "Resting mast cells",
    "Memory B cells", "Mast cells", "Eosinophils", "Activated mast cells", "M0 macrophages", "M2 macrophages",
    "Activated dendritic cells", "Plasma cells", "Naive B cells", "Naive CD4 T Cells", "Activated NK cells",
    "Gammadelta T cells", "Monocytes", "Activated memory CD4 T cells", "Resting NK cells", "Resting memory CD4 T cells", 
    "Regulatory T cells", "SNV Neoantigens", "Indel Neoantigens", "Unknown")
names(fix_name) = c("Macrophages_M1", "T_Cells_Follicular_Helper","T_Cells_CD8","Neutrophils","Mast_Cells_Resting",
    "B_Cells_Memory","Mast_Cells","Eosinophils","Mast_Cells_Activated","Macrophages_M0","Macrophages_M2",
    "Dendritic_Cells_Activated","Plasma_Cells","B_Cells_Naive","T_Cells_CD4_Naive","NK_Cells_Activated",
    "T_Cells_gamma_delta","Monocytes","T_Cells_CD4_Memory_Activated","NK_Cells_Resting","T_Cells_CD4_Memory_Resting",
    "T_Cells_Regulatory_Tregs", "SNV_Neoantigens_rate_log", "Indel_Neoantigens_rate_log", "Unknown")

plot_sankey_networkD3 = function(m, add_cancer = F, add_neoantigen = F){
    select_name = c("variant", "Immune_cell_type", "mutational_propensity", "group")
    m$Immune_cell_type = fix_name[m$Immune_cell_type]
    m$group = m$mutational_propensity
    m$Immune_cell_type[-which(m$Immune_cell_type %in% c("T follicular helper cells","M1 macrophages", "CD8 T cells", "Neutrophils", "Resting mast cells", "Memory B cells"))] = "Others"
    m$variant[which(m$variant %in% c("Fusion","Germline_missense", "Germline_truncation", "Germline_structural_varinat"))] = "Others"
    
    m = m %>% ungroup %>% group_by(cancer_type) %>% mutate(cancer_type=paste0(cancer_type, "(", length(unique(gene)), ")")) %>%
        ungroup %>% group_by(variant) %>% mutate(variant=paste0(variant,"(", length(unique(gene)), ")")) %>%
        ungroup %>% group_by(Immune_cell_type) %>% mutate(Immune_cell_type=paste0(Immune_cell_type,"(", length(unique(gene)), ")")) %>%
        ungroup %>% group_by(mutational_propensity) %>% mutate(mutational_propensity=paste0(mutational_propensity, "(", length(unique(gene)), ")")) %>%
        ungroup %>% select(all_of(select_name))
        # select(feature, cancer, neoantigens_rate_log, cell, sig, group)
    nodes = unique(gather(m[,1:(ncol(m)-1)], type, name))
    nodes = nodes[nodes$type!='gene',]
    nodes = cbind(nodes, as.data.frame(str_split(nodes$name, pattern = "[()]",simplify = T))[,1:2])
    colnames(nodes) = c("type","name", "ID", "Num")
    nodes$group = nodes$type
    nodes$group[which(nodes$type=="mutational_propensity")] = nodes$ID[which(nodes$type=="mutational_propensity")]
    a = 1:(ncol(m)-2)
    b = 2:(ncol(m)-1)
    edges = m[,c(a[1]:b[1], ncol(m))]
    colnames(edges) = c(colnames(m)[a[2]], colnames(m)[b[2]], "group")
    for(i in 2:length(a)){
        edges = rbind(edges, m[,c(a[i]:b[i], ncol(m))])
        if(i < (ncol(m)-2)){
            colnames(edges) = c(colnames(m)[a[i+1]], colnames(m)[b[i+1]], "group")
        }else{
            colnames(edges) = c("from", "to", "group")
        }
    }
    edges = edges %>% group_by(from, to, group) %>% summarise(value = length(to)) %>% as.data.frame
    edges$from = match(edges$from, nodes$name) - 1
    edges$to = match(edges$to, nodes$name) - 1
    my_color <- 'd3.scaleOrdinal() .domain(["MP1", "MP2", "MP3", "MP4", "MP5", "MP6", "Immune_cell_type", "variant"]) .range(["#1F77B4FF", "#FF7F0EFF", "#2CA02CFF", "#D62728FF", "#9467BDFF", "#8C564BFF", "#7F7F7FFF", "#BCBD22FF"])'
    sankeyNetwork(Links = edges, Nodes = nodes, Source = "from", NodeID = "name", LinkGroup="group", fontFamily = "arial",
        NodeGroup = "group", Target = "to", Value = "value", fontSize = 12, nodeWidth = 100, sinksRight = F,
        colourScale=my_color, width = 1500, height = 600)
}

saveNetwork(plot_sankey_networkD3(I_genes), "Figure4D.html")
```
