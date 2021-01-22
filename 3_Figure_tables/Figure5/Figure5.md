# Figure 5

## Figure 5A

```R
file_Candidate_gene = "../TableS2/Table_S2A.csv"
file_E_gene = "../TableS2/Table_S2B.csv"
file_I_gene = "../TableS2/Table_S2C.csv"
file_cosmic_gene = "../Figure4/Cancer_Gene_Census_allThu_Mar_14_06_02_23_2019.csv"
file_susceptibility_gene = "../Figure4/Cancer-susceptibility-genes.txt"
file_Compendium_Cancer_Gene = "../Figure4/Compendium_Cancer_Genes.tsv"
file_Nucleotide_context_gene = "../Figure4/Nucleotide_context_gene.csv"

eff_RData = "DepMap_19q2_dependency_score.RData"

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

library(GeneOverlap)
get_bin_gene_enrich = function(eff_raw, gene_list, type, index){
    genome_size = length(unique(eff_raw$variable))
    gene_list = intersect(eff_raw$variable, gene_list)
    bin_gene_enrich = mclapply(index, function(x){
        print(x)
        eff_raw_bin = eff_raw[eff_raw$value < x[2] & eff_raw$value > x[1], ]
        res = testGeneOverlap(newGeneOverlap(unique(as.character(eff_raw_bin$variable)), gene_list, genome.size = genome_size))
        FC = (res@cont.tbl[4]/(res@cont.tbl[2] + res@cont.tbl[4]))/(res@cont.tbl[3]/(res@cont.tbl[1] + res@cont.tbl[3]))
        m = c(p = res@pval, OR = res@odds.ratio, FC = FC, nBnG_nBG_BnG_BG = paste(res@cont.tbl, collapse="_"), type = type)
        return(m)
    }, mc.cores = 20)
    bin_gene_enrich = as.data.frame(do.call(rbind, bin_gene_enrich), stringsAsFactors=F)
    bin_gene_enrich$ID = gsub("c","",paste(index))
    bin_gene_enrich$p = as.numeric(bin_gene_enrich$p)
    bin_gene_enrich$OR = as.numeric(bin_gene_enrich$OR)
    bin_gene_enrich$FC = as.numeric(bin_gene_enrich$FC)
    return(bin_gene_enrich)
}

plot_bin_gene_enrich <- function(eff_raw, index, xlabs = NULL){
  index = lapply(seq(floor(min(eff_raw$value)), ceiling(max(eff_raw$value)), by = index), function(x){return(c(x, x+index))})
  bin_gene_enrich_IV_sig = get_bin_gene_enrich(eff_raw, E_genes$gene, "E-genes", index)
  bin_gene_enrich_immu_sig = get_bin_gene_enrich(eff_raw, I_genes$gene, "I-genes", index)
  a = as.data.frame(do.call(rbind, apply(CGCs, 1, function(x){
    if(x["Role.in.Cancer"] == ""){
      type = ""
    }else{
      type = strsplit(x["Role.in.Cancer"], split = ",")[[1]]
    }
    return(cbind(gene = x["Gene.Symbol"], type = unlist(type)))
    })), stringsAsFactors = F)
  a$type = gsub(" ", "", a$type)
  a = a[a$type %in% c("TSG", "oncogene"), ]
  bin_gene_enrich_COSMIC = get_bin_gene_enrich(eff_raw, CGCs$Gene.Symbol, "COSMIC", index)
  bin_gene_enrich_genes152 = get_bin_gene_enrich(eff_raw, Susceptibility_genes$Gene..Symbol, "Susceptibility_Genes", index)
  bin_gene_enrich_Compendium_Cancer_Genes = get_bin_gene_enrich(eff_raw, unique(Compendium_Cancer_Genes$SYMBOL), "Compendium_Cancer_Genes", index)
  bin_gene_enrich_TableS4_460_gene = get_bin_gene_enrich(eff_raw, Nucleotide_context_genes$Gene, "Nucleotide_context", index)
  bin_gene_enrich = rbind(bin_gene_enrich_IV_sig, bin_gene_enrich_immu_sig,
      bin_gene_enrich_COSMIC, bin_gene_enrich_Compendium_Cancer_Genes,
      bin_gene_enrich_TableS4_460_gene, bin_gene_enrich_genes152)

  bin_gene_enrich = bin_gene_enrich %>% mutate(lab = ifelse(p < 0.001,"*\n",""))
  bin_gene_enrich = bin_gene_enrich %>% mutate(lab = ifelse(p < 0.01, paste0(lab, "*\n"),""))
  bin_gene_enrich = bin_gene_enrich %>% mutate(lab = ifelse(p < 0.05, paste0(lab, "*"),""))
  bin_gene_enrich = bin_gene_enrich %>% group_by(ID) %>% mutate(sum = sum(OR)) %>% filter(sum > 0)
  bin_gene_enrich$ID = factor(bin_gene_enrich$ID, levels = unique(bin_gene_enrich$ID))
  bin_gene_enrich$type = factor(bin_gene_enrich$type, levels = unique(bin_gene_enrich$type))
  colors = data.frame(color = c("#ED751B", "#F4BC1B", "#B0C33D", "#6BA732", "#18855B", "#104287"),
                      method = unique(bin_gene_enrich$type), stringsAsFactors = F)
  p = ggplot(bin_gene_enrich, aes(x = ID, y = FC, group = type)) +
    geom_bar(aes(fill = type), stat = "identity", position = position_dodge(width=1)) +
    geom_text(aes(label = lab), position = position_dodge(width=0.9), vjust=1, hjust = 0.5, size = 6) +
    facet_grid(~ID, scales = "free", space = "free") +
    labs(x = xlabs, y = "Fold enrichment ", fill = "") +
    scale_fill_manual(breaks = colors$method, values = colors$color) +
    theme_gray() +
    theme(panel.grid = element_blank(),
          strip.text = element_blank(),
          legend.position = c(0.5, 0.8),
          legend.background = element_blank(),
          strip.background = element_blank(),
          plot.background = element_blank())
  return(p)
}

load(eff_RData)
eff_raw = melt(depmap_eff, id.vars="cancer")
index_0.2 = lapply(seq(-4, -0.2, by = 0.2), function(x){return(c(x, x+0.2))})
eff_raw_median = eff_raw %>% group_by(variable) %>% summarise(value = median(value))

pdf("Figure5A.pdf", width = 10, height = 3)
plot_bin_gene_enrich(eff_raw_median, 0.2, xlabs = "Median gene dependence scores")
dev.off()
```

## Figure 4B

```R
load("IC50_associated_genes.RData")

get_gene_num = function(genes, panCan = T){
    genes = unique(genes)
    if(panCan){
        m = CCLE_lm_IC50_mRNA_sig %>% filter(cancer == "panCan", gene %in% genes)
        cancer = "panCan"
    }else{
         m = CCLE_lm_IC50_mRNA_sig %>% filter(cancer != "panCan", gene %in% genes)
         cancer = "cancer_specific"
    }
    int = intersect(m$gene, genes)
    return(cbind(length(int), length(genes), cancer))
}

gene_num = as.data.frame(rbind(cbind(get_gene_num(E_genes$gene[E_genes$cancer_type == "pan-cancer"]), type = "E-genes(pan-Cancer)"),
    cbind(get_gene_num(I_genes$gene[I_genes$cancer_type == "pan-cancer"]), type = "I-genes(pan-Cancer)"),
    cbind(get_gene_num(CGCs$Gene.Symbol), type = "CGCs"),
    cbind(get_gene_num(Compendium_Cancer_Genes$SYMBOL), type = "Compendium_Cancer_Genes"),
    cbind(get_gene_num(Susceptibility_genes$Gene..Symbol), type = "Susceptibility_genes"),
    cbind(get_gene_num(Nucleotide_context_genes$Gene), type = "Nucleotide_context_genes")), stringsAsFactors=F)
gene_num$V1 = as.numeric(gene_num$V1)
gene_num$V2 = as.numeric(gene_num$V2)
gene_num$rate = gene_num$V1/gene_num$V2

pdf("Figure5B.pdf", width = 7, height = 5)
as.data.frame(gene_num) %>% arrange(rate) %>% mutate(type = factor(type, levels = rev(unique(type)))) %>%
    ggplot(aes(type, rate)) + geom_bar(stat = "identity", position = "dodge", width = 0.8) +
    labs(y="Fraction of genes", x = "", title = "Genes associated with drug IC50\nin all cell lines of CCLE") +
    coord_cartesian(ylim = c(0.63, 0.71)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = -30, hjust = 0, vjust = 1, size = 14),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 14),
          legend.position = "none",
          aspect.ratio = 1)
dev.off()
```

# Figure 5C

```R
load("IC50_associated_genes.RData")

gdsc = read.table(gzfile('depmap_GDSC1_fitted_dose_response_17Jul19.txt.gz'),header=T,sep="\t",stringsAsFactors=F)
drug_info = unique(gdsc[, c("DRUG_ID","DRUG_NAME","PUTATIVE_TARGET","PATHWAY_NAME")])
library(wesanderson)

pdf("Figure5C.pdf", width = 11, height = 5)
CCLE_lm_IC50_mRNA_IV %>% filter(cancer=="panCan") %>% group_by(PATHWAY_NAME, sig) %>%
    mutate(gene_rate = length(unique(gene))/length(unique(E_genes$gene[E_genes$mutational_propensity == unique(sig)])),
           drug_rate = length(unique(DRUG_NAME))/length(unique(drug_info$DRUG_NAME[drug_info$PATHWAY_NAME==unique(PATHWAY_NAME)])),
           pathway = paste0(unique(PATHWAY_NAME), "(", length(unique(drug_info$DRUG_NAME[drug_info$PATHWAY_NAME==unique(PATHWAY_NAME)])), ")"),
           stanBeta = mean(abs(stanBeta))) %>% ungroup %>%
    select(sig, pathway, drug_rate, stanBeta) %>% unique() %>%
    arrange(stanBeta) %>% mutate(pathway = factor(pathway, levels = rev(unique(pathway)))) %>%
    ggplot(aes(pathway, sig, color = stanBeta, size = drug_rate)) + geom_point() +
    labs(x="",y="", color = "Average absolute effect size\nof MP-genes on drugs\nof targeted pathways",
         size = "Fraction of drugs\nassociated with\nE-genes' expression\nin each targeted pathways") +
    scale_color_fermenter(direction = 1, n.breaks = 4) +
    facet_grid(sig ~ pathway, scales = "free", space = "free") +
    scale_size(range = c(1, 8)) +
    theme(axis.text.x = element_text(angle = -30, hjust = 0, vjust = 1, size = 12),
        axis.text.y = element_text(size = 10),
          axis.ticks = element_blank(),
          strip.text = element_blank(),
          strip.background = element_blank(),
          legend.direction = "horizontal",
          legend.box = "vertical",
          panel.spacing = unit(0.1, 'lines'),
          panel.grid.major = element_blank(),
          aspect.ratio = 1)
dev.off()
```

## Figure 5D and 5E

```R
# Melanoma: Ratio of the interferon-γ signature to the immunosup-pression signature predicts anti-PD-1 therapy response in melanoma
# Gastric: Comprehensive molecular characterization of clinical responses to PD-1 inhibition in metastatic gastric cancer
ICI <- read.table("ICI_result.csv", header = T, sep=",", stringsAsFactors = F)
ICI_sig = ICI[which(ICI$p_value < 0.05), ]

get_gene_num_ICI = function(genes, panCan = T){
    genes = unique(genes)
    if(panCan){
        cancer = "panCan"
    }else{
         cancer = "cancer_specific"
    }
    return(t(sapply(unique(ICI_sig$cancer), function(x){
        int = intersect(ICI_sig$gene[ICI_sig$cancer == x], genes)
        return(cbind(length(int), length(genes), cancer, x))
    })))
}

gene_num_ICI = as.data.frame(rbind(cbind(get_gene_num_ICI(E_genes$gene[E_genes$cancer_type == "pan-cancer"]), type = "E-genes(pan-Cancer)"),
    cbind(get_gene_num_ICI(I_genes$gene[I_genes$cancer_type == "pan-cancer"]), type = "I-genes(pan-Cancer)"),
    cbind(get_gene_num_ICI(CGCs$Gene.Symbol), type = "CGCs"),
    cbind(get_gene_num_ICI(Compendium_Cancer_Genes$SYMBOL), type = "Compendium_Cancer_Genes"),
    cbind(get_gene_num_ICI(Susceptibility_genes$Gene..Symbol), type = "Susceptibility_genes"),
    cbind(get_gene_num_ICI(Nucleotide_context_genes$Gene), type = "Nucleotide_context_genes")), stringsAsFactors=F)
gene_num_ICI$V1 = as.numeric(gene_num_ICI$V1)
gene_num_ICI$V2 = as.numeric(gene_num_ICI$V2)
gene_num_ICI$rate = gene_num_ICI$V1/gene_num_ICI$V2

pdf("Figure5D_E.pdf", width = 8, height = 5)
as.data.frame(gene_num_ICI) %>% filter(V3 == "panCan", V4 == "Melanoma") %>%
            arrange(rate) %>% mutate(type = factor(type, levels = rev(unique(type)))) %>%
            ggplot(aes(type, rate)) + geom_bar(stat = "identity", position = "dodge", width = 0.8) +
            labs(y="Fraction of genes", x = "", title = "Genes associated with anti-PD-1\ntherapy response in melanoma") +
            coord_cartesian(ylim = c(0.092, 0.1)) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = -30, hjust = 0, vjust = 1, size = 14),
                axis.text.y = element_text(size = 12),
                axis.title.y = element_text(size = 14),
                legend.position = "none",
                aspect.ratio = 1)
as.data.frame(gene_num_ICI) %>% filter(V3 == "panCan", V4 == "Gastric") %>%
            arrange(rate) %>% mutate(type = factor(type, levels = rev(unique(type)))) %>%
            ggplot(aes(type, rate)) + geom_bar(stat = "identity", position = "dodge", width = 0.8) +
            labs(y = "Fraction of genes", x = "", title = "Genes associated with anti-PD-1\ntherapy response in metastatic gastric cancer") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = -30, hjust = 0, vjust = 1, size = 14),
                axis.text.y = element_text(size = 12),
                axis.title.y = element_text(size = 14),
                legend.position = "none",
                aspect.ratio = 1)
dev.off()
```
