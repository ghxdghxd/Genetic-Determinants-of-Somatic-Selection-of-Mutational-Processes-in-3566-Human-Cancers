# intergrated regression

## prepare input data

```R
#!/usr/bin/env Rscript
# @Date    : 2018/11/19 14:43:58
# @Author  : JT Guo
# @Email   : guojt-4451@163.com

# 00_preAnalysis.R for ~/Projects/TCGA/hg19/workflow/12rare

source("~/Projects/TCGA/hg19/Functions.R")
library(parallel)
library(data.table)
library(GeneGeneInteR)
library(methods)
library(stringr)
## genotypes

load("/share/data4/TCGA/Germline/split_vcf/norm/exonic/burden.50.RData")

## genetic burden in panCancer

miss = sapply(names(burden), function(x){
    m = burden[[x]]
    return(m[,grep("miss", colnames(m))])
    })
miss = do.call(cbind, miss)

trun = sapply(names(burden), function(x){
    m = burden[[x]]
    return(m[,grep("trun", colnames(m))])
    })
trun = do.call(cbind, trun)

large = sapply(names(burden), function(x){
    m = burden[[x]]
    return(m[,grep("large", colnames(m))])
    })
large = do.call(cbind, large)

## genetic burden in cancer-specific

miss_split = as.data.frame(do.call(rbind, mclapply(unique(mem$cancer), function(x){
    return(apply(miss[mem$cancer == x, ], 2, rank))
}, mc.cores = 15)))
miss_split = miss_split[rownames(miss), ]
trun_split = as.data.frame(do.call(rbind, mclapply(unique(mem$cancer), function(x){
    return(apply(trun[mem$cancer == x, ], 2, rank))
}, mc.cores = 15)))
trun_split = trun_split[rownames(trun), ]
large_split = as.data.frame(do.call(rbind, mclapply(unique(mem$cancer), function(x){
    return(apply(large[mem$cancer == x, ], 2, rank))
}, mc.cores = 15)))
large_split = large_split[rownames(large), ]
miss = as.data.frame(apply(miss, 2, rank))
trun = as.data.frame(apply(trun, 2, rank))
large = as.data.frame(apply(large, 2, rank))

# somatic Nsy-mutations
load("~/Projects/TCGA/hg19/workflow/01maf/mem.snv.RData")
mem = mem[which(mem$pop == "CEU"), ]
mem$cancer[grep('BRCA', mem$cancer)] = "BRCA"
int = intersect(rownames(mem), rownames(miss))
mem = mem[int, ]
snv = as.data.frame(snv[int, ])
rownames(snv) = int

## Somatic mutational propensity
load("~/Projects/TCGA/hg19/pmSignatre/SNV10/sig7.RData")
logSigMat = pmLogSigMat[int, ]

## SCNAs

load("~/Projects/TCGA/hg19/CNV/Gistic2_CopyNumber_by_genes/thresholded/all.cnv.RData")
cnvT = as.data.frame(cnv[match(int, rownames(cnv)), ])
rownames(cnvT) = int

cnvT3 = data.matrix(cnvT)
cnvT3[cnvT3 == 2] = 1
cnvT3[cnvT3 == -2] = -1
cnvT3 = as.data.frame(cnvT3, stringsAsFactors = F)

## TSS methylation levels

load("~/Projects/TCGA/hg19/workflow/04CNV_methy_mRNA/methy/396065_probes/all.gene.RData")
methy = as.data.frame(methy[match(int, rownames(methy)), ])
rownames(methy) = int

methy0 = methy
methy0[is.na(methy0)] = 0

methyT = as.matrix(methy)
methyT[which(methyT >= 0.8)] = 2
methyT[which(methyT >= 0.2 & methyT < 0.8)] = 1
methyT[which(methyT < 0.2)] = 0
methyT = as.data.frame(methyT)

methyT0 = methyT
methyT0[is.na(methyT0)] = 0

## gene fusion statuses

load("/share/data4/TCGA/TCGA_fusion/fusion_mat.RData")
fusion_mat = fusion_mat[match(rownames(miss), rownames(fusion_mat)),]
rownames(fusion_mat) = rownames(miss)
fusion_mat[is.na(fusion_mat)] = 0

## clinical infomations

load("~/Projects/TCGA/hg19/workflow/12rare/clin.RData")
clin = clin[match(int, rownames(clin)), c("age_at_initial_pathologic_diagnosis", "gender", "stage", "purity")]
colnames(clin) = c("age","gender","stage", "purity")
clin$gender = as.factor(clin$gender)
clin$stage = as.factor(clin$stage)

## the fraction of immune cells

immu <- read.table("~/Projects/signature_immune/The_Immune_Landscape_of_Cancer_Figure1.csv",header = T,sep="\t",stringsAsFactors = F)
rownames(immu) = immu$TCGA.Participant.Barcode
immu$TCGA.Participant.Barcode = NULL
immu = immu[int, ]

save(logSigMat, snv, cnvT3, methyT0, fusion_mat,
    miss, miss_split, trun, trun_split, large, large_split,
    mem, clin, immu, file = "all.50.flipped.run.RData")
```

## run regression

```sh
Rscript run_regression.R --all_RData all.50.flipped.run.RData --out_name mat_rank --processNumber 45
```

## selecting candidated genes

```R
library(dplyr)
library(ggplot2)
load("mat_rank.all.RData")
all_mat[[1]] = all_mat[[1]] %>% filter(cancer=="panCan")
mat = do.call(rbind, all_mat)
rownames(mat) = NULL

mat = mat %>% group_by(gene, sig, cancer) %>%
    filter(length(intersect(c("miss", "trun", "large"), feature)) > 0 & length(intersect(c("snv", "scna", "methy", "fusion"), feature)) > 0) %>%
    ungroup() %>% group_by(cancer) %>% mutate(FDR = p.adjust(P, method = "fdr"), lmFDR = p.adjust(lmP, method = "fdr"))

mat_FDR = mat[mat$FDR < 0.1 & mat$lmFDR < 0.1, ]

mat_sig = mat_FDR %>% group_by(gene, sig, cancer) %>%
    filter((length(intersect(c("miss", "trun", "large"), feature)) > 0 & length(intersect(c("snv", "scna", "methy", "fusion"), feature)) > 0))

write.table(mat_FDR,"mat_FDR.txt",r=F,c=T,sep="\t",quote=F)
write.table(mat_sig,"mat_sig.txt",r=F,c=T,sep="\t",quote=F)
```
