# Gene somatic functional statuses

## filter samples

```R
#!/usr/bin/env Rscript
# @Date    : 2018/11/8 下午5:12:50
# @Author  : JT Guo
# @Email   : guojt-4451@163.com

# run.R

# Args <- commandArgs(T)

# methylation
library(data.table)
library(parallel)
library(robustbase)
files = list.files("~/Projects/TCGA/hg19/CpG", pattern="*.CpG.RData", full.names=T)
mc=30
info <- read.csv("/share/data0/reference/methylation/Methyl450CpgIslandDetails.hg19.txt", header = T, 
    sep = ",", stringsAsFactor = F)
probes <- read.csv("/share/data0/reference/methylation/probes_to_gene_hg19.txt", header = T, 
    sep = "\t", stringsAsFactor = F)
probes = probes[which(probes$UCSC_RefGene_Group %in% c("TSS1500", "TSS200")), ]
probes = probes[-which(info$CHR[match(probes$Name, info$Name)] %in% c("X","Y")), ]

gene = unique(probes$UCSC_RefGene_Name)

loadMethy <- function(file) {
    load(file)
    met = data.matrix(met)
    methy <- mclapply(gene, function(y){
      print(y)
      int = intersect(colnames(met), probes$Name[which(probes$UCSC_RefGene_Name == y)])
      if(length(int) == 1){
        return(met[, int])
      }else{
        # return(rowMedians(met[, int], na.rm=T))
        return(apply(met[, int], 1, function(x){return(max(x, na.rm = T))}))
      }}, mc.cores = mc)
    methy <- do.call('rbind', methy)
    rownames(methy) <- gene
    methy <- t(data.frame(methy, stringsAsFactors = F, check.names = F))
    # rownames(methy) = gsub("-01$", "", rownames(methy))
    rownames(methy) = substr(rownames(methy), 1, 12)
    save(methy, file = gsub("\\.CpG", ".gene.TSS.CpG", file))
}

for(file in files){
    loadMethy(file)
}


```
