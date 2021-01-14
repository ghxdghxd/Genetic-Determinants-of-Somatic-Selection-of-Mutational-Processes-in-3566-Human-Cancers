#!/usr/bin/env Rscript
# @Date    : 2018/12/6 10:16:04
# @Author  : JT Guo
# @Email   : guojt-4451@163.com

# 00_preanalysis.R

# Args <- commandArgs(T)
Args <- commandArgs(T)
cancer = Args[1]
chrom = Args[2]
mc = 48
# dir.create(as.character(chrom), showWarnings=F)
# setwd(as.character(chrom))

library(parallel)
library(methods)
library(stringr)
library(data.table)
library(GeneGeneInteR)

if(cancer == "panCancer"){
    genoFile = paste0("~/Projects/TCGA/hg19/workflow/03genotype/02imputeALL/new_info0.9_maxPro0.9/4CEU_com_geno/", chrom, ".geno.com.gz")
}else{
    genoFile = paste0("~/Projects/TCGA/hg19/workflow/03genotype/02imputeALL/new_info0.9_maxPro0.9/4CEU_com_geno/", cancer, "/", chrom, ".geno.com.gz")
}

geno = as.data.frame(fread(genoFile, header=T,sep="\t",stringsAsFactors=F, check.names=F), stringsAsFactors = F)
rownames(geno) = geno$ID
geno$ID = NULL

load("~/Projects/TCGA/hg19/workflow/01maf/mem.RData")
mem = mem[which(mem$pop == "CEU"), ]
mem$cancer[grep('BRCA', mem$cancer)] = "BRCA"

# load mutational propensity (MP)
load("~/Projects/TCGA/hg19/pmSignatre/SNV10/sig7.RData")


sample_int = intersect(colnames(geno), rownames(MP))
geno = geno[, sample_int]
geno = t(geno)

MP = MP[sample_int, ]

mem = mem[sample_int, ]

rmWHE <- function(geno, MAF = 0.05, HWE = 1e-5, CallRate = 0.95){
    new_geno <- new("SnpMatrix", as.matrix(geno) + 1)
    new_geno1 <- snpMatrixScour(new_geno, min.maf = MAF, min.eq = HWE, call.rate = CallRate)
    rs = colnames(as.matrix(new_geno1$snpX))
    return(rs)
}

snp = rmWHE(geno, MAF = 0.05, HWE = 1e-5, CallRate = 0.95)
geno = geno[, which(colnames(geno) %in% snp)]
geno = as.data.frame(geno)
gc()
save(geno, logSigMat, file = gsub("gz","RData", basename(genoFile)))

