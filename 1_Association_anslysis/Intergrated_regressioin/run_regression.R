library("optparse")

option_list = list(
  make_option(c("--all_RData"), type = "character", default = NULL, help = "the all RData",
    metavar = "file", dest = "all_RData"),
  make_option(c("--out_name"), type = "character", default = NULL, help = "the output name",
    metavar = "str", dest = "out_name"),
  make_option(c("--processNumber"), type = "character", default = 1, help = "the number of process",
    metavar = "int", dest = "mc")
    )

args_parser = OptionParser(usage = "usage: 01_run_all [options]",
    option_list = option_list,
    add_help_option = TRUE, prog = NULL, description = "", epilogue = "")
args = parse_args(args_parser, print_help_and_exit = TRUE)

if (length(args) == 1){
  print_help(args_parser)
  q()
}
print(args)

load(args$all_RData)

genes = unique(c(colnames(miss), colnames(trun),colnames(large), colnames(snv), colnames(cnvT3), colnames(methyT0)))
genes = genes[!grepl("^HLA", genes)]
genes = genes[!grepl("^OR", genes)]
mRNA_gene = read.table("/share/data4/TCGA/TCGA_The_Immune_Landscape_of_Cancer/mRNA_gene.txt", header = F, stringsAsFactors = F)
genes = intersect(genes, mRNA_gene$V1)

library(parallel)
library(stringr)
library(tidyverse)

filter_feature <- function(data, feature){
    if(feature %in% colnames(data) & feature != "cancer"){
        t = table(data[, feature] != min(abs(data[, feature])))
        if(length(which(t >= 5)) == 1){
            data[, feature] = NULL
        }
    }
    return(data)
}

run_all = function(cancer, sig, gene, logSigMat, miss, trun, large, snv, scna, methy, fusion, mem, stepwise = F){
    dat = data.frame(sig = logSigMat[, sig],
        miss = miss[, which(colnames(miss) == gene)],
        trun = trun[, which(colnames(trun) == gene)],
        large = large[, which(colnames(large) == gene)],
        snv = snv[, which(colnames(snv) == gene)],
        scna = scna[, which(colnames(scna) == gene)],
        methy = methy[, which(colnames(methy) == gene)],
        fusion = fusion[, which(colnames(fusion) == gene)],
        cancer = mem$cancer, stringsAsFactors = F)
    if(cancer != "panCan"){
        dat = na.omit(dat[which(mem$cancer == cancer), ])
        dat$cancer = NULL
    }
    for(i in colnames(dat)[-1]){
        dat = filter_feature(dat, i)
    }
    if(!is.null(peer_factors)){
        dat = na.omit(cbind(dat, peer_factors[rownames(dat), ]))
        covariate = paste("+", paste(colnames(peer_factors), collapse = "+"))
        peer_num = ncol(peer_factors)
    }else{
        peer_num = 0
        covariate = NULL
    }
    if(nrow(dat) < 30){
        return(NULL)
    }
    if(cancer == "panCan"){
        fm = paste0("sig~", paste(colnames(dat)[2:(ncol(dat)- peer_num -1)], "cancer", sep="*", collapse="+"), covariate)
    }else{
        fm = "sig~."
    }
    if(stepwise){
        lm_res = summary(step(lm(fm, na.omit(dat)), trace = F))
    }else{
        lm_res = summary(lm(fm, dat))
    }
    # print(lm_res$call)
    res = as.data.frame(coefficients(lm_res))[, c(1, 4)]
    colnames(res) = c("Est", "P")
    res$r2 = lm_res$adj.r.squared
    index = which(rownames(res) %in% c("miss","trun","large","snv","scna","methy", "fusion"))
    if(length(index) > 0){
        res1 = res[index, ]
        res1$cancer = cancer
        res1$feature = rownames(res1)
    }else{
        res1 = NULL
    }
    res2 = res[grepl(":", rownames(res)), ]
    if(nrow(res2) != 0){
        a = t(sapply(rownames(res2), function(x){return(sort(strsplit(x,split=":")[[1]]))}))
        res2$cancer = gsub("cancer", "", a[, 1])
        res2$feature = a[,2]
        res = rbind(res1, res2)
    }else{
        res = res1
    }
    if(is.null(res)){
        return(res)
    }
    res = data.frame(gene = gene, sig = sig, res, lmP = 1 - pf(lm_res$fstatistic[1], lm_res$fstatistic[2], lm_res$fstatistic[3]), 
        stringsAsFactors = F, row.names = NULL)
    return(res)
}

all_mat = lapply(c("panCan", unique(mem$cancer)), function(cancer){
    if(cancer == "panCan"){
        mat = lapply(colnames(logSigMat), function(x){
            res = mclapply(genes, function(y){
                print(paste(cancer, x, y))
                # return(run_all(x, y))
                return(run_all(cancer, x, y, logSigMat, miss, trun, large, snv, cnvT3, methyT0, fusion_mat, mem, args$stepwise))
            }, mc.cores = args$mc)
            res = do.call(rbind, res)
            return(res)
        })
    }else{
        mat = lapply(colnames(logSigMat), function(x){
            res = mclapply(genes, function(y){
                print(paste(cancer, x, y))
                # return(run_all(x, y))
                return(run_all(cancer, x, y, logSigMat, miss_split, trun_split, large_split, snv, cnvT3, methyT0, fusion_mat, mem, args$stepwise))
            }, mc.cores = args$mc)
            res = do.call(rbind, res)
            return(res)
        })
    }
    mat = do.call(rbind, mat)
    mat = na.omit(mat)
    return(mat)
})
names(all_mat) = c("panCan", unique(mem$cancer))
save(all_mat, file = paste0(args$out_name, ".all.RData"))

