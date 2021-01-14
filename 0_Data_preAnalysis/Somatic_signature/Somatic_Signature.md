# Get Mutational signature (MS)

## prepare input file

```R
#!/usr/bin/env Rscript
# @Date    : 2018/11/8 下午4:21:25
# @Author  : JT Guo
# @Email   : guojt-4451@163.com

# 03_pmSingature.R

# Args <- commandArgs(T)
load("~/Projects/TCGA/hg19/workflow/01maf/all.maf.filtered.RData")
write.table(maf[maf$Variant_Type == "SNP", c("Tumor_Sample_Barcode", "Chromosome",
  "Start_position", "Reference_Allele", "Tumor_Seq_Allele2")], "DP30.AD10.SNP10.input.pm",
  row.names = F, col.names = F, sep = "\t", quote = F)
```

## get mutational signatures

```sh
Rscript ~/Projects/TCGA/hg19/pmSignatre/getSig.all.R \
   ~/Projects/TCGA/hg19/pmSignatre/DP30.AD10.SNV10.input.pm \
   10 20 ~/Projects/TCGA/hg19/pmSignatre/SNV10

Rscript est.R ~/Projects/TCGA/hg19/pmSignatre/SNV10
```

# get mutational propensity

```R
load("~/service/hpc/Projects/TCGA/hg19/workflow/02pmsignature/SNV10/sig7.RData")
getLogMat = function(sigMat, BG, trans = "log"){
    sigMat = as.matrix(sigMat)
    sigMat[which(sigMat == 0)] = NA
    if(trans == 'log'){
        sigMat = log(sigMat/sigMat[, BG])
    }else if(trans == 'log2'){
        sigMat = log2(sigMat/sigMat[, BG])
    }
    sigMat = sigMat[, -grep(BG, colnames(sigMat))]
    sigMat = sigMat[which(rowSums(sigMat, na.rm = T) != 0), which(colSums(sigMat, na.rm = T) != 0)]
    return(as.data.frame(sigMat))
}

MS[MS < 1e-3] = 0
MS = MS[which(MS$BG > 0), ]
MP = getLogMat(MS, "BG", trans = 'log')
```
