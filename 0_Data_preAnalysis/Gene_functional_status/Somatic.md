# Gene somatic functional statuses

## Somatic mutations

### filter samples

```R
#!/usr/bin/env Rscript
# @Date    : 2018/11/8 下午3:15:57
# @Author  : JT Guo
# @Email   : guojt-4451@163.com

# merge_maf.R

# Args <- commandArgs(T)

library(stringr)
library(parallel)
files = list.files("~/Projects/TCGA/hg19/00BroadMaf", pattern="*maf$", full.names=T)

readMaf = function(x) {
    maf <- read.csv(x, header=T, sep="\t", stringsAsFactors=F, comment.char="#")
    maf = data.frame(Cancer=strsplit(basename(x), "\\.")[[1]][1], maf, stringsAsFactors=F)
    return(maf)
}

maf = mclapply(files, readMaf, mc.cores=12)
maf = do.call(rbind, maf)

filterSample = function(){
    Sample = unique(maf$Tumor_Sample_Barcode)
    TriSample = substr(Sample, 1, 12)
    # duplicated sample，TCGA-E2-A1IE
    DupTriSample = unique(TriSample[duplicated(TriSample)], 1,12)
    # unique samples
    ID = Sample[-which(TriSample %in% DupTriSample)]
    # duplicated full-names
    DupSample = sapply(DupTriSample, function(x){list(Sample[grep(x, Sample)])})
    # selected Primary Solid Tumor
    DupSample = sapply(DupSample, function(x){
        if(length(x) > 1){
            y = x[which(substr(x, 14,15) == "01")]
            if(length(y) < 1){
                return(x)
            }else{
                return(y)
            }
        }else{
            return(x)
        }
    })
    # selected unique samples
    ID1 = as.character(DupSample[which(sapply(DupSample, length)==1)])
    DupSample = DupSample[-which(sapply(DupSample, length)==1)]
    # selected the Vial
    DupSample = sapply(DupSample, function(x){
        if(length(x) > 1){
            y = x[which(substr(x, 16,16) == "A")]
            if(length(y) < 1){
                return(x)
            }else{
                return(y)
            }
        }else{
            return(x)
        }
    })
    ID2 = as.character(DupSample[which(sapply(DupSample, length)==1)])
    DupSample = DupSample[-which(sapply(DupSample, length)==1)]
    # selected DNA samples
    DupSample = sapply(DupSample, function(x){
        if(length(x) > 1){
            y = x[which(substr(x, 20, 20) == "D")]
            if(length(y) < 1){
                return(x)
            }else{
                return(y)
            }
        }else{
            return(x)
        }
    })

    ID3 = as.character(DupSample[which(sapply(DupSample, length)==1)])
    DupSample = DupSample[-which(sapply(DupSample, length)==1)]
    # selected Portion 18-19: 01
    DupSample = sapply(DupSample, function(x){
        if(length(x) > 1){
            y = x[which(substr(x, 18, 19) == "01")]
            if(length(y) < 1){
                return(x)
            }else{
                return(y)
            }
        }else{
            return(x)
        }
    })
    ID4 = as.character(DupSample[which(sapply(DupSample, length)==1)])
    DupSample = DupSample[-which(sapply(DupSample, length)==1)]

    # selected the higer plate 22-25
    DupSample = sapply(DupSample, function(x){
        return(sort(x)[length(x)])
        })
    return(c(ID, ID1,ID2,ID3,ID4,DupSample))
}

sample = filterSample() # filter hg19.maf

maf = maf[which(maf$Tumor_Sample_Barcode %in% sample), ]
maf$t_depth = maf$t_ref_count + maf$t_alt_count
maf = maf[which(maf$t_depth > 30 & maf$t_alt_count > 10), ]
maf = maf[which(maf$Chromosome != "MT"), ]
maf$ID = substr(maf$Tumor_Sample_Barcode, 1, 12)
```

### get Somatic mutation

```R
getMutSpetra = function(maf, strand = F){
  m_raw = as.data.frame.matrix(table(maf[, c("Tumor_Sample_Barcode", "Variant_Type")]))
  colnames(m_raw) = gsub("SNP","SNV",colnames(m_raw))
  m_raw$varSum = rowSums(m_raw)
  m <- maf[which(maf$Variant_Type=="SNP"), c("Tumor_Sample_Barcode","ref_context", "Tumor_Seq_Allele2")]
  m$ref_context <- substr(m$ref_context, 10, 12)
  m$ref_context <- toupper(paste(m$ref_context, m$Tumor_Seq_Allele2, sep = ""))
  changeBase <- function(x){
    if(x == "A"){return("T")}else if(x == "C"){return("G")}else if(x == "G"){return("C")}else{return("A")}
  }
  m$ref_context <- sapply(m$ref_context, function(x){
    name = unlist(strsplit(x, ""))
    if(name[2] %in% c("C", "T")){
      return(paste(paste(name[2], name[4], sep = ""), paste(name[1], name[3], sep = "."), "t"))
    }else{
      name = sapply(name, changeBase)
      return(paste(paste(name[2], name[4], sep = ""), paste(name[3], name[1],sep = "."), "u"))
    }
  })
  m$alteration = sub("([ACGTN])([ACGTN]) .+", "\\1>\\2", m$ref_context)
  if(strand){
    m$context = sub("([ACGTN])([ACGTN]) ([ACGTN]).([ACGTN]) ([ut])", "\\3[\\1>\\2]\\4\\_\\5", m$ref_context)
  }else{
    m$context = sub("([ACGTN])([ACGTN]) ([ACGTN]).([ACGTN]) ([ut])", "\\3[\\1>\\2]\\4", m$ref_context)
  }
  m_sub <- as.data.frame.matrix(table(m[, c("Tumor_Sample_Barcode", "alteration")]), stringsAsFactors = F)
  m_trim <- as.data.frame.matrix(table(m[, c("Tumor_Sample_Barcode", "context")]), stringsAsFactors = F)
  return(cbind(sample = rownames(m_raw), m_raw, 
    m_sub[match(rownames(m_raw), rownames(m_sub)),], 
    m_trim[match(rownames(m_raw), rownames(m_sub)),]))
}

mem = na.omit(getMutSpetra(maf))
# colnames(mem) = gsub("SNP","SNV",colnames(mem))
mem$sample = substr(mem$sample, 1, 12)

mem = mem[which(mem$SNV > 10), ] ################ SNV > 10

maf = maf[which(maf$ID %in% mem$sample), ]

save(maf, mem, file="~/Projects/TCGA/hg19/workflow/01maf/all.maf.filtered.RData")
```

### get Somatic Nsy-mutations status

```R
maf$variant = sapply(maf$Variant_Classification, function(x) {
    if(x %in% c("Silent", "RNA")){
        return("Silent")
    }else{
        return("nonSilent")
    }
})

snv_count = as.data.frame.matrix(table(maf[which(maf$variant=="nonSilent"), c("ID","Hugo_Symbol")]))

# CEU population
pca <- read.table("~/Projects/TCGA/hg19/pop/PCA.txt",
    sep="\t", header=T, stringsAsFactors=F)
pca = pca[which(pca$pop == "Case"), ]
pca$pop[which(pca$PC2 > 27)] = "CHB"
pca$pop[which(pca$PC1 < 2 & pca$PC2 < 2)] = "CEU"
pca$pop[which(pca$PC1 > 15 & pca$PC1 < 38 & pca$PC2 < 0)] = "ASW"
pca$pop[which(pca$PC1 > 38 & pca$PC2 < 0)] = "YRI"
pca$pop[which(pca$PC1 < 5 & pca$PC2 < 27 & pca$PC2 > 2)] = "MXL"
pca$pop[which(pca$pop == "Case")] = NA
rownames(pca) = substr(rownames(pca), 1, 12)

mem = data.frame(cancer = maf$Cancer[match(mem$sample, maf$ID)], 
    pop = pca$pop[match(mem$sample, rownames(pca))], mem, 
    stringsAsFactors=F, check.names=F)
mem$ID = rownames(mem)
rownames(mem) = mem$sample

escc <- read.table("~/Projects/TCGA/TCGA.ESCC.txt", header = T, sep = "\t", stringsAsFactors = F)
mem$cancer[which(mem$sample %in% escc$barcode)] = 'ESCC'

save(mem, snv, pca, file = "mem.snv.pca.RData")
```

## gene fusion status

> download the data from https://tumorfusions.org/

```R
library(stringr)
library(dplyr)
fusion <- read.table("TCGA/TCGA_fusion/TCGA_cancer_fusion.txt",header=T,sep="\t",stringsAsFactors=F)

a = as.data.frame(t(apply(fusion, 1, function(x){
    gene = sort(c(x["Gene_A"], x["Gene_B"]))
    x["Gene_A"] = gene[1]
    x["Gene_B"] = gene[2]
    return(x)
})), stringsAsFactors = F)

fusion_sub = unique(a[,c("Tissue","Sample","Gene_A","Gene_B")])

filterSample = function(){
    Sample = unique(fusion$Sample)
    TriSample = substr(Sample, 1, 12)
    # duplicate shor-name，TCGA-E2-A1IE
    DupTriSample = unique(TriSample[duplicated(TriSample)], 1,12)
    # unique samples
    ID = Sample[-which(TriSample %in% DupTriSample)]
    # dupilcated full-name
    DupSample = sapply(DupTriSample, function(x){list(Sample[grep(x, Sample)])})
    # selected Primary Solid Tumor
    DupSample = sapply(DupSample, function(x){
        if(length(x) > 1){
            y = x[which(substr(x, 14,15) == "01")]
            if(length(y) < 1){
                return(x)
            }else{
                return(y)
            }
        }else{
            return(x)
        }
    })

    ID1 = as.character(DupSample[which(sapply(DupSample, length)==1)])
    DupSample = DupSample[-which(sapply(DupSample, length)==1)]
    # selected the Vial
    DupSample = sapply(DupSample, function(x){
        if(length(x) > 1){
            y = x[which(substr(x, 16,16) == "A")]
            if(length(y) < 1){
                return(x)
            }else{
                return(y)
            }
        }else{
            return(x)
        }
    })
    ID2 = as.character(DupSample[which(sapply(DupSample, length)==1)])
    DupSample = DupSample[-which(sapply(DupSample, length)==1)]
    # selected the DNA sample
    DupSample = sapply(DupSample, function(x){
        if(length(x) > 1){
            y = x[which(substr(x, 20, 20) == "")]
            if(length(y) < 1){
                return(x)
            }else{
                return(y)
            }
        }else{
            return(x)
        }
    })
    ID3 = as.character(DupSample[which(sapply(DupSample, length)==1)])
    return(c(ID, ID1,ID2,ID3))
}

samples = filterSample()
fusion_sub = fusion_sub[fusion_sub$Sample %in% samples, ]
fusion_sub$ID = substr(fusion_sub$Sample, 1, 12)
genes = unique(c(fusion_sub$Gene_A, fusion_sub$Gene_B))
fusion_sub$Gene_A = factor(fusion_sub$Gene_A, levels = genes)
fusion_sub$Gene_B = factor(fusion_sub$Gene_B, levels = genes)

fusion_mat = as.data.frame.matrix(table(fusion_sub[,c("ID","Gene_A")])) + as.data.frame.matrix(table(fusion_sub[,c("ID","Gene_B")]))

save(fusion_mat, "/share/data4/TCGA/TCGA_fusion/fusion_mat.RData")
```
