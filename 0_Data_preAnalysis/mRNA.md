# filter the mRNA expressions

```R
####### TCGA immune paper
mRNA <- read.table("EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv.gz",header=T,sep="\t",stringsAsFactors=F)
genes = str_split(mRNA$gene_id, pattern="\\|", simplify=T)
mRNA = mRNA[which(genes[, 1] != "?"), ]
genes = genes[which(genes[, 1] != "?"), ]

colnames(mRNA) = gsub("\\.", "-", colnames(mRNA))
samples = colnames(mRNA)[-1]

filterSample = function(sample){
    Sample = unique(sample)
    TriSample = substr(Sample, 1, 12)
    DupTriSample = unique(TriSample[duplicated(TriSample)], 1,12)
    ID = Sample[-which(TriSample %in% DupTriSample)]
    DupSample = sapply(DupTriSample, function(x){list(Sample[grep(x, Sample)])})
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
    DupSample = DupSample[-c(1:2)]
    DupSample = sapply(DupSample, function(x){
        return(sort(x)[length(x)])
        })
    return(c(ID, ID1,ID2,DupSample))
}

samples = filterSample(samples)

mRNA = mRNA[, c("gene_id", samples)]

mRNA$gene_id = gsub("SLC35E2\\|728661", "SLC35E2B|728661", mRNA$gene_id)
mRNA$gene_id = gsub("SLC35E2\\|9906", "SLC35E2A|9906", mRNA$gene_id)

genes = str_split(mRNA$gene_id, pattern="\\|", simplify=T)
mRNA$gene_id = genes[,1]

rownames(mRNA) = mRNA$gene_id
mRNA$gene_id = NULL

mRNA = t(mRNA)

mRNA = log2(mRNA + 1)
save(mRNA, file = "pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.log2.RData")
```
