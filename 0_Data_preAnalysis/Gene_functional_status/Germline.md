# Download TCGA WES germline variant

> Pathogenic Germline Variants in 10,389 Adult Cancers
> /share/data4/TCGA/Germline/PCA.r1.TCGAbarcode.merge.tnSwapCorrected.10389.vcf.gz

10389 cancers contains 21,121,543 germline variants.

## split vcf by Chr, extract genotype

> /share/data4/TCGA/Germline/split_vcf

```sh
for i in 5;
do
{
tabix -h ../../PCA.r1.TCGAbarcode.merge.tnSwapCorrected.10389.vcf.gz $i | \
bcftools norm -m -both -f /share/data4/TCGA/Germline/GRCh37-lite/references_GRCh37lite_GRCh37-lite.fa.gz \
-Ou -|bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' |sed 's/\.\/\././g;s/0\/0/0/g;s/0\/1/1/g;s/1\/1/2/g' > $i.geno
} &
done

for i in {1..22} X Y;
do
{
awk '{OFS="\t";print $1,$2,".",$3,$4,".",".","."}' $i.geno > $i.vcf
}&
done
```

## annotation variants

```sh
/share/apps/annovar_2018Apr16/annotate_variation.pl -downdb -buildver hg19 -webfrom annovar esp6500siv2_all humandb/
for i in `seq 1 22` X Y;
do
{/share/apps/annovar_2018Apr16/table_annovar.pl /share/data4/TCGA/Germline/split_vcf/$i.vcf \
    -vcfinput /share/data0/reference/annovar_humandb/hg19 \
    -build hg19 -out $i -otherinfo -remove -nastring . \
    -protocol refGene,avsnp150,1000g2015aug_eur,gnomad211_genome,exac03 -operation g,f,f,f,f
}&
done

# zcat *.hg19_multianno.txt.gz|wc
# 21795172 1048243494 6147601213
```

## filter non-exonic variants

```sh

for i in {10..19};
do
julia get_exonic.jl $i exonic
done

```

## get the variant of CEU samples

```R
load("~/Projects/TCGA/hg19/workflow/01maf/mem.RData")

chrom = 22
library(stringr)
library(parallel)
get_pop_variant = function(chrom){
    anno <- read.table(gzfile(paste0(chrom, ".exonic.anno.gz")), header=T,sep="\t",stringsAsFactors=F)
    geno <- read.table(gzfile(paste0(chrom, ".exonic.geno.gz")), header=T,sep="\t",stringsAsFactors=F, check.names = F)
    colnames(geno) = substr(colnames(geno), 1, 12)
    pop = c("CEU", "CHB", "ASW")
    geno = lapply(pop, function(x){
        if(x == "CEU"){
            AF = "ExAC_NFE"
        }else if(x == "CHB"){
            AF = "ExAC_EAS"
        }else{
            AF = "ExAC_AFR"
        }
        geno = data.matrix(geno[,intersect(mem$sample[mem$pop == x], colnames(geno))])
        index = which(rowSums(is.na(geno)) < ncol(geno))
        geno = geno[index, ]
        anno = anno[index, c("Chr", "Start", "End", "Ref", "Alt", "Gene_refGene", "ExonicFunc_refGene", "AAChange_refGene",
                    "avsnp150", AF)]
        return(list(geno = geno, anno = anno))
    })
    names(geno) = pop
    return(geno)
}

all = mclapply(1:22, function(chrom){
    print(chrom)
    return(get_pop_variant(chrom))
}, mc.cores = 22)

library(data.table)

geno_CEU = rbindlist(lapply(1:22, function(x){
    return(as.data.frame(all[[x]]$CEU$geno))
}))
anno_CEU = rbindlist(lapply(1:22, function(x){
    return(all[[x]]$CEU$anno)
}))

library(stringr)
library(parallel)
library(reshape2)

filter_dup_variant = function(geno, anno){
    gc()
    geno[is.na(geno)] = 0
    geno = as.data.frame(geno)
    ID = paste(anno$Chr, anno$Start, anno$End, anno$Ref, anno$Alt, sep = ":")
    dupID = unique(ID[duplicated(ID)])
    dup_geno = mclapply(dupID, function(x){
        # print(grep(x, dupID))
        g = geno[which(ID == x), ]
        g = apply(g, 2, max)
        return(g)
    }, mc.cores=45)
    dup_geno = do.call(rbind, dup_geno)
    rownames(dup_geno) = dupID
    geno = geno[-which(ID %in% dupID), ]
    rownames(geno) = ID[-which(ID %in% dupID)]
    geno = rbind(geno, dup_geno)
    anno = unique(cbind(ID, anno))
    geno = geno[anno$ID, ]
    anno$type = ifelse(anno$ExonicFunc_refGene == "nonsynonymous SNV", "miss", "trun")
    # Mutational Signatures in Breast Cancer: The Problem at the DNA Level
    anno$type[which(str_length(anno$Alt) > 50 | str_length(anno$Ref) > 50)] = "structural"
    anno[, c("Chr","Start","End","Ref","Alt")] = NULL
    genes = data.frame(ID = anno$ID, str_split(anno$Gene_refGene, pattern=";", simplify=T), stringsAsFactors = F)
    genes = melt(genes, id.vars = "ID")
    genes = genes[which(genes$value != ""), ]
    genes = unique(genes[, -2])
    anno = anno[match(genes$ID, anno$ID), ]
    anno$Gene_refGene = genes$value
    anno = anno[!grepl("^OR",anno$Gene_refGene), ]
    anno = anno[!grepl("^HLA",anno$Gene_refGene), ]
    return(list(geno = geno, anno = anno))
}

CEU = filter_dup_variant(geno_CEU, anno_CEU)

save(CEU, file = "CEU.RData")
```

## get genes' genetic burdens of missense, truncation, structural

```R
load("CEU.RData")

library(tidyverse)
library(parallel)
library(stringr)

anno = CEU$anno
geno = CEU$geno

anno$MAF = as.numeric(anno[[grep("ExAC", colnames(anno))]])

anno$flipped = ifelse(anno$MAF > 0.5, 1, 0)

anno$MAF[which(anno$MAF > 0.5)] = 1 - anno$MAF[which(anno$MAF > 0.5)]

# dim(anno)
# 557872      9
# length(which(anno$MAF > 0.05))
# 19796

anno = anno[which(is.na(anno$MAF) | (anno$MAF < 0.05 & anno$MAF > 0)), ]
anno$MAF[is.na(anno$MAF)] = min(anno$MAF, na.rm = T)

# dim(anno)
# 524912      9

geno = geno[unique(anno$ID), ]
geno[intersect(rownames(geno), anno$ID[anno$flipped == 1]),] = 2 - geno[intersect(rownames(geno), anno$ID[anno$flipped == 1]),]

bed = as.data.frame(str_split(anno$ID, pattern=":",simplify=T), stringsAsFactors=F)
colnames(bed) = c("Chr", "Start", "End", "Ref", "Alt")
bed$ref_len = str_length(bed$Ref)
bed$alt_len = str_length(bed$Alt)
bed$type = ifelse(anno$ExonicFunc_refGene == "nonsynonymous SNV", "miss", "trun")
bed$type[which(bed$alt_len > 50 | bed$ref_len > 50)] = "structural"
anno$type = bed$type

anno = anno[, c("ID", "Gene_refGene", "type", "MAF", "flipped")]
anno = anno %>% group_by(Gene_refGene, type) %>% mutate(weight = 1/sqrt(MAF *(1 - MAF))) %>% as.data.frame
geno = t(geno)
genes = unique(anno$Gene_refGene)
burden = mclapply(genes, function(x){
    a = anno[which(anno$Gene_refGene == x),] %>% group_by(type) %>% summarise(burden = list(data.matrix(geno[, ID]) %*% weight))
    burden = do.call(cbind, a$burden)
    colnames(burden) = a$type
    rownames(burden) = rownames(geno)
    return(burden)
}, mc.cores = 10)
names(burden) = genes
save(geno, anno, burden, file = "burden.RData")
```
