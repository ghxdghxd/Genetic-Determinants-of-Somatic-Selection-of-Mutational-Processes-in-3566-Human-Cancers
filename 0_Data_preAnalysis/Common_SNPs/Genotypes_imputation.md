# Genotype imputation

## get raw genotype of each population

```R
library(dplyr)
load("raw.flipped.geno.RData")
geno[SNP_With_NA, sample_With_NA] = -1
geno = as.data.frame(geno)
cancer <- read.table('cancer.txt', header = T, sep = "\t", stringsAsFactors = F)
cancer = unique(cancer[which(cancer$samplename %in% colnames(geno)),])

pop = "CEU"
library(parallel)
get_pop_geno = function(pop){
    POP <- read.table(paste0("TCGA/SNP6_Genotype/samples/", pop, ".sample"),header=T,sep=" ",stringsAsFactors=F)
    cancer = cancer[which(cancer$samplename %in% POP$ID_1),]
    res = mclapply(unique(cancer$cancer), function(c){
        assign(c, cbind(SNPs[,1:6], geno[,cancer$samplename[which(cancer$cancer == c)]]))
        write.table(get(c), gzfile(paste0(c, ".geno.txt.gz")), r=F,c=T,sep="\t",quote=F)
        return(c)
    }, mc.cores=15)
}
```

## raw genotype To shapeit

> TCGA/SNP6_Genotype/0_raw_genotypes/CEU
> +++ To === >
> TCGA/SNP6_Genotype/1_for_shapeit/CEU/\*/\*.geno.txt.gz

```sh
python ToShapeit.py
```

## get genotype of each chromsome

>TCGA/SNP6_Genotype/1_for_shapeit/CEU

```sh
for j in `ls -1|grep -v ESCA`;
do
    for i in `seq 1 22` X;
    do
        zgrep -w chr$i $j/$j.geno.txt.gz| sed '1d;s/,/\t/g'|cut -f 1-2,4-|gzip > $j/$j.CEU.chr$i.geno.gz
    done
done
```

## shapeit + impute2

>TCGA/SNP6_Genotype/2_shapeit_phased
>TCGA/SNP6_Genotype/3_imputed_geno

```sh
#!/bin/bash
chrom=$1
cancer=$2

python GenoImpute.py phase --qsub -p 8 -g TCGA/SNP6_Genotype/1_for_shapeit/CEU/$cancer/$cancer.CEU.chr$chrom.geno.gz \
-s TCGA/SNP6_Genotype/1_for_shapeit/CEU/$cancer/$cancer.CEU.sample \
-c chr$chrom \
-o impute2/TCGA/$cancer


start=(`grep -w ^$chrom impute.region.txt|cut -f 2`)
end=(`grep -w ^$chrom impute.region.txt|cut -f 3`)

for((i=0;i < ${#start[@]};i++));
do
    python GenoImpute.py impute -g TCGA/SNP6_Genotype/2_shapeit_phased/$cancer/$cancer\_CEU_chr$chrom\_phased.haps.gz \
        -c chr$chrom -s ${start[$i]} -e ${end[$i]} \
        -o impute2/TCGA/$cancer/chr$chrom
done
```

## fcgene

> 5_fcGene_geno_thresh_***
> TCGA/SNP6_Genotype/5_fcGene_geno_thresh_maxProb0.9_info0.9

```sh
for y in BRCA COAD ESCA GBM KIRC LIHC LUAD OV PRAD STAD THCA UCEC;
do {
    mkdir $y
    for i in TCGA/SNP6_Genotype/3_imputed_geno/$y/*/*.gz; do
    name=`basename ${i%.*}`
    fcgene --gens $i \
    --thresh 0.9 \
    --info ${i%.*}\_info \
    --info-thresh 0.9 \
    --oformat r \
    --transpose \
    --force ref-allele=allele2 \
    --out TCGA/SNP6_Genotype/5_fcGene_geno_thresh_maxProb0.9_info0.9/$y/$name >/dev/null 2>&1
    done
} & done
```

## filter CN variants from maxProb0.9_info0.9

> 6_filterCN_INS_fcGene_geno_thresh_maxProb0.9_info0.9

```sh
for y in BRCA COAD ESCA GBM KIRC LIHC LUAD OV PRAD STAD THCA UCEC;
do
{
    for i in {1..22};
    do
        sample=`cut -f 1 TCGA/SNP6_Genotype/2_shapeit_phased/$y/$y\_CEU_chr10_phased.sample| \
            sed '1,2d;3irsid'|cut -d '-' -f 1-3| paste -s -d '\t'`
        cat TCGA/SNP6_Genotype/5_fcGene_geno_thresh_maxProb0.9_info0.9/$y/$y\_CEU_chr$i\_*_genotype.txt| \
        grep -v -e "^rsid" -e "<" | sed -e 's/ \+/\t/g' -e 1i$sample | \
        bgzip -@ 3 -c > TCGA/SNP6_Genotype/6_filterCN_INS_fcGene_geno_thresh_maxProb0.9_info0.9/$y\_CEU_chr$i.geno.gz
    done
} &
done

# BRCA 514492
# COAD 510835
# ESCA 524669
# GBM 509795
# KIRC 511514
# LIHC 517081
# LUAD 511055
# OV 511850
# PRAD 516644
# STAD 522380
# THCA 510224
# UCEC 512643
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 509795  511000  512246  514432  516753  524669
```

# maxProb0.9_info0.9 To Rdata

## get all imputed loci, filter CN variants

```sh
# TCGA/hg19/workflow/03genotype/02imputeALL/new_info0.9_maxPro0.9/0ID
# 80112006
for i in `seq 1 22`;
do
{
    cut -d ' ' -f 2-5 TCGA/SNP6_Genotype/3_imputed_geno/ESCA/chr$i/*_info| \
    grep -v -e "rs_id" -e "<"|awk -F " " -v chr="$i" '{OFS=":";print chr,$2,$3,$4"\t"$1}' > $i.info.txt
} &
done
```

## get CEU AF from TGP3

> TCGA/hg19/workflow/03genotype/02imputeALL/new_info0.9_maxPro0.9/1AF

```R
library(parallel)
library(data.table)
get_AF = function(chr){
    tgp_file = paste0("/share/data0/reference/ALL_1000G_phase3integrated_v3_impute/1000GP_Phase3_chr", chr, ".legend.gz")
    info_file = paste0("~/Projects/TCGA/hg19/workflow/03genotype/02imputeALL/new_info0.9_maxPro0.9/0ID/", chr, ".info.txt.gz")
    tgp <- fread(tgp_file, header = T, stringsAsFactors = F, sep = " ")
    tgp$ID = paste(chr, tgp$position, tgp$a0, tgp$a1, sep = ":")
    info <- fread(info_file, header = F, stringsAsFactors = F, sep = "\t")
    info$CEU_MAF = tgp$EUR[match(info$V1, tgp$ID)]
    info$CEU_MAF[which(info$CEU_MAF > 0.5)] = 1 - info$CEU_MAF[which(info$CEU_MAF > 0.5)]
    colnames(info)[1:2] = c("ID1", "ID2")
    write.table(info, gzfile(paste0(chr, ".AF.txt.gz")), r = F, c = T, sep = "\t", quote = F)
}

mclapply(1:21, get_AF, mc.cores = 22)
```

## classify the SNPs into common or rare

>TCGA/hg19/workflow/03genotype/02imputeALL/new_info0.9_maxPro0.9/2CEU_com

```sh
for i in `seq 1 22`;
do
{zcat $i.AF.txt.gz | awk '{if($3>0.05||NR==1){print}}' |bgzip > $i.CEU.COM.AF.gz}&
done
# 6,868,961
```

>TCGA/hg19/workflow/03genotype/02imputeALL/new_info0.9_maxPro0.9/3CEU_rare

```sh
for i in `seq 1 22`;
do
{zcat $i.AF.txt.gz | awk '{if(($3>0 && $3<0.05)||NR==1){print}}' |bgzip > $i.CEU.Rare.AF.gz}&
done
#16,783,851
```

## get SNPs genotypes

```sh
for i in BRCA COAD GBM KIRC LIHC LUAD OV PRAD STAD THCA UCEC; do python getGeno.py $i;done
```

```R
library(parallel)

info_files = list.files("TCGA/hg19/workflow/03genotype/02imputeALL/new_info0.9_maxPro0.9/2CEU_com", pattern="*.COM.AF.gz", full.names=T)
info = mclapply(info_files, function(x){a = read.table(x, header=T,sep="\t",stringsAsFactors=F);return(a)}, mc.cores=22)
info = do.call(rbind, info)
geno_files = list.files("TCGA/SNP6_Genotype/6_filterCN_INS_fcGene_geno_thresh_maxProb0.9_info0.9", pattern="*geno.gz", full.names = T)

geno = mclapply(geno_files, function(x){
    geno = read.table(gzfile(x), header=T,sep="\t",stringsAsFactors=F,check.names=F)
    int = intersect(geno$rsid, info$ID2)
    geno1 = geno[match(int, geno$rsid), ]
    geno1$rsid = info$ID1[match(int, info$ID2)]
    write.table(geno1, gzfile(basename(x)), c=T,r=F,sep="\t", quote=F)
    return(geno1)
}, mc.cores = 12)

int = unique(unlist(lapply(geno, function(x){return(x$rsid)})))

geno = mclapply(geno, function(x){
    return(x[match(int, x$rsid), -1])
}, mc.cores = 12)

geno = do.call(cbind, geno)
rownames(geno) = int
write.table(geno, gzfile("all.geno.gz"), c=T,r=T,sep="\t", quote=F)
save(geno, file="all.geno.RData")
```

## annotation

> TCGA/hg19/workflow/03genotype/02imputeALL/new_info0.9_maxPro0.9/1Annovar

```sh
for i in `seq 1 22`;
do
{
    zcat ~/Projects/TCGA/hg19/workflow/03genotype/02imputeALL/new_info0.9_maxPro0.9/0ID/$i.info.txt.gz | \
        awk '{split($1,a,":");OFS="\t";print a[1],a[2],$2,a[3],a[4],".",".","."}' > $i.info.vcf
}&
done

for i in `seq 1 22`;
do
{/share/apps/annovar_2018Apr16/table_annovar.pl $i.info.vcf \
    -vcfinput /share/data0/reference/annovar_humandb/hg19 \
    -build hg19 -out $i.info -otherinfo -remove -nastring . \
    -protocol refGene,avsnp150,1000g2015aug_eur,gnomad211_genome,exac03 -operation g,f,f,f,f
}&
done
```

## get CEU AF from TGP3

```bash
zcat /share/data1/TGP/TGP.AF.txt.gz|cut -f 1-5,9 |awk '{print $1":"$2":"$4":"$5"\t"$6}' > TGP.AF.CEU.txt
gzip TGP.AF.CEU.txt
```

## merge ID, cytoband, EUR AF

> TCGA/hg19/workflow/03genotype/02imputeALL/new_info0.9_maxPro0.9/1AF

```r
library(parallel)
info_files = list.files("", pattern="*.txt$",full.names=T)
AF = read.table(gzfile("TGP.AF.CEU.txt.gz"), header=T, sep="\t", stringsAsFactors = F)
info = mclapply(info_files, function(x){
    a = read.table(x, header = F, stringsAsFactors = F, sep="\t")
    colnames(a) = c("chr","pos","ref","alt","cytoband")
    a$ID = paste(a$chr, a$pos, a$ref, a$alt, sep = ":")
    a$EUR = AF[match(a$ID, AF$chr.pos.ref.alt), "EUR"]
    return(a)
    }, mc.cores = 22)

info = mclapply(1:length(info), function(x){
    a = info[[x]]
    a = a[, c("ID", "cytoband", "EUR")]
    write.table(a, gzfile(paste0(x, ".AF.gz")), c=T,r=F,sep="\t",quote=F)
    }, mc.cores = 22)
```
