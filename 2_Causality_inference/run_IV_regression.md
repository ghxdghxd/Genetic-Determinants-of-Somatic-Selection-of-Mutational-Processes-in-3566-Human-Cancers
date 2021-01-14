# Causality inference

## candidate genes

### hypothesis 1 (E-gene)

```sh
julia Gene_hypothesis1.jl \
    all.50.flipped.run.RData \
    mat_sig.txt \
    out_dir \
    mat_sig_IV_Egenes
```

### hypothesis 2 (I-gene)

```sh
julia Gene_hypothesis2.jl \
    all.50.flipped.run.RData \
    mat_sig.txt \
    out_dir \
    mat_sig_IV_Igenes
```

## mpQTLs

### get the cis genes of mpQTLs

```sh
for i in `seq 1 22`; do
{
    cut -f 3 mat_sig.txt|grep "^$i:"|sort -u| \
        xargs -i zgrep -w {} ~/Projects/TCGA/hg19/workflow/03genotype/02imputeALL/new_info0.9_maxPro0.9/0ID/$i.info.txt.gz >> ID
} &
done

cut -f 3 mat_sig.txt| xargs -i grep -w {} ID |cut -f 2 > ID.txt
paste -d "\t" sig.mat.txt ID.txt > sig.mat.txt.1
cut -f 3 sig.mat.txt.1 |awk -F ':' 'NR>1{OFS="\t";print $1,$2,$2,$3,$4}'| \
    intersectBed -a - -b /share/data0/reference/Genome/hg19/cytoband.txt -wa -wb|cut -f 9|paste -d "\t" sig.mat.txt.1 - > sig.mat.txt.2
sed "1icancer\tsig\tID\texp\tP\tID2\tcytoband" sig.mat.txt.2 > sig.mat.txt.1
rm ID ID.txt sig.mat.txt.2
```

### hypothesis 1 (E-mpQTLs)

```sh
for j in panCancer BRCA COAD ESCA ESCC GBM KIRC LIHC LUAD OV PRAD STAD THCA UCEC;
do
    for i in {1..22};
    do
        julia mpQTLs_hypothesis1.jl $i $j out_dir_E mat_sig.txt.1
    done
done
```

### hypothesis 2 (I-mpQTLs)

```sh
for j in panCancer BRCA COAD ESCA ESCC GBM KIRC LIHC LUAD OV PRAD STAD THCA UCEC;
do
    for i in {1..22};
    do
        julia mpQTLs_hypothesis2.jl $i $j out_dir_I mat_sig.txt.1
    done
done
```
