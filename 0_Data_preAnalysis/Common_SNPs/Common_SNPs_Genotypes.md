# Filter sample with genotypes

```R
#!/usr/bin/env Rscript
# @Date    : 2018/11/8 11:21:18
# @Author  : JT Guo
# @Email   : guojt-4451@163.com


# Args <- commandArgs(T)

library(stringr)
files = list.files("TCGA/SNP6_Genotype/00_gtypRData", pattern="*RData",
    recursive=T, full.names=T)

read_gt <- function(file){
    cancer = tolower(str_split(basename(file), pattern="[_|.]", simplify=T)[2])
    load(file)
    gc()
    return(get(paste0("gtyp.normal.", cancer)))
}

geno = mclapply(files, read_gt, mc.cores=14)
geno = do.call(cbind, geno)

geno = read_gt(files[1])
for(i in files[-1]){
    geno = cbind(geno, read_gt(i))
    gc()
}

load("TCGA/SNP6_Genotype/map_SNP6_1KG_dbSNP138.RData")

int = intersect(rownames(geno), SNPs$Affy_Probe_ID)
geno <- geno[int, ]
SNPs <- SNPs[int, c("Affy_Probe_ID", "TGP_RS_ID", "Chromsome_hg19", "Physical_Position_hg19", 
    "Allele_A_1KG", "Allele_B_1KG", "flipped")]


# remove duplicated colnums
geno = geno[, !duplicated(colnames(geno))]

# short-names
Triname = substr(colnames(geno), 1, 12)

# duplicated short-names, exp. TCGA-E2-A1IE
dupTriname = unique(Triname[duplicated(Triname)], 1,12)

# unique samples
name = colnames(geno)[-which(Triname %in% dupTriname)]

# duplicated full-names
dupname = sapply(dupTriname, function(x){colnames(geno)[grep(x, colnames(geno))]})
# selected the DNA samples
dupname = sapply(dupname, function(x){
    return(x[which(substr(x,14,14) == 1 & substr(x,20,20) == "D")])
    })
# selected "Blood Derived Normal"
dupname = sapply(dupname, function(x){
    if(length(x) > 1){
        y = x[which(substr(x, 14,15) == 10)]
        if(length(y) < 1){
            return(x)
        }else{
            return(y)
        }
    }else{
        return(x)
    }
})
# selected the higher plate, https://confluence.broadinstitute.org/pages/viewpage.action?pageId=67404946
dupname = sapply(dupname, function(x){
    if(length(x) > 1){
        return(sort(x)[length(x)])
    }else{
        return(x)
    }
})

write.table(c(name, dupname),"Sample.5676.txt", col.names=F, row.names=F, quote=F)

geno = geno[, c(name, dupname)]

for(i in 1:nrow(geno)){
    geno[i, ][which(geno[i, ] == -1)] = NA
    if(SNPs$flipped[i]){
        geno[i, ] = 2 - geno[i, ]
    }
}

### remove dup rsid

save.image("raw.flipped.geno.RData")
```

# population classification

```R

load("~/TCGA/Genotyping/TGP/906600_all/TGP.SNP6.geno.RData")
load("~/TCGA/Genotyping/01RawSNP6/raw.flipped.geno.RData")

SNPs$ID = paste(gsub("chr","",SNPs$Chromsome_hg19), 
  SNPs$Physical_Position_hg19, SNPs$Allele_A_1KG, SNPs$Allele_B_1KG, sep=":")

int <- intersect(tgp$ID, SNPs$ID)

geno = geno[match(int, SNPs$ID), ]
SNPs = SNPs[match(int, SNPs$ID), ]
tgp = tgp[match(int, tgp$ID), ]
rownames(tgp) = SNPs$TGP_RS_ID
SNPs$INFO = tgp$INFO
tgp[, 1:3] = NULL

rsNA = apply(geno, 2, function(x){length(which(is.na(x)))})
sampleNA = apply(geno, 1, function(x){length(which(is.na(x)))})

tgp1 <- matrix(as.numeric(as.matrix(tgp)), nrow = nrow(tgp))
rownames(tgp1) = rownames(tgp)
colnames(tgp1) = colnames(tgp)
tgp = tgp1
rm(tgp1)
gc()

panel <- read.table("~/Projects/TCGA/Genotyping/TGP/TGP.sample",
  stringsAsFactors=F, sep="\t")
colnames(panel) = c("sample", "gender", "ID", "pop1", "pop1_name", "pop2", "pop2_name", "phase")
panel_tgp <- panel[which(panel$pop1 %in% c("CHB","CEU","ASW","YRI","MXL")), c("sample", "pop1")]

SNPs$sample_AF = rowMeans(geno)/2

# save.image("pop.all.RData")

################# PCA
panel_tgp <- panel[which(panel$pop1 %in% c("CHB","CEU","ASW","YRI","MXL")), c("sample", "pop1")]
tgp1 <- tgp[, panel_tgp[,1]]
geno1 <- cbind(tgp1, geno)


pre_rs = as.character(SNPs$TGP_RS_ID[which(SNPs$sample_AF > 0.49 & SNPs$sample_AF < 0.51)])

library(GeneGeneInteR)
new_geno <- new("SnpMatrix", as.matrix(t(geno[pre_rs, ])) + 1)
new_geno1 <- snpMatrixScour(new_geno, min.maf = 0.49, min.eq = 0.05, call.rate = 0.95)
rs = colnames(as.matrix(new_geno1$snpX))

geno_com <- geno1[which(rownames(geno1) %in% rs), ]
geno_com1 <- matrix(as.numeric(as.matrix(geno_com)), nrow = nrow(geno_com))
rownames(geno_com1) = rownames(geno_com)
colnames(geno_com1) = colnames(geno_com)

print("PCA...")
res <- prcomp(t(geno_com1), center = TRUE, scale. = TRUE)
pop <- panel_tgp[match(colnames(geno_com1), panel_tgp[,1]), 2]
pop[is.na(pop)] <- "Case"

library(ggplot2)
library(ggbiplot)
library(RColorBrewer)

pdf("Figure_S1A.pdf")
plot(res, type = "l")
res_plot = ggbiplot(res, choices = 1:2, obs.scale = 1, var.scale = 1, groups = pop, 
  ellipse = TRUE, circle = TRUE, var.axes = F) +
scale_color_manual(name=pop, values = brewer.pal(8, "Set1"))
plot(res_plot)
pca = as.data.frame(predict(res)[,1:2])
pca$pop = panel_tgp[, 2][match(rownames(pca), panel_tgp[, 1])]
pca$pop[is.na(pca$pop)] = "Case"
pca_plot = ggplot() + 
  geom_point(data=pca[which(pca$pop=="Case"), ], 
    aes(x = PC1, y = PC2), size = 0.5, colour = "gray", alpha = 1) +
  geom_point(data = pca[which(pca$pop!="Case"), ],  
    aes(x = PC1, y = PC2, colour = pop), size = 1, alpha = 0.8) +
  scale_color_brewer(palette = "Set1") +
  theme(legend.title = element_blank()) +
  coord_equal() +
  geom_rect(aes(xmin=-12, xmax=10, ymin=27, ymax=40), color="gray70", alpha=0.1) +
  geom_text(aes(x=-5, y=35, label="CHB"), size=4) +
  geom_rect(aes(xmin=-12, xmax=5, ymin=2, ymax=27), color="gray70", alpha=0.1) +
  geom_text(aes(x=-5, y=15, label="MXL"), size=4) +
  geom_rect(aes(xmin=-12, xmax=2, ymin=-10, ymax=2), color="gray70", alpha=0.1) +
  geom_text(aes(x=-5, y=0, label="CEU"), size=4) +
  geom_rect(aes(xmin=15, xmax=38, ymin=-12, ymax=2), color="gray70", alpha=0.1) +
  geom_text(aes(x=30, y=0, label="ASW"), size=4) +
  geom_rect(aes(xmin=38, xmax=45, ymin=-12, ymax=2), color="gray70", alpha=0.1) +
  geom_text(aes(x=41, y=0, label="YRI"), size=4) +
  theme_classic()
plot(pca_plot)
dev.off()

save(geno_com1, res, pca, res_plot, pca_plot, file = "PCA.res.RData")

write.table(pca, "PCA.txt", col.names=T, row.names=T, sep="\t", quote=F)
```
