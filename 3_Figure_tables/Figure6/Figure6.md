# Figure 6

## Figure 6A

```R
library(tidyverse)
library(circlize)
library(tidyverse)
library(reshape2)
library(ggsci)

mat_sig <- read.table("mpQTLs_1e-5.txt", header = T, sep = "\t", stringsAsFactors = F)
matMH = dcast(mat_sig %>% group_by(snp, sig) %>% summarise(P = min(P)), snp~sig)

matMH$BG = rowSums(matMH[,2:7], na.rm = T)
matMH$BG[which(matMH$BG < 5e-8)] = NA
matMH[which(!is.na(matMH$BG)),2:7] = NA

matMH = cbind(data.frame(str_split(matMH$snp, pattern = ":", simplify = T), stringsAsFactors = F), matMH)
matMH = matMH[, c(1:2,2,6:12)]
colnames(matMH)[1:3] = c("chr","start","end")
# matMH = arrange(matMH, chr, start, end)
matMH$chr = as.numeric(matMH$chr)
matMH$start = as.numeric(matMH$start)
matMH$end = as.numeric(matMH$end)
matMH[,4:10] = -log10(matMH[,4:10])
matMH[is.na(matMH)] = 0
matMH$chr = paste0("chr", matMH$chr)
matMH = matMH[, c("chr","start","end","BG", paste0("MP",1:6))]

pdf("Figure6A.pdf", width = 5, height = 5)
circos.par(gap.after=c(rep(1,21), 10), "start.degree" = 90, track.height = 0.3)
circos.initializeWithIdeogram(chromosome.index = paste0("chr", 1:22), plotType = c("ideogram", "labels"))
circos.genomicTrack(matMH, ylim = c(4.5, 9), panel.fun = function(region, value, ...) {
    col = c("grey60", pal_d3("category10")(10)[1:6])
    cex_list = c(0.3, rep(0.5, 6))
    for(i in 1:ncol(value)){
        value1 = value[which(value[i] != 0), i]
        region1 = region[which(value[i] != 0), ]
        circos.genomicPoints(region1, value1, pch = 16, cex = cex_list[i], col = col[i], transparency = 0.8)
    }
})
circos.yaxis(side = "left", at = seq(4.5, 9, by = 0.5),
             sector.index = get.all.sector.index()[1], labels.cex = 0.5)
circos.clear()
text(-0.1, 0.75, "-log10(pvalue)", cex = 0.6, srt = 95)
legend("center", title = "", pch = 16, cex = 1, bty = "n",
       col = pal_d3("category10")(10)[1:6],
       legend = colnames(matMH)[-c(1:4)])
dev.off()
```

## Figure 6 BCDE

```R
library(rstatix)
IV_sig <- read.table("../TableS4/SupplementaryTable_S4B.csv", header=T,sep="\t", stringsAsFactors = F)
logSigMat <- read.table("../TableS1/SupplementaryTable_S1.csv", header = T, sep= "\t", stringsAsFactors = F)
rownames(logSigMat) = logSigMat$SampleID
load("mRNA_and_geno.RData")

logSigMat = logSigMat[rownames(geno), paste0("MP", 1:6)]

com_plot2 = function(IV_sig, cancer, chrom, nrows){
    p_list = lapply(nrows, function(x){
        snp = IV_sig$chr_pos_ref_alt[x]
        ID = str_split(IV_sig$rsid[x], pattern=":")[[1]][1]
        signature = IV_sig$mutational_propensity[x]
        gene = IV_sig$cis_gene[x]
        a = data.frame(sig = logSigMat[, signature], snp = geno[, snp])
        int = intersect(rownames(geno), rownames(mRNA))
        b = data.frame(snp = geno[int, snp], mRNA = mRNA[int, gene])
        c = data.frame(mRNA = mRNA[int, gene], sig = logSigMat[int, signature])
        pa = ggboxplot(a[!is.na(a$snp), ], x="snp", y="sig", shape = 1, bxp.errorbar= T, outlier.shape = NA, width = 0.5,#add.params = list(size=0.5),
            ylab = signature, xlab = paste0(snp, "\n", ID)) +
            stat_pvalue_manual(a %>% t_test(sig ~ snp) %>% add_xy_position(x = "dose"),
            size = 2, tip.length = 0.01, label = "p.adj")
        pd = ggboxplot(b[!is.na(b$snp), ], x="snp", y="mRNA", shape = 1, bxp.errorbar= T, outlier.shape = NA, width = 0.5,#add.params = list(size=0.5),
            ylab = paste(gene, "expression"), xlab = paste0(snp, "\n", ID)) +
            stat_pvalue_manual(b %>% t_test(mRNA ~ snp) %>% add_xy_position(x = "dose"),
                tip.length = 0.01, label = "p.adj", size = 2, angle = 0)
        p = ggarrange(pa, pd, ncol = 2, align='hv')
        return(p)
    })
    return(p_list)
}

p1 = com_plot2(IV_sig, cancer = "panCan", chrom = "chr20", nrows = 6)
p2 = com_plot2(IV_sig, cancer = "panCan", chrom = "chr6", nrows = 13)

pdf("Figure6BCDE.pdf", width = 6, height = 8)
plot_grid(p1[[1]], p2[[1]], ncol = 1)
dev.off()
```
