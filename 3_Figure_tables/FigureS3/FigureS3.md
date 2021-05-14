# FigureS3

```R
library(pmsignature)

table_pca = "../FigureS1/Table_PCA.csv"

pca <- read.table(table_pca, header=T, sep="\t", stringsAsFactors = F)
pca$pop = pca$population
pca$pop[which(pca$source=='TCGA')]="TCGA"


sig_RData="../FigureS2/signature_7.RData"

loadPmSig <- function(pmSigRData){
    load(pmSigRData)
    return(list(G = G, BG_prob = BG_prob, Param = Param))
}

res_pmSig = loadPmSig(sig_RData)

## get the signature activity in each sample
pmSigMat = getMembershipValue(res_pmSig$Param)
colnames(pmSigMat) = c("EPP1","Tobacco", "APOBEC", "NpCpG", "dMMR", "EPP2", "BG")
pmSigMat = pmSigMat[, c("NpCpG", "APOBEC", "dMMR", "Tobacco", "EPP1", "EPP2", "BG")]
rownames(pmSigMat) = substr(rownames(pmSigMat), 1, 12)
colnames(pmSigMat) = paste0("MS", 1:7)

getLogMat <- function(sigMat, BG, threshold = 1e-3, trans = "log") {
    sigMat <- as.matrix(sigMat)
    sigMat[sigMat < threshold] <- 0
    sigMat <- sigMat[which(sigMat[, BG] > 0), ]
    sigMat[which(sigMat == 0)] <- NA
    if (trans == "log") {
        sigMat <- log(sigMat / sigMat[, BG])
    } else if (trans == "log2") {
        sigMat <- log2(sigMat / sigMat[, BG])
    }
    sigMat <- sigMat[, -grep(BG, colnames(sigMat))]
    sigMat <- sigMat[which(rowSums(sigMat, na.rm = T) != 0), which(colSums(sigMat, na.rm = T) != 0)]
    return(as.data.frame(sigMat))
}
pmLogSigMat = getLogMat(pmSigMat, "MS7", threshold = 1e-3, trans = 'log')
colnames(pmLogSigMat) = paste0("MP", 1:6)

pmLogSigMat$cancer = pca$cancer_type[match(rownames(pmLogSigMat), pca$sampleID)]

subCloneNum = read.table("journal.pgen.1007669.csv", header=T,sep="\t",stringsAsFactors = F)
subCloneNum$sample_name = substr(subCloneNum$sample_name, 1, 12)

pmLogSigMat$tree = subCloneNum$Tree.score[match(rownames(pmLogSigMat), subCloneNum$sample_name)]
pmLogSigMat$clone_num = subCloneNum$number.of.clones[match(rownames(pmLogSigMat), subCloneNum$sample_name)]

library(vegan)
mat_shannon = exp(pmLogSigMat[,paste0("MP", 1:6)])
mat_shannon[is.na(mat_shannon)] = 0
shannon_index = diversity(mat_shannon, index="shannon")
pmLogSigMat$shannon = shannon_index[rownames(pmLogSigMat)]
pmLogSigMat$clone_num_status = NA
pmLogSigMat$clone_num_status[pmLogSigMat$clone_num > 0] = "low"
pmLogSigMat$clone_num_status[pmLogSigMat$clone_num > 2.5] = "middle"
pmLogSigMat$clone_num_status[pmLogSigMat$clone_num > 3.6] = "high"

pmLogSigMat$cancer = factor(pmLogSigMat$cancer, levels = c("THCA","PRAD","BRCA","KIRC","ESCC","GBM","OV","UCEC","ESCA","LIHC","STAD","COAD","LUAD"))

library(ggsci)

p1 = ggviolin(na.omit(pmLogSigMat[, c("cancer","clone_num")]), x="cancer", y="clone_num", xlab = "", ylab = "Clone counts", fill = "gray") + theme_presentation()
p2 = ggplot(pmLogSigMat, aes(x=clone_num)) +
    geom_histogram(aes(y = ..density..), colour="black", fill="gray", bins = 100, position="identity") +
    geom_density(alpha=.2) +
    geom_vline(xintercept=c(2.5, 3.65), linetype="dashed") + 
    labs(x="Clone counts") +
    theme_presentation()
p3 = ggboxplot(na.omit(pmLogSigMat[, c("clone_num_status", "shannon")]), fill = "clone_num_status",
          outlier.shape = NA, add.params = list(size=0.8), add = c("jitter"), order = c("low", "middle","high"),
          x= "clone_num_status", xlab = "Clone count levels", y ="shannon", ylab = "Shannon's index", 
          width = 0.98, ggtheme = theme_presentation()) + scale_fill_d3() +
    stat_compare_means(comparisons = list(c("low","middle"), c("middle","high"),c("low","high"))) +
    theme(aspect.ratio = 1, legend.position = "none")

pdf("FigureS3.pdf", width = 15, height = 6)
plot_grid(plot_grid(p1, p2, ncol = 1, labels = c("A", "B"), align = 'hv'), p3, nrow = 1, rel_widths = c(3,1.5), labels = c("","C"))
dev.off()
```
