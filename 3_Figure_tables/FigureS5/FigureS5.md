# Figure S5

```R
library(stringr)
library(survival)
library("survminer")
library(gridExtra)
library(cowplot)

plotSurv <- function(gene, snv, clin, cancer="panCan", plot = T, plot_table = T, legend.title = NULL, legend.labs = NULL, title = NULL,
    xlab = "time", ylab = "Survival probability"){
        dat = cbind(clin, snv = snv[clin$PATIENT_ID, gene])
        if(cancer=="panCan"){
            fit1 = try(coxph(Surv(OS_MONTHS, status) ~ snv + DRUG_TYPE + CANCER_TYPE, data = dat))
        }else{
            dat = dat[dat$CANCER_TYPE == cancer, ]
            fit1 = try(coxph(Surv(OS_MONTHS, status) ~ snv + DRUG_TYPE, data = dat))
        }
        if(class(fit1) == "try-error"){
            return()
        }
        res = summary(fit1)
        hr = cbind(res$conf.int, res$coef)[, c(1,3,4,9)]
        if(!is.null(names(hr))){
            names(hr) = c("HR", "HR.CI.lower", "HR.CI.upper", "pvalue")
            hr = t(as.data.frame(hr))
            rownames(hr) = NULL
        }else{
            colnames(hr) = c("HR", "HR.CI.lower", "HR.CI.upper", "pvalue")
        }
        hr = as.data.frame(signif(hr, 3))
        hr = hr[grep("snv", rownames(hr)), ]
        if(!plot){
            hr$gene = gene
            hr$cancer = cancer
            return(hr)
        }
        fit2 = surv_fit(as.formula("Surv(OS_MONTHS, status) ~ snv"), data = dat)
        p = ggsurvplot(fit2, pval = F, risk.table = TRUE,
            font.title = c(20, "plain", "bold"), font.x = c(20, "plain", "black"), font.y = c(20, "plain", "black"),
            font.tickslab = c(15, "plain", "black"), font.legend = c(15, "plain", "black"),
            palette = "Dark2", legend.title = legend.title, legend.labs = legend.labs, ggtheme = theme_minimal())
        p1 = p$plot + labs(title = title, x = xlab, y = ylab) +
            theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position=c(0.8, 0.9), legend.direction="vertical", aspect.ratio = 1) +
            annotation_custom(tableGrob(hr, theme = ttheme_minimal(base_size = 15)),
            xmin = min(p$data.survplot$time), xmax = max(p$data.survplot$time), ymin = 0, ymax = 0.2)
        if(plot_table){
            p2 = p$table + theme(legend.position="none", axis.text.y = element_blank()) + labs(x = xlab)
            p = plot_grid(p1, p2, rel_heights=c(3,1), ncol = 1, nrow=2, align = 'v')
            return(p)
        }else{
            return(p1)
        }
}


I_genes <- read.table("../TableS2/Table_S2C.csv", header=T,sep="\t",stringsAsFactors = F)

# Tumor_mutational_load_predicts_survival_after_immunotherapy_across_multiple_cancer_types
clin <- read.table("data_clinical_patient.txt", header=T,sep="\t",stringsAsFactors = F)
samples <- read.table("data_clinical_sample.txt", header=T,sep="\t",stringsAsFactors = F)
clin$status = as.numeric(str_split(clin$OS_STATUS, pattern = ":", simplify = T)[,1])
clin$CANCER_TYPE = samples$CANCER_TYPE
maf <- read.table("data_mutations_mskcc.txt", header=T,sep="\t",stringsAsFactors = F)
snv = as.data.frame.matrix(table(maf[-which(maf$Variant_Classification %in% c("Intron", "5'Flank", "3'Flank", "5'UTR", "Silent")), c("Tumor_Sample_Barcode", "Hugo_Symbol")]))
snv[snv>0] = 1
rownames(snv) = paste0("P-", str_split(rownames(snv), pattern = "-", simplify = T)[,2])

snv = snv[,intersect(I_genes$gene, colnames(snv))]

snv_surv = do.call(rbind, lapply(colnames(snv), function(x){
    return(plotSurv(x, snv, clin, plot = F))
}))

p_list_snv_surv = lapply(snv_surv$gene[which(snv_surv$pvalue<0.05)], function(x){
    return(plotSurv(x, snv, clin, plot = T, legend.title = x, xlab = "Month"))
})

pdf("FigureS5.pdf", width = 12, height = 18)
plot_grid(plotlist = p_list_snv_surv, ncol = 2, align = "hv", axis = "tbrl")
dev.off()
```
