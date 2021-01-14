# mpQTLs analysis

## prepare input data

```sh
for c in panCancer BRCA COAD ESCA ESCC GBM KIRC LIHC LUAD OV PRAD STAD THCA UCEC;
do
    for i in {1..22};
    do
        Rscript prepare_input_data.R $c $i
    done
done
```

## mpQTLs association

```sh

for c in panCan BRCA COAD ESCA ESCC GBM KIRC LIHC LUAD OV PRAD STAD THCA UCEC;
do
    for i in APOBEC EPP1 EPP2 MMR Tobacco NpCpG;
    do
    {
        for j in {22..1};
        do
            julia run_lm.jl $j $i $c mat_out
        done
    } &
    done
done
```

## select mpQTLs

```R
library(dplyr)
library(parallel)
library(ggplot2)

mat_files = list.files(".", pattern = "*.txt.gz", recursive=T, full.names=T)

mat = mclapply(mat_files, function(x){
    print(x)
    mat <- read.table(gzfile(x), header=T,sep="\t",stringsAsFactors=F)
    return(mat)
}, mc.cores = 1)
mat = do.call(rbind, mat)
mat_sig = mat[which(mat$P < 5e-8), ]

ggQQplot <- function(pvalues, title = NULL, random = NULL){
    pvalues = sort(na.omit(as.numeric(pvalues)))
    pvalues = cbind(op = pvalues, ep = ppoints(length(pvalues)))
    if(!is.null(random)){
        pvalues = pvalues[sample(1:nrow(pvalues), random), ]
    }
    lab <- quantile(pvalues[,"op"], 0.475)/quantile(pvalues[,"ep"], 0.475)
    pvalues <- -log10(pvalues)
    maxP = ceiling(max(pvalues))
    if(is.null(title)){
        title = ""
    }
    p = ggplot(as.data.frame(pvalues)) +
        geom_point(aes(ep, op), size = 1) +
        scale_x_continuous(expand = c(0, 0), limits = c(0, maxP)) +
        scale_y_continuous(expand = c(0, 0), limits = c(0, maxP)) +
        geom_abline(intercept = 0, slope = 1, alpha = 1) +
        labs(x = "Expected (-log10(P))", y = "Observed (-log10(P))", title = title) +
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5, size = 15),
            axis.text = element_text(hjust = 0, size = 10),
             plot.margin=unit(c(1,2,1,1),"lines")) +
        annotate("text", x = 0.7 * maxP, y = 0.2 * maxP, size = 5,
            label = paste0("lambda = ", signif(lab, 3)))
    return(p)
}

p_list = lapply(c("BRCA","COAD","ESCA","ESCC","GBM","KIRC","LIHC","LUAD","OV","PRAD","STAD","THCA","UCEC","panCan"), function(x){
    print(x)
    p1 = lapply(c("APOBEC", "EPP1", "EPP2", "MMR", "NpCpG", "Tobacco"), function(y){
        print(y)
        return(ggQQplot(mat$P[mat$cancer == x & mat$sig == y], title = paste0(x, '-', y)) )
    })
    return(cowplot::plot_grid(plotlist = p1, nrow=1))
})

png("QQplot.png", width = 1500, height = 3000)
gridExtra::marrangeGrob(p_list, nrow = 14, ncol = 1, top = NULL, padding = unit(2, "line"),
    layout_matrix = matrix(seq_len(14 * 1), nrow = 14, ncol = 1, byrow = T))
dev.off()
```
