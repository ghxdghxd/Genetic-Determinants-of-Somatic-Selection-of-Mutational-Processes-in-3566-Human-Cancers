# FigureS5

```R
load("drug_logSigMat.RData")
mat <- read.table("cell_drug_logSigMat.txt", header=T,sep="\t",stringsAsFactors=F)

data.frame(sig = cell_logSigMat[, "MMR"], drug = cell_drug[, "Dabrafenib"]) %>%
    ggscatter(x = "sig", y = "drug", xlab = "MP3", ylab = "IC50 of BRAF inhibitor(dabrafenib)",
              add = "reg.line",  # Add regressin line
              add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
              conf.int = TRUE, # Add confidence interval
              ggtheme = theme_bw()) +
    annotate("text", x=1,y=1,label="R=")

res = cor.test(cell_logSigMat$MMR, cell_drug$Dabrafenib, method = "pearson")
FDR = signif(mat$FDR[mat$sigs == "MMR" & mat$DRUG_NAME == "Dabrafenib"], 3)
p1_cell = data.frame(sig = cell_logSigMat[, "MMR"], drug = cell_drug[, "Dabrafenib"]) %>%
    ggscatter(x = "sig", y = "drug", xlab = "MP3", ylab = "IC50 of BRAF inhibitor(dabrafenib)",
              add = "reg.line",  # Add regressin line
              add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
              conf.int = TRUE, # Add confidence interval
              ggtheme = theme_bw()) +
    annotate("text", x=2.5,y=-2.5,label = paste0("R = ", signif(res$estimate, 2), "\n", "FDR = ", FDR))

res = cor.test(cell_logSigMat$MMR, cell_drug[,"PLX-4720"], method = "pearson")
FDR = signif(mat$FDR[mat$sigs == "MMR" & mat$DRUG_NAME == "PLX-4720"][2], 3)
p2_cell = data.frame(sig = cell_logSigMat[, "MMR"], drug = cell_drug[, "PLX-4720"]) %>%
    ggscatter(x = "sig", y = "drug", xlab = "MP3", ylab = "IC50 of BRAF inhibitor(PLX-4720)",
        add = "reg.line",  # Add regressin line
        add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
        conf.int = TRUE, # Add confidence interval
        ggtheme = theme_bw()) +
    annotate("text", x = 2.5,y = -1,label = paste0("R = ", signif(res$estimate, 2), "\n", "FDR = ", FDR))

res = cor.test(cell_logSigMat$MMR, cell_drug$Dabrafenib, method = "pearson")
FDR = signif(mat$FDR[mat$sigs == "Tobacco" & mat$DRUG_NAME == "CCT-018159"], 3)
p3_cell = data.frame(sig = cell_logSigMat[, "Tobacco"], drug = cell_drug[, "CCT-018159"]) %>%
    ggscatter(x = "sig", y = "drug", xlab = "MP4", ylab = "IC50 of HSP90 inhibitor(CCT-018159)",
        add = "reg.line",  # Add regressin line
        add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
        conf.int = TRUE, # Add confidence interval
        ggtheme = theme_bw()) +
    annotate("text", x = 2.5,y = 6.1,label = paste0("R = ", signif(res$estimate, 2), "\n", "FDR = ", FDR))

pdf("FigureS5.pdf", width = 9, height = 3)
cowplot::plot_grid(p1_cell,p2_cell, p3_cell, nrow = 1, align = 'hv')
dev.off()
```
