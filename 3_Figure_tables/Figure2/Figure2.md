# Figure2

## Figure 2A

```R
sig_RData="../FigureS2/signature_7.RData"

loadPmSig <- function(pmSigRData){
    load(pmSigRData)
    return(list(G = G, BG_prob = BG_prob, Param = Param))
}

res_pmSig = loadPmSig(sig_RData)

## get signatures matrix
getMutSignatures <- function(Param, sigIndex, lenMotif = 3) {
  numBases = Param@flankingBasesNum
  bases = c("A","C","G","T")
  subs = c("CA", "CG", "CT", "TA", "TC", "TG")
  vF <- Param@signatureFeatureDistribution[sigIndex, , ]
  sub6 = vF[1, ]
  names(sub6) = subs
  context = vF[2:(numBases), 1:4]
  colnames(context) = bases
  
  motif_16 = context[3,] %*% t(context[2,])
  rownames(motif_16) = bases
  motif_16 = melt(motif_16, varnames = c("A","B"))
  rownames(motif_16) = paste0(motif_16$B,".", motif_16$A)
  motif_16[,1:2] = NULL
  
  motif_96 = melt(data.matrix(motif_16) %*% t(sub6), varnames = c("A","B"))
  rownames(motif_96) = paste(motif_96$B, motif_96$A)
  motif_96[,1:2] = NULL
  colnames(motif_96) <- paste0("S", sigIndex)
  
  if(lenMotif == 3){
    return(as.matrix(motif_96))
  }else{
    motif_384 = data.matrix(motif_96) %*% context[1,]
    colnames(motif_384) = bases
    motif_384 = melt(motif_384, varnames = c("A","B"))
    motif_1536 = data.matrix(motif_384[,3]) %*% context[4, ]
    colnames(motif_1536) = bases
    rownames(motif_1536) = paste(motif_384$A, motif_384$B)
    motif_1536 = melt(motif_1536, varnames = c("A","B"))
    motif_1536 = cbind(motif_1536, str_split(motif_1536$A, pattern = " ", simplify = T))
    rownames(motif_1536) = paste(motif_1536$`1`, paste0(motif_1536$`3`, motif_1536$`2`, motif_1536$B))
    motif_1536[,c(1:2,4:6)] = NULL
    colnames(motif_1536) <- paste0("S", sigIndex)
    return(as.matrix(motif_1536))
  }
}

getMatSignatures <- function(Param, BG = NULL, reorderSig = NULL, Signame = NULL, lenMotif = 3){
  signatureNum <- Param@signatureNum
  subs = c("CA", "CG", "CT", "TA", "TC", "TG")
  bases = c("A", "C", "G", "T")
  strands = c("+", "-")
  if (Param@isBackGround){
    BG_signatureNum = signatureNum - 1
  }else{
    BG_signatureNum = signatureNum
  }
  
  if (is.null(reorderSig)) {
    reorderSig <- 1:BG_signatureNum
  }
  sigOrder <- reorderSig
  
  sig_mm <- getMutSignatures(Param, sigOrder[1],lenMotif)
  if(BG_signatureNum > 1){
    for(i in 2:BG_signatureNum){
      sig_mm <- cbind(sig_mm, getMutSignatures(Param, sigOrder[i], lenMotif))
    }
  }
  if(!is.null(BG)){
    bg_name = unique(substr(names(BG), 1, 9))
    bg_name_2 = sort(c(paste(bg_name, "1", sep=","), paste(bg_name, "2", sep=",")))
    BG_new = BG[match(bg_name_2, names(BG))]
    names(BG_new) = bg_name_2
    BG <- matrix(BG_new, ncol=2, byrow=T)
    rownames(BG) <- bg_name
    BG[is.na(BG)] <- 0
    BG <- BG/sum(BG) # normalize to 1
    BG <- as.data.frame(rowSums(BG))
    colnames(BG) = "BG"
    if(lenMotif == 3){
      BG$name <- sapply(rownames(BG),function(x){
        name = as.numeric(unlist(strsplit(x, ",")))
        return(paste(subs[name[1]], paste0(bases[name[3]], ".", bases[name[4]])))
      })
      BG = BG %>% group_by(name) %>% summarise(BG=sum(BG)) %>% as.data.frame()
    }else{
      BG$name <- sapply(rownames(BG),function(x){
        name = as.numeric(unlist(strsplit(x, ",")))
        return(paste(subs[name[1]], paste0(paste(bases[name[2:3]], collapse = ""), ".",
                                        paste(bases[name[4:5]], collapse = ""))))
      })
    }
    sig_mm <- cbind(sig_mm, Sbg=BG$BG[match(rownames(sig_mm), BG$name)])
  }
  if (!is.null(Signame)){
    colnames(sig_mm) <- Signame
  }
  return(sig_mm)
}

sig_mat5 <- getMatSignatures(res_pmSig$Param, BG = res_pmSig$BG_prob, lenMotif = 5,
                             Signame = c("EPP1","Tobacco", "APOBEC", "NpCpG", "dMMR", "EPP2", "BG"))
sig_mat5 = sig_mat5[,c("NpCpG", "APOBEC","dMMR","Tobacco","EPP1","EPP2", "BG")]
colnames(sig_mat5) = paste0("MS", 1:7)

## plot
plotSignatures <- function(sig_mat, lenMotif = 3, charSize = 10, scales = "fixed",
                           scale_y_sqrt = FALSE, showMotif = T, SNVs = T, plot_margin = c(0,0,0,0)){
    if(SNVs){
        level = c("SNVs", colnames(sig_mat))
        sig_mat <- cbind(sig_mat, SNVs = rowSums(sig_mat))
    }else{
        level = colnames(sig_mat)
    }
    w_df = melt(sig_mat, varnames = c("motif", "signature"))
    w_df$alteration = sub("([ACGTN])([ACGTN]) .+", "\\1>\\2", w_df$motif)
    # w_df$context = sub("[ACGTN][ACGTN] (.+)", "\\1", w_df$motif)
    w_df$context = sub("([ACGTN])([ACGTN]) ([ACGTN]*).([ACGTN]*)", "\\3\\1\\4", w_df$motif)
    w_df$signature = factor(w_df$signature, levels = level)
    if(lenMotif == 3){
        xlab = "Tri-nucleotide Sequence Motifs"
    }else{
        xlab = "Five nucleotide Sequence Motifs"
    }
    p = ggplot(w_df)
    p = p + geom_bar(aes_string(x = "context", y = "value", fill = "alteration"),
                     stat = "identity", position = "identity")
    # stat = "identity", position = "identity",  colour="black",size=0.1)
    p = p + facet_grid(signature ~ alteration, scales = scales)
    p = p + theme_bw()
    p = p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = charSize*0.6),
                  axis.text.y = element_text(hjust = 0.5, size = charSize),
                  axis.title.y = element_text(size = charSize*1.2),
                  panel.background = element_blank(),
                  panel.border = element_rect(color = "gray50"),
                  panel.spacing = unit(0.1,"lines"),
                  strip.text.x = element_text(size = charSize),
                  strip.text.y = element_text(size = charSize, angle = 0),
                  strip.placement = "outside",
                  strip.background.y = element_blank(),
                  axis.ticks.x = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  legend.position = "none",
                  plot.margin=margin(plot_margin, unit = "lines"))
    p = p + scale_fill_d3()
    # scale_fill_manual(values = alpha(c("Sky Blue","black", "firebrick3",
    #                                            "gray","lightgreen", "lightpink")))
    p = p + labs(x = xlab, y = "Activities")
    p = p + scale_y_continuous(n.breaks=3)
    if(scale_y_sqrt == TRUE){
        p = p + scale_y_sqrt()
    }
    if(!showMotif){
        p = p + theme(axis.text.x = element_blank())
    }
    return(p)
}


p_sigs5 = plotSignatures(sig_mat5, lenMotif = 5, charSize = 10, scales = "free",
                         scale_y_sqrt = FALSE, showMotif = F, SNVs = F, plot_margin = c(0,1,1,0))

pdf("Figure2A.pdf", width = 6, height = 3)
p_sigs5
dev.off()
```

## Figure 2B and 2C

```R
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

plot_qqnorm = function(pmSigMat, pmLogSigMat){
  p_theme = theme_bw() + theme(plot.title = element_text(hjust = 0.5, size = 14),
      axis.text = element_text(size = 12),
      legend.position = c(0.2, 0.8),
      legend.background = element_blank(),
      strip.placement = "outside",
      strip.text = element_text(size = 14),
      strip.background = element_blank(),
      axis.title = element_text(size = 14))
  p_raw = ggplot(melt(pmSigMat), aes(sample= value, color = variable)) + stat_qq() + stat_qq_line() +
    scale_color_d3() + labs(x = "Theoretical Quantiles", y = "MS activities", color = "Mutational signatures") + p_theme
  p_relative = ggplot(melt(pmLogSigMat), aes(sample= value, color = variable)) + stat_qq() + stat_qq_line() +
    scale_color_d3() + labs(x = "Theoretical Quantiles", y = "MP activities", color = "Mutational propensities") + p_theme
  return(list(p_raw = p_raw, p_relative = p_relative))
}

p_qqnorm = plot_qqnorm(pmSigMat, pmLogSigMat)
```

## Figure 2D

```R
mem = read.table("TCGA_5102_TMB.txt", header=T,sep="\t",stringsAsFactors = F)
rownames(mem) = mem$sample

plot_sample_signature = function(mem, Mat, pmSigMat, plot=T, plot_margin = c(0,0,0,0), ylab = "Activities"){
    mem$BG = pmSigMat[mem$sample, "MS7"]
    mat <- mem[, c("sample", "cancer", "TMB", "BG")] %>% group_by(cancer) %>%
      mutate(TMB_median = median(TMB)) %>% arrange(TMB_median, desc(BG)) %>% ungroup %>%
      mutate(cancer = factor(cancer, levels = unique(cancer)),
            sample = factor(sample, levels = sample))
    mat$TMB_group = "middle"
    mat$TMB_group[mat$TMB>20] = "high"
    mat$TMB_group[mat$TMB<5] = "low"
    mat$TMB_group[mat$TMB<1] = "ultralow"
    mat$TMB_group = factor(mat$TMB_group, levels = c("ultralow","low","middle","high"))
    mat <- cbind(mat, Mat[as.character(mat$sample), ])
    if(!plot){
      return(mat)
    }
    mat = mat[, -c(grep("TMB",colnames(mat)), grep("BG", colnames(mat)))]
    m = na.omit(melt(mat))
    m$source = "TCGA"
    p <- ggplot(m) + geom_col(aes(x = sample, y = value, fill = variable), width = 1) +
      labs(x = "", y = ylab) +
      scale_fill_d3() +
      scale_y_continuous(expand = c(0,0)) +
      facet_grid(source~cancer, scales = "free_x") +
      guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
      theme_bw() +
      theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none",
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0.1, "lines"),
        plot.background = element_blank(),
        plot.margin = margin(plot_margin, unit = "lines"))
    return(p)
}

p_sample_MPs = plot_sample_signature(mem, pmLogSigMat, pmSigMat, ylab = "Relative activities")

mat_BG_TMB = plot_sample_signature(mem, pmLogSigMat, pmSigMat, ylab = "Relative activities", plot = F)

p_TMB = ggplot(mat_BG_TMB, aes(sample,y=1, fill=TMB_group)) + geom_tile() +
      scale_fill_grey(start = 0.8, end=0.2) + labs(y = "TMB", fill = "TMB levels") +
      facet_grid(~cancer, scales = "free_x") +
      theme(axis.text = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(angle=0,size=12,vjust=0.5),
          axis.line = element_blank(),
          panel.background = element_blank(),
          plot.background = element_blank(),
          legend.background = element_blank(),
          panel.spacing = unit(0.1, "lines"),
          axis.ticks = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          legend.title=element_blank(),
          legend.position="top",
          legend.direction = "horizontal",
          legend.justification = "center",
          legend.text=element_text(size = 10),
          plot.margin=margin(c(0,0,0,0), unit = "lines"))


cell_pmLogSigMat = read.table("CCLE_pmLogSigMat.txt", header = T, sep = "\t", stringsAsFactors = F)

CCLE_TCGA = rbind(cbind(pmLogSigMat, cancer = mem[rownames(pmLogSigMat), "cancer"], type = "TCGA"),
                cbind(cell_pmLogSigMat[,c(1:7)], type = "CCLE"))
CCLE_TCGA = melt(CCLE_TCGA, id.vars = c("cancer", "type"))

get_CCLE_TCGA_mat = function(CCLE_TCGA, mem){
    CCLE_TCGA_meanSig = CCLE_TCGA %>% group_by(type, cancer, variable) %>% mutate(num = sum(!is.na(value))) %>%
        filter(num > 2) %>% summarise(value=mean(value, na.rm = T)) %>%
        acast(cancer + variable ~ type) %>% na.omit() %>% as.data.frame()
    CCLE_TCGA_meanSig = cbind(CCLE_TCGA_meanSig,
        as.data.frame(str_split(rownames(CCLE_TCGA_meanSig), pattern = "_", simplify = T), stringsAsFactors = F))
    colnames(CCLE_TCGA_meanSig)[3:4] = c("cancer","sig")
    res_cosine = PharmacoGx::cosinePerm(CCLE_TCGA_meanSig$TCGA, CCLE_TCGA_meanSig$CCLE,include.perm=T, nperm = 1e6)
    res_cosine_cancers = sapply(unique(CCLE_TCGA_meanSig$cancer), function(x){
      res = cosine(CCLE_TCGA_meanSig$TCGA[CCLE_TCGA_meanSig$cancer==x],
                CCLE_TCGA_meanSig$CCLE[CCLE_TCGA_meanSig$cancer==x])
      return(res)
      })
    CCLE_TCGA$label = paste0(CCLE_TCGA$cancer,"\n", signif(res_cosine_cancers[CCLE_TCGA$cancer], digits=3))
    CCLE_TCGA$label = gsub("NA","",CCLE_TCGA$label)

    mat <- mem[, c("sample", "cancer", "TMB")] %>% group_by(cancer) %>%
      mutate(TMB_median = median(TMB)) %>% arrange(TMB_median, TMB) %>% ungroup %>%
      mutate(cancer = factor(cancer, levels = unique(cancer)),
          sample = factor(sample, levels = sample))
    CCLE_TCGA$cancer = factor(CCLE_TCGA$cancer, levels = levels(mat$cancer))
    CCLE_TCGA$label = factor(CCLE_TCGA$label, levels = unique(CCLE_TCGA$label[order(CCLE_TCGA$cancer)]))
    CCLE_TCGA_meanSig$cancer = factor(CCLE_TCGA_meanSig$cancer, levels = levels(mat$cancer))
    return(list(CCLE_TCGA = CCLE_TCGA, CCLE_TCGA_meanSig = CCLE_TCGA_meanSig, res_cosine = res_cosine))
}

CCLE_TCGA_mat = get_CCLE_TCGA_mat(CCLE_TCGA, mem)

p_TCGA_boxplot_MP = ggplot(CCLE_TCGA_mat$CCLE_TCGA %>% filter(type=="TCGA"), aes(x = variable, y = value, fill = variable, color = variable)) +
    geom_violin(outlier.size = 0.5) +
    facet_grid(type ~ cancer, scales = "free_y") +
    guides(fill = guide_legend(nrow = 1)) +
    ylab("Relative activities") + scale_fill_d3() +
    theme_light() + scale_color_d3() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size = 12),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14),
          legend.position = 'bottom',
          legend.direction = 'horizontal',
          legend.title = element_blank(),
          legend.key.size=unit(1,'cm'),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.spacing = unit(0.1, 'lines'),
          plot.background = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank())

pdf("Figure2D.pdf", width = 10, height = 6)
plot_grid(p_TMB, p_sample_MPs, p_TCGA_boxplot_MP,
        ncol = 1,rel_heights = c(0.5,1.5,1.8), align = "v", axis = 'rl')
dev.off()
```

## Figure 2E

```R
mat_BG_TMB_raw = plot_sample_signature(mem, pmSigMat, pmSigMat, plot = F)
colnames(mat_BG_TMB_raw)[grep("MS", colnames(mat_BG_TMB_raw))] = gsub("MS", "MP",colnames(mat_BG_TMB_raw)[grep("MS", colnames(mat_BG_TMB_raw))])
mat_BG_TMB_raw$MP7 = NULL

mat_merge = rbind(cbind(mat_BG_TMB_raw, type = "raw"),
      cbind(mat_BG_TMB, type = "log"))

sig_TMB_cor = lapply(paste0("MP",1:6), function(s){
  b = mat_merge[,c("cancer","TMB", s, "type")]
  colnames(b)[3] = "sig"
  res = lapply(unique(b$cancer), function(x){
    b1 = b[which(b$cancer == x & b$type == "raw"), ]
    res1 = cor.test(b1$sig, b1$TMB)
    b2 = b[which(b$cancer == x & b$type == "log"), ]
    res2 = cor.test(b2$sig, b2$TMB)
    res = rbind(data.frame(sig = s, cancer = x, est_raw = res1$estimate, p_raw = res1$p.value,
      est_log = res2$estimate, p_log = res2$p.value))
    return(res)
  })
  res = do.call(rbind, res)
  return(res)
})
sig_TMB_cor = do.call(rbind, sig_TMB_cor)

library(ggpubr)
sig_TMB_cor_res = cor.test(sig_TMB_cor$est_raw, sig_TMB_cor$est_log, method = "spearman")
p_sig_TMB_cor_cor = ggplot(sig_TMB_cor, aes(est_raw, est_log)) + geom_point(aes(color = cancer), size=3) +
  geom_smooth(method = "lm",linetype = "dashed") +
  labs(x="MS effect sizes on TMB", y = "MP effect sizes on TMB",
      color = "", shape="") +
  # title = "Signature relative activities\nin TCGA and CCLE"
  theme_bw() + scale_color_simpsons() +
  guides(color = guide_legend(ncol = 3)) +
  theme(legend.position = "right",
    legend.text = element_text(size=12),
    legend.box.margin = margin(0,0,0,0,"lines"),
    axis.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(size = 14),
    aspect.ratio = 1) +
  annotate("text", x = -Inf, y = Inf, size = 5, hjust = -0.2, vjust = 2,
          label= paste("rho =", signif(sig_TMB_cor_res$estimate, digits = 3), "\n p =",
                        signif(sig_TMB_cor_res$p.value, digits = 3)))

pdf("Figure2E.pdf", width = 8, height = 4)
p_sig_TMB_cor_cor
dev.off()
```

## Figure 2F

```R
immu = read.table("TCGA_ITH.txt", header=T,sep="\t",stringsAsFactors = F)

mat_merge$ITH = immu$Intratumor.Heterogeneity[match(mat_merge$sample, immu$TCGA.Participant.Barcode)]
mat_merge = mat_merge %>% group_by(cancer, type) %>% mutate(ITH_group = ifelse(ITH < median(ITH[ITH > 0], na.rm = T), "low", "high"))

sig_ITH_cor = lapply(paste0("MP",1:6), function(s){
  b = mat_merge[,c("cancer","ITH", s, "type")]
  colnames(b)[3] = "sig"
  res = lapply(unique(b$cancer), function(x){
    b1 = b[which(b$cancer == x & b$type == "raw"), ]
    res1 = cor.test(b1$sig, b1$ITH)
    b2 = b[which(b$cancer == x & b$type == "log"), ]
    res2 = cor.test(b2$sig, b2$ITH)
    res = rbind(data.frame(sig = s, cancer = x, est_raw = res1$estimate, p_raw = res1$p.value,
      est_log = res2$estimate, p_log = res2$p.value))
    return(res)
  })
  res = do.call(rbind, res)
  return(res)
})
sig_ITH_cor = do.call(rbind, sig_ITH_cor)

sig_ITH_cor_res = cor.test(sig_ITH_cor$est_raw, sig_ITH_cor$est_log, method = "spearman")
p_sig_ITH_cor_cor = ggplot(sig_ITH_cor, aes(est_raw, est_log)) + geom_point(aes(color = cancer), size=3) +
    geom_smooth(method = "lm",linetype = "dashed") +
    labs(x="MS effect sizes on ITH", y = "MP effect sizes on ITH",
         color = "", shape="") +
    # title = "Signature relative activities\nin TCGA and CCLE"
    theme_bw() + scale_color_simpsons() +
    guides(color = guide_legend(ncol = 3)) +
    theme(legend.position = "right",
        legend.text = element_text(size=12),
        legend.box.margin = margin(0,0,0,0,"lines"),
        axis.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 14),
        aspect.ratio = 1) +
    annotate("text", x = -Inf, y = Inf, size = 5, hjust = -0.2, vjust = 2,
             label= paste("rho =", signif(sig_ITH_cor_res$estimate, digits = 3), "\n p =",
                          signif(sig_ITH_cor_res$p.value, digits = 3)))

pdf("Figure2F.pdf", width = 8, height = 4)
p_sig_ITH_cor_cor
dev.off()
```

## Figure 2G

```R
p_CCLE_TCGA_meanSig = ggplot(CCLE_TCGA_mat$CCLE_TCGA_meanSig, aes(CCLE, TCGA)) + geom_point(aes(color = cancer), size=3) +
    geom_smooth(method = "lm",linetype = "dashed") +
    labs(x="Mean MP activities in CCLE", y = "Mean MP activities in TCGA",
         color = "Cancer types", shape="") +
    theme_bw() + scale_color_simpsons() +
    guides(color = guide_legend(ncol = 3)) +
    theme(legend.position = "right",
        legend.text = element_text(size=12),
        legend.box.margin = margin(0,0,0,0,"lines"),
        axis.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 14),
        aspect.ratio = 1) +
    annotate("text", x=-3, y=0, size = 5,
             label= paste("r =", signif(CCLE_TCGA_mat$res_cosine$estimate, digits = 3), "\n p =",
                          signif(CCLE_TCGA_mat$res_cosine$p.value, digits = 3)))

pdf("Figure2G.pdf", width = 8, height = 4)
p_CCLE_TCGA_meanSig
dev.off()
```
