# Figure S2

## Figure S2ABC

```R
library(GenomicRanges, BSgenome.Hsapiens.UCSC.hg19)
library(pmsignature)
library(ggplot2)
library(cowplot)
library(parallel)

RData_list = list.files(".", "*.RData$", full.names=T)

estimate_mat <- as.data.frame(t(sapply(RData_list, function(x){
  load(x)
   return(c(Signature_num, likelihood, error, max.corr))
  })))

colnames(estimate_mat) <- c("K", "likelihood", "error", "max.corr")

p1 = ggplot(estimate_mat) + geom_line(aes(x = K, y = likelihood)) +
  geom_point(aes(x = K, y = likelihood), size = 4, shape = 20) +
  xlab("Signature Number") + ylab("Log Likelihood") +
  scale_x_continuous(breaks = 2:13) +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 12, face = "bold"),
    aspect.ratio = 1)

p2 = ggplot(estimate_mat) + geom_line(aes(x = K, y = error)) +
  geom_point(aes(x = K, y = error), size = 4, shape = 20) +
  xlab("Signature Number") + ylab("Bootstrap Error") +
  scale_x_continuous(breaks = 2:13) +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 12, face = "bold"),
    aspect.ratio = 1)

p3 = ggplot(estimate_mat) + geom_line(aes(x = K, y = max.corr)) +
  geom_point(aes(x = K, y = max.corr), size = 4, shape = 20) +
  xlab("Signature Number") + ylab("Maximum correlation coefficient") +
  scale_x_continuous(breaks = 2:13) +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 12, face = "bold"),
    aspect.ratio = 1)

pdf("FigureS2ABC.pdf", width = 10,height = 3)
plot_grid(p1, p2, p3, nrow = 1, align = "hv")
dev.off()
```

## Figure S2D

```R
file_cosmicSig = "sigProfiler_SBS_signatures_2019_05_22_Vignettes.csv"
cosmicSig <- read.table(file_cosmicSig, header = T, sep = ",", stringsAsFactors = F)
rownames(cosmicSig) = paste0(substr(cosmicSig$SubType,1,1),"[", cosmicSig$Type, "]", substr(cosmicSig$SubType,3,3))
cosmicSig = cosmicSig[, grep("\\d",colnames(cosmicSig))]
colnames(cosmicSig) = gsub("^[A-Za-z]+","", colnames(cosmicSig))

sig_RData="signature_7.RData"

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

plotHeatmap = function(sig_mat, cosmicSig, cosmicSig_list = NULL, plot=T, fontsize = 10, CS = 0.75, sortSig = NULL, plot_margin = c(0,0,0,0)){
  Signature_num = ncol(sig_mat)
  rownames(sig_mat)<- sub("([ACGTN])([ACGTN]) ([ATCG]).([ATCG])", "\\3[\\1>\\2]\\4", rownames(sig_mat))
  sig_mat <-cbind(sig_mat, cosmicSig[match(rownames(sig_mat), rownames(cosmicSig)),])
  corr.sig <- sapply(colnames(sig_mat)[1:Signature_num], function(x){
    return(sapply(colnames(cosmicSig), function(y){
      return(cosine(sig_mat[, x], sig_mat[, y]))
    }))
  })
  if(!plot){
    return(corr.sig)
  }
  m <- melt(corr.sig, varnames = c("A","B"))
  m$lab = signif(m$value, 2)
  m$lab[which(m$lab < CS)] <- ""
  if(!is.null(sortSig)){
    m$B <- factor(m$B, levels = sortSig)
  }
  p_heat <- ggplot(m, aes(x = B, y = A, fill = value)) + geom_tile(color="white", size=0.1) +
    geom_text(aes(label=lab), size = fontsize * 4/14, color = "white") +
    scale_fill_gradient2(low = "white", high = "grey10") +
    labs(y = "COSMIC SBS Signatures", x = "") +
    scale_y_discrete(position = "right", breaks = levels(factor(m$A)),
                     labels = levels(factor(m$A))) +
    coord_flip() +
    theme(axis.text.x = element_text(hjust = 0, size = fontsize, angle = 0),
          axis.text.y = element_text(hjust = 0.5, size = fontsize),
          axis.title = element_text(size = fontsize * 1.2),
          axis.line = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          legend.title=element_blank(),
          legend.position="top",
          legend.justification = "center",
          legend.text=element_text(size = fontsize),
          plot.margin=margin(c(0,0,0,0), unit = "lines"))
  return(p_heat)
}

sig_mat <- getMatSignatures(Param = res_pmSig$Param, BG = res_pmSig$BG_prob, lenMotif = 3,
                            Signame = c("EPP1","Tobacco", "APOBEC", "NpCpG", "dMMR", "EPP2", "Background"))
sig_mat = sig_mat[,c("NpCpG", "APOBEC","dMMR","Tobacco","EPP1","EPP2", "Background")]
colnames(sig_mat) = paste0("MS", 1:7)

p_heat = plotHeatmap(sig_mat, cosmicSig, fontsize = 10, plot_margin = c(0,0,0,0),
            CS = 0.75, sortSig = rev(paste0("MS", 1:7)))

pdf("FigureS2D.pdf", width = 13.5, height = 3.2)
p_heat
dev.off()
```
