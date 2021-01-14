Args <- commandArgs(T)

input_pm = Args[1]                    # {Cancer_Type}/{Cancer_Type}.pm
Signature_num = as.numeric(Args[2])   # the number of signatures
mc.cores = as.numeric(Args[3])        # the number of cores
outdir = Args[4]                      # the output dir

if(length(Args)==0){
  print("Rscript getSig.all.R input_pm Signature_num process_num output_dir")
  q()
}

out_RData = gsub("pm", paste0(paste0("est.sig.",Signature_num), ".RData"), basename(input_pm))
out_RData = paste0(outdir, "/", out_RData)

library(GenomicRanges, BSgenome.Hsapiens.UCSC.hg19)
library(pmsignature)
library(ggplot2)
library(cowplot)
library(reshape)
library(lsa)
library(parallel)

G <- readMPFile(input_pm, numBases = 5, trDir = TRUE)
BG_prob <- readBGFile(G)

# Deciphering Signatures
for(i in list.files("~/R/packages/pmsignature-devel/R", "*.R", full.names=T)){source(i)}

getPMSignature_parallel = function (mutationFeatureData, K, BG = NULL, numInit = 10, tol = 1e-04, 
    maxIter = 10000) 
{
    if (!is.null(BG)) {
        isBG <- TRUE
        varK <- K - 1
    }else {
        isBG <- FALSE
        varK <- K
        BG <- 0
    }
    sampleNum <- length(slot(mutationFeatureData, "sampleList"))
    fdim <- slot(mutationFeatureData, "possibleFeatures")
    tempL <- -Inf
    tempPar <- c()
    ress = mclapply(1:numInit, function(kkk){
        F <- array(0, c(varK, length(fdim), max(fdim)))
        for (k in 1:varK) {
            for (kk in 1:length(fdim)) {
                F[k, kk, 1:fdim[kk]] <- rgamma(fdim[kk], rep(1, 
                  fdim[kk]))
                F[k, kk, 1:fdim[kk]] <- F[k, kk, 1:fdim[kk]]/sum(F[k, 
                  kk, 1:fdim[kk]])
            }
        }
        Q <- matrix(rgamma(sampleNum * K, 1, 1), K, sampleNum)
        Q <- sweep(Q, 2, apply(Q, 2, sum), `/`)
        p0 <- c(convertToTurbo_F(as.vector(F), fdim, K, isBG), 
            convertToTurbo_Q(as.vector(t(Q)), K, sampleNum))
        Y <- list(list(sampleNum, fdim, slot(mutationFeatureData, 
            "featureVectorList"), slot(mutationFeatureData, "countData")), 
            K, isBG, BG)
        res1 <- mySquareEM(p0, Y, tol = tol, maxIter = maxIter)
        cat(paste("#trial: ", sprintf("%2d", kkk), "; #iteration: ", 
            sprintf("%4d", as.integer(res1$itr)), "; time(s): ", 
            sprintf("%4.2f", res1$elapsed.time), "; convergence: ", 
            res1$convergence, "; loglikelihood: ", sprintf("%.4f", 
                res1$value.objfn), "\n", sep = ""))
        return(c(res1$value.objfn, res1$par))
    }, mc.cores = mc.cores)

    tempL = max(sapply(ress, function(x){x[1]}))
    tempPar = unlist(ress[which.max(sapply(ress, function(x){x[1]}))])[-1]
    
    lenF <- varK * (sum(fdim) - length(fdim))
    lenQ <- sampleNum * (K - 1)
    F <- convertFromTurbo_F(tempPar[1:lenF], fdim, K, isBG)
    Q <- convertFromTurbo_Q(tempPar[(lenF + 1):(lenF + lenQ)], 
        K, sampleNum)
    dim(F) <- c(varK, length(fdim), max(fdim))
    dim(Q) <- c(sampleNum, K)
    return(new(Class = "EstimatedParameters", type = slot(mutationFeatureData, 
        "type"), flankingBasesNum = slot(mutationFeatureData, 
        "flankingBasesNum"), transcriptionDirection = slot(mutationFeatureData, 
        "transcriptionDirection"), possibleFeatures = slot(mutationFeatureData, 
        "possibleFeatures"), sampleList = slot(mutationFeatureData, 
        "sampleList"), signatureNum = as.integer(K), isBackGround = isBG, 
        signatureFeatureDistribution = F, sampleSignatureDistribution = Q, 
        loglikelihood = tempL))
}

bootPMSignature_parallel <- function(mutationFeatureData, Param0, bootNum = 10, BG = NULL, 
    tol = 1e-2, maxIter = 10000) {
  K <- slot(Param0, "signatureNum")
  isBG <- slot(Param0, "isBackGround")

  if (isBG == TRUE) {
    if (!is.null(BG)) {
      varK <- K - 1
    } else {
      stop(paste("The input parameter is estimated using a background signature.\n",
                 "Please specify the same background signature."))
    }
  } else {
    if (!is.null(BG)) {      
      warning(paste("The input parameter is estimated without using a background signature.\n",
                    "Specified background signature is ignored."))
    }
    varK <- K
    BG <- 0
  }

  sampleNum <- length(slot(mutationFeatureData, "sampleList"))
  fdim <- slot(mutationFeatureData, "possibleFeatures")
  countData_org <- slot(mutationFeatureData, "countData")
  bootData <- countData_org

  F0 <- slot(Param0, "signatureFeatureDistribution")
  Q0 <- slot(Param0, "sampleSignatureDistribution")

  tempL <- -Inf
  tempPar <- c()

  sqF <- array(0, c(bootNum, varK, length(fdim), max(fdim)))
  sqQ <- array(0, c(bootNum, nrow(Q0), ncol(Q0)))

  sq_list = mclapply(1:bootNum, function(bbb){
      ##########
      # This part is under construction!!!!
      # bootData violates the validity rules of the mutation feature class... I don't like this..
      tempG <- table(sample(1:length(countData_org[3,]), sum(countData_org[3,]), replace=TRUE, prob= countData_org[3,] / sum(countData_org[3,]) ))
      bootData[3, ] <- 0
      bootData[3, as.integer(names(tempG))] <- tempG
      ##########
      
      p0 <- c(convertToTurbo_F(as.vector(F0), fdim, K, isBG), convertToTurbo_Q(as.vector(t(Q0)), K, sampleNum))
      Y <- list(list(sampleNum, fdim, slot(mutationFeatureData, "featureVectorList"), bootData), K, isBG, BG)
      
      res1 <- mySquareEM(p0, Y, tol = tol, maxIter = maxIter)
      cat(paste("#trial: ", sprintf("%2d", bbb), 
                "; #iteration: ", sprintf("%4d", as.integer(res1$itr)), 
                "; time(s): ", sprintf("%4.2f", res1$elapsed.time), 
                "; convergence: ", res1$convergence,
                "; loglikelihood: ", sprintf("%.4f", res1$value.objfn), "\n", sep=""
      ))
      
      tempPar <- res1$par
      lenF <- varK * (sum(fdim) - length(fdim))
      lenQ <- sampleNum * (K - 1)
      F <- convertFromTurbo_F(res1$par[1:lenF], fdim, K, isBG)
      Q <- convertFromTurbo_Q(res1$par[(lenF + 1):(lenF + lenQ)], K, sampleNum)
      dim(F) <- c(varK, length(fdim), max(fdim))
      dim(Q) <- c(sampleNum, K)
      
      sqF_l = lapply(1:varK, function(k){(F[k,,] - F0[k,,])^2})
      sqQ_l = lapply(1:sampleNum, function(n){(Q[n,] - Q0[n,])^2})
      return(list(bbb, sqF_l, sqQ_l))
      }, mc.cores = mc.cores)
  
  for(i in 1:bootNum){
    a = sq_list[[i]]
    for(j in 1:length(a[[2]])){
      sqF[a[[1]], j, , ] <- a[[2]][[j]]
    }
    sqQ[a[[1]],,] <- matrix(unlist(a[[3]]), nrow=length(a[[3]]), byrow=T)
  }
  return(list(sqF, sqQ))
}

Param <- getPMSignature_parallel(G, K = Signature_num, numInit = 10, BG = BG_prob)
error = mean(bootPMSignature_parallel(G, Param, bootNum = 10, BG = BG_prob)[[1]])
likelihood = Param@loglikelihood
max.corr = max(cor(Param@sampleSignatureDistribution) - 2 * diag(rep(1, Signature_num)))

save.image(out_RData)
