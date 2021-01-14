Args <- commandArgs(T)

input_dir = Args[1]        # dir

library(GenomicRanges, BSgenome.Hsapiens.UCSC.hg19)
library(pmsignature)
library(ggplot2)
library(cowplot)
library(parallel)

RData_list = list.files(input_dir, "*.\\d.RData$", full.names = T)

estimate_mat <- as.data.frame(t(sapply(RData_list, function(x){
  load(x)
   return(c(Signature_num, likelihood, error, max.corr))
  })))

colnames(estimate_mat) <- c("K", "likelihood", "error", "max.corr")

p1 = ggplot(estimate_mat) + geom_line(aes(x = K, y = likelihood)) + 
  geom_point(aes(x = K, y = likelihood), size = 4, shape = 20) +
  # labs(title="The Estimated Log Likelihood") + 
  xlab("Signature Number") + ylab("Log Likelihood") +
  scale_x_continuous(breaks = 2:13) +
  theme(axis.text = element_text(size = 15),
    axis.title = element_text(size = 20, face = "bold"),
    plot.title = element_text(size = 20, face = "bold"))

p2 = ggplot(estimate_mat) + geom_line(aes(x = K, y = error)) + 
  geom_point(aes(x = K, y = error), size = 4, shape = 20) +
  # labs(title = "The Estimated Bootstrap Error") +
  xlab("Signature Number") + ylab("Bootstrap Error") +
  scale_x_continuous(breaks = 2:13) +
  theme(axis.text = element_text(size = 15),
    axis.title = element_text(size = 20, face = "bold"),
    plot.title = element_text(size = 20, face = "bold"))

p3 = ggplot(estimate_mat) + geom_line(aes(x = K, y = max.corr)) + 
  geom_point(aes(x = K, y = max.corr), size = 4, shape = 20) +
  # labs(title="The Estimated Maximum correlation coefficient") +
  xlab("Signature Number") + ylab("Maximum correlation coefficient") +
  scale_x_continuous(breaks = 2:13) +
  theme(axis.text = element_text(size = 15),
    axis.title = element_text(size = 20, face = "bold"),
    plot.title = element_text(size = 20, face = "bold"))


pdf(paste0(input_dir, "/est.pdf"), width = 20,height = 6)
plot_grid(p1, p2, p3, nrow = 1, align = "hv") 
dev.off()
