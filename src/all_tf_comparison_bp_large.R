library(multiseq)

nbase = 2^14

normalize = function(x) x/sum(x)

setwd("~/projects/multiseq_analysis")

source("src/rocplot.R")

file.names=list.files("data/", pattern = "*.txt.gz")
chip.file.names=file.names[seq(5,310,5)]

for(j in 1:length(chip.file.names)){
  factor.name.short = strsplit(chip.file.names[j],"_")[[1]][1]
  factor.name = paste(strsplit(chip.file.names[j],"_")[[1]][1], strsplit(chip.file.names[j],"_")[[1]][2],sep="_")
  dnase.file.name = paste(factor.name, "dnase_rowsums_large.txt", sep = "_")
  dnase.file.path = file.path("data",dnase.file.name)
  if(!file.exists(dnase.file.path)) next
  dcut=as.matrix(read.table(dnase.file.path))
  cinfo=read.table(file.path("data",chip.file.names[j]),header=TRUE)
  mapsites = read.table(file.path("data", paste(factor.name, "mappability_flags_large.txt.gz", sep = "_")))
  dcut = dcut[mapsites[, 1] == 1 & cinfo[, 5] >= 10, ]
  cinfo = cinfo[mapsites[, 1] == 1 & cinfo[, 5] >= 10, ]
  ccount=cinfo[,6]
  peak.info = read.table(file.path("data","macs_peaks", paste(factor.name.short, "q1pc", "peaks.bed", sep = "_")))
  peak.ind = 0
  for(i in 1:dim(cinfo)[1]){
    peak.ind[i] = sum((as.character(peak.info[, 1]) == as.character(cinfo[i, 1])) & (peak.info[, 2] < cinfo[i, 2]) & (peak.info[, 3] > cinfo[i, 3])) > 0
  }
  
  unbound.ind = !(peak.ind) & cinfo[, 6] < median(cinfo[, 6])
  bound.ind = as.logical(peak.ind)

  na.ind = is.na(rowSums(dcut))
  dcut = dcut[!na.ind, ]
  cinfo = cinfo[!na.ind, ]
  unbound.ind = unbound.ind[!na.ind]
  bound.ind = bound.ind[!na.ind]

  print("Read data")
  
  setwd("results/roc/large")

  
  pdf(paste("roc", factor.name, "binary_large.pdf", sep = "_"))
  
  roc.res.dcut = rocplot(dcut[, 1][unbound.ind], dcut[, 1][bound.ind])
  plot(roc.res.dcut$fpr,roc.res.dcut$tpr, type = 'l', xlim = c(0, 1), ylim = c(0, 1), xlab = "FPR", ylab = "TPR", main = "#Dnase cuts, bound vs unbound")
  
  roc.res.dcut = rocplot(dcut[, 2][unbound.ind], dcut[, 2][bound.ind])
  lines(roc.res.dcut$fpr,roc.res.dcut$tpr, col = 2, lty = 1)
  
  roc.res.dcut = rocplot(dcut[, 3][unbound.ind], dcut[, 3][bound.ind])
  lines(roc.res.dcut$fpr,roc.res.dcut$tpr, col = 3, lty = 1)
  
  roc.res.dcut = rocplot(dcut[, 4][unbound.ind], dcut[, 4][bound.ind])
  lines(roc.res.dcut$fpr,roc.res.dcut$tpr, col = 4, lty = 1)
  
  roc.res.dcut = rocplot(dcut[, 5][unbound.ind], dcut[, 5][bound.ind])
  lines(roc.res.dcut$fpr,roc.res.dcut$tpr, col = 5, lty = 1)
  
  legend("bottomright", legend = c("2^10-Dcuts", "2^11-Dcuts", "2^12-Dcuts", "2^13-Dcuts", "2^14-Dcuts"), lty = rep(1, 5), col = 1:5)
  dev.off()
  
  
  #   pdf(paste("roc", factor.name, "strand_comp_1024bp.pdf", sep = "_"))
  #   roc.res.ms = rocplot(ccount.post.logmean.list[[4]][unbound.ind], ccount.post.logmean.list[[2]][bound.ind])
  #   plot(roc.res.ms$fpr,roc.res.ms$tpr, type = 'l', xlim = c(0, 1), ylim = c(0, 1), xlab = "FPR", ylab = "TPR", main = "Strand comparison, 1024 bp, multiseq")
  #   
  #   roc.res.ms = rocplot(ccount.post.logmean.for.list[[4]][unbound.ind], ccount.post.logmean.for.list[[2]][bound.ind])
  #   lines(roc.res.ms$fpr,roc.res.ms$tpr, col = 2)
  # 
  #   legend("bottomright", legend = c("Both", "Forward"), lty = rep(1, 2), col = 1:2)
  #   dev.off()
  
  unbound.ind.counts = cinfo[, 6] <= median(cinfo[, 6])
  bound.ind.counts = cinfo[, 6] > quantile(cinfo[, 6], probs = 0.95)  
  
  pdf(paste("roc", factor.name, "binary_counts_large.pdf", sep = "_"))
  
  roc.res.dcut = rocplot(dcut[, 1][unbound.ind.counts], dcut[, 1][bound.ind.counts])
  plot(roc.res.dcut$fpr,roc.res.dcut$tpr, type = 'l', xlim = c(0, 1), ylim = c(0, 1), xlab = "FPR", ylab = "TPR", main = "#Dnase cuts, bound vs unbound using counts only")
  
  roc.res.dcut = rocplot(dcut[, 2][unbound.ind.counts], dcut[, 2][bound.ind.counts])
  lines(roc.res.dcut$fpr,roc.res.dcut$tpr, col = 2, lty = 1)
  
  roc.res.dcut = rocplot(dcut[, 3][unbound.ind.counts], dcut[, 3][bound.ind.counts])
  lines(roc.res.dcut$fpr,roc.res.dcut$tpr, col = 3, lty = 1)
  
  roc.res.dcut = rocplot(dcut[, 4][unbound.ind.counts], dcut[, 4][bound.ind.counts])
  lines(roc.res.dcut$fpr,roc.res.dcut$tpr, col = 4, lty = 1)
  
  roc.res.dcut = rocplot(dcut[, 5][unbound.ind.counts], dcut[, 5][bound.ind.counts])
  lines(roc.res.dcut$fpr,roc.res.dcut$tpr, col = 5, lty = 1)
  
  legend("bottomright", legend = c("2^10-Dcuts", "2^11-Dcuts", "2^12-Dcuts", "2^13-Dcuts", "2^14-Dcuts"), lty = rep(1, 5), col = 1:5)
  dev.off()
  
  

  corr.dcut = c(cor(log(dcut[, 1] + 1), log(cinfo[, 6] + 1)),
                cor(log(dcut[, 2] + 1), log(cinfo[, 6] + 1)),
                cor(log(dcut[, 3] + 1), log(cinfo[, 6] + 1)),
                cor(log(dcut[, 4] + 1), log(cinfo[, 6] + 1)),
                cor(log(dcut[, 5] + 1), log(cinfo[, 6] + 1)))
  
  pdf(paste("corr_pearson", factor.name, "large.pdf", sep = "_"))
  plot(10:14, corr.dcut, type = "b", ylim = c(0, 1), xlab = "window size (log2)", ylab = "correlation", main = "#Dnase cuts, correlation")
  dev.off()
  
  corr.rank.dcut = c(cor(log(dcut[, 1] + 1), log(cinfo[, 6] + 1), method = "spearman"),
                     cor(log(dcut[, 2] + 1), log(cinfo[, 6] + 1), method = "spearman"),
                     cor(log(dcut[, 3] + 1), log(cinfo[, 6] + 1), method = "spearman"),
                     cor(log(dcut[, 4] + 1), log(cinfo[, 6] + 1), method = "spearman"),
                     cor(log(dcut[, 5] + 1), log(cinfo[, 6] + 1), method = "spearman"))
  
  pdf(paste("corr_spearman", factor.name, "large.pdf", sep = "_"))
  plot(10:14, corr.rank.dcut, type = "b", ylim = c(0, 1), xlab = "window size (log2)", ylab = "correlation", main = "#Dnase cuts, rank correlation")
  dev.off()

  
  setwd("../../..")
  
  print(j)
}

