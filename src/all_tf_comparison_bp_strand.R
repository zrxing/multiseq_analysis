library(multiseq)
library(CENTIPEDE)
library(pROC)

nbase = 2^10

normalize = function(x) x/sum(x)

setwd("~/projects/multiseq_analysis")

source("src/rocplot.R")

file.names=list.files("data/", pattern = "*.txt.gz")
chip.file.names=file.names[seq(5,310,5)]
dnase.file.names=file.names[seq(2,307,5)]
map.file.names=file.names[seq(3,308,5)]

for(j in 1:length(chip.file.names)){
  factor.name.short = strsplit(chip.file.names[j],"_")[[1]][1]
  factor.name = paste(strsplit(chip.file.names[j],"_")[[1]][1], strsplit(chip.file.names[j],"_")[[1]][2],sep="_")
  if(as.numeric(file.info(file.path("data",dnase.file.names[j]))[1])>(12*1e+6)) next
  try.test=try(dprof<-read.table(file.path("data",dnase.file.names[j])))
  if(class(try.test)=="try-error") next
  cinfo=read.table(file.path("data",chip.file.names[j]),header=TRUE)
  mapsites = read.table(file.path("data", map.file.names[j]))
  dprof = dprof[mapsites[, 1] == 1 & cinfo[, 5] >= 10, ]
  cinfo = cinfo[mapsites[, 1] == 1 & cinfo[, 5] >= 10, ]
  ccount = cinfo[,6]
  peak.info = read.table(file.path("data","macs_peaks", paste(factor.name.short, "q1pc", "peaks.bed", sep = "_")))
  peak.ind = 0
  for(i in 1:dim(cinfo)[1]){
    peak.ind[i] = sum((as.character(peak.info[, 1]) == as.character(cinfo[i, 1])) & (peak.info[, 2] < cinfo[i, 2]) & (peak.info[, 3] > cinfo[i, 3])) > 0
  }
  
  unbound.ind.all = !(peak.ind) & cinfo[, 6] < median(cinfo[, 6])
  bound.ind.all = as.logical(peak.ind)
 
  print("Read data")

  train.size = round(0.8*dim(dprof)[1])
  set.seed(417)
  train.ind = sample(1:dim(dprof)[1],train.size,replace = FALSE)
  train.ind = (1:dim(dprof)[1])%in%train.ind
  
  dprof.train = dprof[train.ind,]
  dprof.test = dprof[!train.ind,]
  
  ccount.train = ccount[train.ind]
  ccount.test = ccount[!train.ind]
  
  nbase = 2^6
  nbase[2] = 2^7
  nbase[3] = 2^8
  nbase[4] = 2^9
  nbase[5] = 2^10
  
  base.ind.for = list()
  base.ind.rev = list()
  
  base.ind.for[[1]] = 4065:4128 - 3584
  base.ind.rev[[1]] = 12257:12320 - 10752
  
  base.ind.for[[2]] = 4033:4160 - 3584
  base.ind.rev[[2]] = 12225:12352 - 10752
  
  base.ind.for[[3]] = 3969:4224 - 3584
  base.ind.rev[[3]] = 12161:12416 - 10752
  
  base.ind.for[[4]] = 3841:4352 - 3584
  base.ind.rev[[4]] = 12033:12544 - 10752
  
  base.ind.for[[5]] = 3585:4608 - 3584
  base.ind.rev[[5]] = 11777:12800 - 10752
  
  
  ccount.post.logmean.for.list = list()
  ccount.post.logmean.rev.list = list()
  
  ccount.post.logmean.list = list()
  
  
  for(l in 1:5){
    ##estimating mean intensity
    nr = 10
    eff.for = matrix(0,nr,nbase[l])
    eff.rev = matrix(0,nr,nbase[l])
    base.for = matrix(0,nr,nbase[l])
    base.rev = matrix(0,nr,nbase[l])
    
    for(i in 1:nr){
      ind = sample(1:dim(dprof.train)[1],min(round(dim(dprof)[1]/3),3000))
      est.for = multiseq(as.matrix(dprof.train[ind,base.ind.for[[l]]]),log(ccount.train[ind]+1),lm.approx = TRUE)
      est.rev = multiseq(as.matrix(dprof.train[ind,base.ind.rev[[l]]]),log(ccount.train[ind]+1),lm.approx = TRUE)
      base.for[i,] = est.for$baseline.mean
      base.rev[i,] = est.rev$baseline.mean
      eff.for[i,] = est.for$effect.mean
      eff.rev[i,] = est.rev$effect.mean
      print(i)
    }
    
    base.mean.for = colMeans(base.for)
    base.mean.rev = colMeans(base.rev)
    eff.mean.for = colMeans(eff.for)
    eff.mean.rev = colMeans(eff.rev)
    
    ##estimating prior
    prior.breaks.ini = function(x0,n){
      x = c(-x0,x0)
      for(i in 3:n){
        x[i] = 2 * log(i - 1) - x[i - 1]
      }
      return(x)
    }
    prior.max = ceiling(10 * log(max(ccount) + 1))/10
    prior.breaks = c(prior.breaks.ini(0.445,11),seq(2.5,prior.max,0.1))
    ccount.prior.info = hist(log(ccount.train+1),breaks = prior.breaks, plot = FALSE)
    ccount.prior.val = ccount.prior.info$mids
    ccount.prior.prob = normalize(smooth.spline(ccount.prior.info$mids,ccount.prior.info$counts)$y)
    
    ##computing the likelihood
    lambda.for = exp(rep(1,length(ccount.prior.val))%o%base.mean.for + ccount.prior.val%o%eff.mean.for)
    lambda.rev = exp(rep(1,length(ccount.prior.val))%o%base.mean.rev + ccount.prior.val%o%eff.mean.rev)
    lik.for = matrix(0,length(ccount.test),length(ccount.prior.val))
    lik.rev = matrix(0,length(ccount.test),length(ccount.prior.val))
    for(i in 1:length(ccount.test)){
      loglik.ini.for = rowSums(t(apply(lambda.for,1,dpois,x = as.numeric(dprof.test[i,base.ind.for[[l]]]),log = TRUE)))
      loglik.ini.for = loglik.ini.for - max(loglik.ini.for)
      loglik.ini.rev = rowSums(t(apply(lambda.rev,1,dpois,x = as.numeric(dprof.test[i,base.ind.rev[[l]]]),log = TRUE)))
      loglik.ini.rev = loglik.ini.rev - max(loglik.ini.rev)
      lik.for[i,] = exp(loglik.ini.for)
      lik.rev[i,] = exp(loglik.ini.rev)
    }
    
      ##computing the posterior
      ccount.post.val = ccount.prior.val
      ccount.post.prob.for = lik.for * (rep(1,length(ccount.test))%o%ccount.prior.prob)
      ccount.post.prob.rev = lik.rev * (rep(1,length(ccount.test))%o%ccount.prior.prob)
      ccount.post.prob.for = t(apply(ccount.post.prob.for,1,normalize))
      ccount.post.prob.rev = t(apply(ccount.post.prob.rev,1,normalize))
      
      ccount.post.mean.for = 0  
      ccount.post.mode.for = 0
      ccount.post.logmean.for = 0
      ccount.post.mean.rev = 0
      ccount.post.mode.rev = 0
      ccount.post.logmean.rev = 0
      for(i in 1:length(ccount.test)){
 #       ccount.post.mode.for[i] = exp(ccount.post.val[which(ccount.post.prob.for[i,]==max(ccount.post.prob.for[i,]))])
 #       ccount.post.mean.for[i] = sum(exp(ccount.post.val)*ccount.post.prob.for[i,])
        ccount.post.logmean.for[i] = sum(ccount.post.val*ccount.post.prob.for[i,])
 #      ccount.post.mode.rev[i] = exp(ccount.post.val[which(ccount.post.prob.rev[i,]==max(ccount.post.prob.rev[i,]))])
 #       ccount.post.mean.rev[i] = sum(exp(ccount.post.val)*ccount.post.prob.rev[i,])
        ccount.post.logmean.rev[i] = sum(ccount.post.val*ccount.post.prob.rev[i,])
        
      }
      
      ccount.post.logmean.for.list[[l]] = ccount.post.logmean.for
      ccount.post.logmean.rev.list[[l]] = ccount.post.logmean.rev
    
    ##computing the posterior
    ccount.post.val = ccount.prior.val
    ccount.post.prob = lik.for * lik.rev * (rep(1,length(ccount.test))%o%ccount.prior.prob)
	  zero.ind = rowSums(ccount.post.prob) == 0 #use prior probabilities if likelihood = 0
    ccount.post.prob[zero.ind, ] = matrix(rep(ccount.prior.prob, sum(zero.ind)), byrow=TRUE, nr = sum(zero.ind))
    ccount.post.prob = t(apply(ccount.post.prob,1,normalize))
    
    ccount.post.mean = 0  
    ccount.post.mode = 0
    ccount.post.logmean = 0
    
    for(i in 1:length(ccount.test)){
 #     ccount.post.mode[i] = exp(ccount.post.val[which(ccount.post.prob[i,]==max(ccount.post.prob[i,]))])
 #     ccount.post.mean[i] = sum(exp(ccount.post.val)*ccount.post.prob[i,])
      ccount.post.logmean[i] = sum(ccount.post.val*ccount.post.prob[i,])
      
    }
    
    ccount.post.logmean.list[[l]] = ccount.post.logmean   
  }
  
  setwd("results/roc")
  
#   unbound.ind = ccount.test <= median(ccount.test)
#   bound.ind = ccount.test > quantile(ccount.test, probs = 0.5)

  unbound.ind = unbound.ind.all[!train.ind]
  bound.ind = bound.ind.all[!train.ind]
  
#   pdf(paste("roc", factor.name, "forward_strand_ms.pdf", sep = "_"))
#   roc.res.ms = rocplot(ccount.post.logmean.for.list[[1]][unbound.ind], ccount.post.logmean.for.list[[1]][bound.ind])
#   plot(roc.res.ms$fpr,roc.res.ms$tpr, type = 'l', xlim = c(0, 1), ylim = c(0, 1), xlab = "FPR", ylab = "TPR", main = "Both strands, multiseq & #dnase cuts")
#   
#   roc.res.ms = rocplot(ccount.post.logmean.for.list[[2]][unbound.ind], ccount.post.logmean.for.list[[2]][bound.ind])
#   lines(roc.res.ms$fpr,roc.res.ms$tpr, col = 2)
#   
#   roc.res.ms = rocplot(ccount.post.logmean.for.list[[3]][unbound.ind], ccount.post.logmean.for.list[[3]][bound.ind])
#   lines(roc.res.ms$fpr,roc.res.ms$tpr, col = 3)
#   
#   roc.res.ms = rocplot(ccount.post.logmean.for.list[[4]][unbound.ind], ccount.post.logmean.for.list[[4]][bound.ind])
#   lines(roc.res.ms$fpr,roc.res.ms$tpr, col = 4)
# 
#   legend("bottomright", legend = c("128", "256", "512", "1024"), lty = rep(1, 4), col = 1:4)
#   dev.off()

  dcut = list()
  dcut.all = list()
  for(l in 1:5){
    dcut[[l]] = rowSums(dprof.test[, c(base.ind.for[[l]], base.ind.rev[[l]])])
    dcut.all[[l]] = rowSums(dprof[, c(base.ind.for[[l]], base.ind.rev[[l]])])
  }  

  pscore = cinfo[, 5]

#   centFit = list()
  centFit.pwm = list()
  centFit.pwm.all = list()
  for(l in 1:5){
#     centFit[[l]] <- fitCentipede(Xlist = list(DNase = as.matrix(dprof)[, c(base.ind.for[[l]], base.ind.rev[[l]])]), Y = matrix(rep(1,dim(dprof)[1], nc = 1)), sweeps = 1000)
#     centFit[[l]]$PostPr = centFit[[l]]$PostPr[!train.ind]
    try.cent = try(centFit <- fitCentipede(Xlist = list(DNase = as.matrix(dprof)[, c(base.ind.for[[l]], base.ind.rev[[l]])]), Y = cbind(rep(1,dim(dprof)[1]),pscore), sweeps = 1000))
    if(class(try.cent)=="try-error"){
      centFit.pwm[[l]] = rep(NA, length(ccount.test))
      centFit.pwm.all[[l]] = rep(NA, length(ccount))
    }else{
      centFit.pwm[[l]] = centFit$PostPr[!train.ind]
      centFit.pwm.all[[l]] = centFit$PostPr
    }
  }
 
  auc.ms = 0
  auc.dcut = 0
  auc.cent = 0
  auc.dcut.all = 0
  auc.cent.all = 0
  for(l in 1:5){
    roc.res.ms = roc(controls = ccount.post.logmean.list[[l]][unbound.ind], cases = ccount.post.logmean.list[[l]][bound.ind])
    auc.ms[l] = roc.res.ms$auc
    roc.res.dcut = roc(controls = dcut[[l]][unbound.ind], cases = dcut[[l]][bound.ind])
    auc.dcut[l] = roc.res.dcut$auc
    roc.res.dcut.all = roc(controls = dcut.all[[l]][unbound.ind.all], cases = dcut.all[[l]][bound.ind.all])
    auc.dcut.all[l] = roc.res.dcut.all$auc
    try.cent = try(roc.res.cent <- roc(controls = centFit.pwm[[l]][unbound.ind], cases = centFit.pwm[[l]][bound.ind]))
    if(class(try.cent)=="try-error"){
      auc.cent[l] = NA
    }else{
      auc.cent[l] = roc.res.cent$auc
    }
    try.cent = try(roc.res.cent.all <- roc(controls = centFit.pwm.all[[l]][unbound.ind.all], cases = centFit.pwm.all[[l]][bound.ind.all]))
    if(class(try.cent)=="try-error"){
      auc.cent.all[l] = NA
    }else{
      auc.cent.all[l] = roc.res.cent.all$auc
    }
  }

  pdf(paste(factor.name, "auc", "binary.pdf", sep = "_"))
  plot(6:10, auc.ms, type = "b", ylim = c(0, 1), xlab = "window", ylab = "auc", main = "multiseq vs dcuts, auc")
  lines(6:10, auc.dcut, col = 2, type = "b")
  lines(6:10, auc.cent, col = 4, type = "b")
  legend("bottomright", legend = c("Multiseq", "Dcuts", "Centipede"), lty = c(1, 1, 1), col = c(1,2,4))
  dev.off()

  pdf(paste(factor.name, "auc", "binary_all.pdf", sep = "_"))
  plot(6:10, auc.dcut.all, type = "b", col = 2, ylim = c(0, 1), xlab = "window", ylab = "auc", main = "multiseq vs dcuts, auc, all sites")
  lines(6:10, auc.cent.all, col = 4, type = "b")
  text(x = 6:10, y = auc.dcut.all, labels = round(auc.dcut.all, 3), pos = 3, cex = 0.7, col = 2)
  text(x = 6:10, y = auc.cent.all, labels = round(auc.cent.all, 3), pos = 1, cex = 0.7, col = 4)
  legend("bottomright", legend = c("Dcuts", "Centipede"), lty = c(1, 1), col = c(2, 4))
  dev.off()

  pdf(paste(factor.name, "roc", "binary.pdf", sep = "_"))
  roc.res.ms = rocplot(ccount.post.logmean.list[[1]][unbound.ind], ccount.post.logmean.list[[1]][bound.ind])
  plot(roc.res.ms$fpr,roc.res.ms$tpr, type = 'l', xlim = c(0, 1), ylim = c(0, 1), xlab = "FPR", ylab = "TPR", main = "multiseq and dcuts, bound vs unbound")
  
  for(l in 2:5){
    roc.res.ms = rocplot(ccount.post.logmean.list[[l]][unbound.ind], ccount.post.logmean.list[[l]][bound.ind])
    lines(roc.res.ms$fpr,roc.res.ms$tpr, col = l)
  }
  
  for(l in 1:5){
    roc.res.dcut = rocplot(dcut[[l]][unbound.ind], dcut[[l]][bound.ind])
    lines(roc.res.dcut$fpr,roc.res.dcut$tpr, lty = 2, col = l)
  }

  for(l in 1:5){
    roc.res.cent = rocplot(centFit.pwm[[l]][unbound.ind], centFit.pwm[[l]][bound.ind])
    lines(roc.res.cent$fpr,roc.res.cent$tpr, lty = 6, col = l)
  }

  legend("bottomright", legend = c("64-MS", "128-MS", "256-MS", "512-MS", "1024-MS", "64-Dcuts", "128-Dcuts", "256-Dcuts", "512-Dcuts", "1024-Dcuts", "64-Centipede", "128-Centipede", "256-Centipede", "512-Centipede", "1024-Centipede"), lty = rep(c(1, 2, 6), each = 5), col = rep(1:5, times = 3))
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
  
#roc plots using counts alone
#   unbound.ind.counts = ccount.test <= median(ccount.test)
#   bound.ind.counts = ccount.test > quantile(ccount.test, probs = 0.95)  
#   
#   pdf(paste("roc", factor.name, "binary_counts.pdf", sep = "_"))
#   roc.res.ms = rocplot(ccount.post.logmean.list[[1]][unbound.ind.counts], ccount.post.logmean.list[[1]][bound.ind.counts])
#   plot(roc.res.ms$fpr,roc.res.ms$tpr, type = 'l', xlim = c(0, 1), ylim = c(0, 1), xlab = "FPR", ylab = "TPR", main = "multiseq and dcuts, bound vs unbound using counts only")
#   
#   roc.res.ms = rocplot(ccount.post.logmean.list[[2]][unbound.ind.counts], ccount.post.logmean.list[[2]][bound.ind.counts])
#   lines(roc.res.ms$fpr,roc.res.ms$tpr, col = 2)
#   
#   roc.res.ms = rocplot(ccount.post.logmean.list[[3]][unbound.ind.counts], ccount.post.logmean.list[[3]][bound.ind.counts])
#   lines(roc.res.ms$fpr,roc.res.ms$tpr, col = 3)
#   
#   roc.res.ms = rocplot(ccount.post.logmean.list[[4]][unbound.ind.counts], ccount.post.logmean.list[[4]][bound.ind.counts])
#   lines(roc.res.ms$fpr,roc.res.ms$tpr, col = 4)
#   
#   roc.res.ms = rocplot(ccount.post.logmean.list[[5]][unbound.ind.counts], ccount.post.logmean.list[[5]][bound.ind.counts])
#   lines(roc.res.ms$fpr,roc.res.ms$tpr, col = 5)
#   
#   roc.res.dcut = rocplot(dcut[[1]][unbound.ind.counts], dcut[[1]][bound.ind.counts])
#   lines(roc.res.dcut$fpr,roc.res.dcut$tpr, lty = 2)
#   
#   roc.res.dcut = rocplot(dcut[[2]][unbound.ind.counts], dcut[[2]][bound.ind.counts])
#   lines(roc.res.dcut$fpr,roc.res.dcut$tpr, col = 2, lty = 2)
#   
#   roc.res.dcut = rocplot(dcut[[3]][unbound.ind.counts], dcut[[3]][bound.ind.counts])
#   lines(roc.res.dcut$fpr,roc.res.dcut$tpr, col = 3, lty = 2)
#   
#   roc.res.dcut = rocplot(dcut[[4]][unbound.ind.counts], dcut[[4]][bound.ind.counts])
#   lines(roc.res.dcut$fpr,roc.res.dcut$tpr, col = 4, lty = 2)
#   
#   roc.res.dcut = rocplot(dcut[[5]][unbound.ind.counts], dcut[[5]][bound.ind.counts])
#   lines(roc.res.dcut$fpr,roc.res.dcut$tpr, col = 5, lty = 2)
#   
#   legend("bottomright", legend = c("64-MS", "128-MS", "256-MS", "512-MS", "1024-MS", "64-Dcuts", "128-Dcuts", "256-Dcuts", "512-Dcuts", "1024-Dcuts"), lty = rep(1:2, each = 5), col = rep(1:5, times = 2))
#   dev.off()

    
  corr.ms = NULL
  corr.dcut = NULL
  corr.rank.ms = NULL
  corr.rank.dcut = NULL
  corr.rank.cent = NULL
  
  for(l in 1:5){
    corr.ms[l] = cor(ccount.post.logmean.list[[l]], log(ccount.test + 1))
    corr.dcut[l] = cor(log(dcut[[l]] + 1), log(ccount.test + 1))
    corr.rank.ms[l] = cor(ccount.post.logmean.list[[l]], log(ccount.test + 1), method = "spearman")
    corr.rank.dcut[l] = cor(log(dcut[[l]] + 1), log(ccount.test + 1), method = "spearman")
    corr.rank.cent[l] = cor(centFit.pwm[[l]], log(ccount.test + 1), method = "spearman")
  }


  pdf(paste(factor.name, "corr_pearson.pdf", sep = "_"))
  plot(6:10, corr.ms, type = "b", ylim = c(0, 1), xlab = "window size (log2)", ylab = "correlation", main = "multiseq vs dcuts, correlation")
  lines(6:10, corr.dcut, col = 2, type = "b")
  legend("bottomright", legend = c("Multiseq", "Dcuts"), lty = c(1, 1), col = 1:2)
  dev.off()

  pdf(paste(factor.name, "corr_spearman.pdf", sep = "_"))
  plot(6:10, corr.rank.ms, type = "b", ylim = c(0, 1), xlab = "window size (log2)", ylab = "correlation", main = "multiseq vs dcuts, rank correlation")
  lines(6:10, corr.rank.dcut, col = 2, type = "b")
  lines(6:10, corr.rank.cent, col = 4, type = "b")
  legend("bottomright", legend = c("Multiseq", "Dcuts", "Centipede"), lty = c(1, 1, 1), col = c(1, 2, 4))
  dev.off()

#   unbound.ind = ccount.test <= median(ccount.test)
#   bound.ind = ccount.test > quantile(ccount.test, probs = 0.95)  
#   pdf(paste("roc", factor.name, "forward_0.95.pdf", sep = "_"))
#   roc.res.ms = rocplot(ccount.post.logmean.list[[1]][unbound.ind], ccount.post.logmean.list[[1]][bound.ind])
#   plot(roc.res.ms$fpr,roc.res.ms$tpr, type = 'l', xlim = c(0, 1), ylim = c(0, 1), xlab = "FPR", ylab = "TPR", main = "Forward strand, cutoff = 0.95")
#   
#   roc.res.ms = rocplot(ccount.post.logmean.list[[2]][unbound.ind], ccount.post.logmean.list[[2]][bound.ind])
#   lines(roc.res.ms$fpr,roc.res.ms$tpr, col = 2)
#   
#   roc.res.ms = rocplot(ccount.post.logmean.list[[3]][unbound.ind], ccount.post.logmean.list[[3]][bound.ind])
#   lines(roc.res.ms$fpr,roc.res.ms$tpr, col = 3)
#   
#   roc.res.ms = rocplot(ccount.post.logmean.list[[4]][unbound.ind], ccount.post.logmean.list[[4]][bound.ind])
#   lines(roc.res.ms$fpr,roc.res.ms$tpr, col = 4)
#   
#   legend("bottomright", legend = c("128", "256", "512", "1024"), lty = rep(1, 4), col = 1:4)
#   dev.off()
#   
#   pdf(paste("roc", factor.name, "both_0.95.pdf", sep = "_"))
#   roc.res.ms = rocplot(ccount.post.logmean.for.list[[1]][unbound.ind], ccount.post.logmean.for.list[[1]][bound.ind])
#   plot(roc.res.ms$fpr,roc.res.ms$tpr, type = 'l', xlim = c(0, 1), ylim = c(0, 1), xlab = "FPR", ylab = "TPR", main = "Both strands, cutoff = 0.95")
#   
#   roc.res.ms = rocplot(ccount.post.logmean.for.list[[2]][unbound.ind], ccount.post.logmean.for.list[[2]][bound.ind])
#   lines(roc.res.ms$fpr,roc.res.ms$tpr, col = 2)
#   
#   roc.res.ms = rocplot(ccount.post.logmean.for.list[[3]][unbound.ind], ccount.post.logmean.for.list[[3]][bound.ind])
#   lines(roc.res.ms$fpr,roc.res.ms$tpr, col = 3)
#   
#   roc.res.ms = rocplot(ccount.post.logmean.for.list[[4]][unbound.ind], ccount.post.logmean.for.list[[4]][bound.ind])
#   lines(roc.res.ms$fpr,roc.res.ms$tpr, col = 4)
#   
#   legend("bottomright", legend = c("128", "256", "512", "1024"), lty = rep(1, 4), col = 1:4)
#   dev.off()
#   
#   
#   pdf(paste("roc", factor.name, "strand_comp_0.95.pdf", sep = "_"))
#   roc.res.ms = rocplot(ccount.post.logmean.list[[2]][unbound.ind], ccount.post.logmean.list[[2]][bound.ind])
#   plot(roc.res.ms$fpr,roc.res.ms$tpr, type = 'l', xlim = c(0, 1), ylim = c(0, 1), xlab = "FPR", ylab = "TPR", main = "Strand comparison, cutoff = 0.95")
#   
#   roc.res.ms = rocplot(ccount.post.logmean.for.list[[2]][unbound.ind], ccount.post.logmean.for.list[[2]][bound.ind])
#   lines(roc.res.ms$fpr,roc.res.ms$tpr, col = 2)
#   
#   legend("bottomright", legend = c("Both", "Forward"), lty = rep(1, 2), col = 1:2)
#   dev.off()
#   
#   
#   pdf(paste("roc", factor.name, "dcut_0.95.pdf", sep = "_"))
#   dcut = rowSums(dprof.test[, c(463:562, 1487:1586)])
#   roc.res.dcut = rocplot(dcut[unbound.ind], dcut[bound.ind])
#   
#   plot(roc.res.dcut$fpr,roc.res.dcut$tpr, type='l', xlim = c(0, 1), ylim = c(0, 1), xlab = "FPR", ylab = "TPR", main = "#DNAse cuts, cutoff = 0.95")
#   
#   
#   dcut = rowSums(dprof.test[, c(413:612, 1437:1636)])
#   roc.res.dcut = rocplot(dcut[unbound.ind], dcut[bound.ind])
#   
#   lines(roc.res.dcut$fpr,roc.res.dcut$tpr, col=2)
#   
#   dcut = rowSums(dprof.test[, c(313:712, 1337: 1736)])
#   roc.res.dcut = rocplot(dcut[unbound.ind], dcut[bound.ind])
#   
#   lines(roc.res.dcut$fpr,roc.res.dcut$tpr, col=3)
#   
#   dcut = rowSums(dprof.test)
#   roc.res.dcut = rocplot(dcut[unbound.ind], dcut[bound.ind])
#   
#   lines(roc.res.dcut$fpr,roc.res.dcut$tpr, col=4)
#   
#   legend("bottomright", legend = c("100", "200", "400", "1024"), lty = rep(1, 4), col = 1:4)
#   dev.off()
#     
#   
#   pdf(paste("roc", factor.name, "ms_vs_dcut_0.95.pdf", sep = "_"))
#   roc.res.ms = rocplot(ccount.post.logmean.list[[2]][unbound.ind], ccount.post.logmean.list[[2]][bound.ind])
#   plot(roc.res.ms$fpr,roc.res.ms$tpr, type = 'l', xlim = c(0, 1), ylim = c(0, 1), xlab = "FPR", ylab = "TPR", main = "Multiseq vs #DNase cuts, cutoff = 0.95")
#   
#   dcut = rowSums(dprof.test[, c(413:612, 1437:1636)])
#   roc.res.dcut = rocplot(dcut[unbound.ind], dcut[bound.ind])
#   
#   lines(roc.res.dcut$fpr,roc.res.dcut$tpr, col=2)
#   
#   legend("bottomright", legend = c("Multiseq", "DNase cuts"), lty = rep(1, 2), col = 1:2)
#   dev.off()
  
  
  setwd("../..")

  print(j)
}
  
