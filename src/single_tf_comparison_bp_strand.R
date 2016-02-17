library(multiseq)
library(CENTIPEDE)
library(pROC)

nbase = 2^10

normalize = function(x) x/sum(x)

setwd("~/projects/multiseq_analysis")

source("src/rocplot.R")

file.names=list.files("data/", pattern = "*.txt.gz")
chip.file.names=file.names[seq(4,248,4)]
dnase.file.names=file.names[seq(2,246,4)]

j = 7

  factor.name.short = strsplit(chip.file.names[j],"_")[[1]][1]
  factor.name = paste(strsplit(chip.file.names[j],"_")[[1]][1], strsplit(chip.file.names[j],"_")[[1]][2],sep="_")
  if(as.numeric(file.info(file.path("data",dnase.file.names[j]))[1])>(20*1e+6)) next
  try.test=try(dprof<-read.table(file.path("data",dnase.file.names[j])))
  if(class(try.test)=="try-error") next
  cinfo=read.table(file.path("data",chip.file.names[j]),header=TRUE)
  ccount=cinfo[,6]
  peak.info = read.table(file.path("data","macs_peaks", paste(factor.name.short, "q1pc", "peaks.bed", sep = "_")))
  peak.ind = 0
  nbase = 2^10
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
  

  unbound.ind = unbound.ind.all[!train.ind]
  bound.ind = bound.ind.all[!train.ind]
  
  dcut = list()
  for(l in 1:5){
    dcut[[l]] = rowSums(dprof.test[, c(base.ind.for[[l]], base.ind.rev[[l]])])
  }  
  
  pscore = cinfo[, 5]
  
  #   centFit = list()
  centFit.pwm = list()
  for(l in 1:5){
    #     centFit[[l]] <- fitCentipede(Xlist = list(DNase = as.matrix(dprof)[, c(base.ind.for[[l]], base.ind.rev[[l]])]), Y = matrix(rep(1,dim(dprof)[1], nc = 1)), sweeps = 1000)
    #     centFit[[l]]$PostPr = centFit[[l]]$PostPr[!train.ind]
    try.cent = try(centFit <- fitCentipede(Xlist = list(DNase = as.matrix(dprof)[, c(base.ind.for[[l]], base.ind.rev[[l]])]), Y = cbind(rep(1,dim(dprof)[1]),pscore), sweeps = 1000))
    if(class(try.cent)=="try-error"){
      centFit.pwm[[l]] = rep(NULL, length(ccount.test))
    }else{
      centFit.pwm[[l]] = centFit$PostPr[!train.ind]
    }
  }
  
  auc.ms = 0
  auc.dcut = 0
  auc.cent = 0
  for(l in 1:5){
    roc.res.ms = roc(controls = ccount.post.logmean.list[[l]][unbound.ind], cases = ccount.post.logmean.list[[l]][bound.ind])
    auc.ms[l] = roc.res.ms$auc
    roc.res.dcut = roc(controls = dcut[[l]][unbound.ind], cases = dcut[[l]][bound.ind])
    auc.dcut[l] = roc.res.dcut$auc
    try.cent = try(roc.res.cent <- roc(controls = centFit.pwm[[l]][unbound.ind], cases = centFit.pwm[[l]][bound.ind]))
    if(class(try.cent)=="try-error"){
      auc.cent[l] = NA
    }else{
      auc.cent[l] = roc.res.cent$auc
    }
  }
  
  pdf(paste(factor.name, "auc", "binary.pdf", sep = "_"))
  plot(6:10, auc.ms, type = "b", ylim = c(0, 1), xlab = "window", ylab = "auc", main = "multiseq vs dcuts, auc")
  lines(6:10, auc.dcut, col = 2, type = "b")
  lines(6:10, auc.cent, col = 4, type = "b")
  legend("bottomright", legend = c("Multiseq", "Dcuts", "Centipede"), lty = c(1, 1, 1), col = 1:3)
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
  legend("bottomright", legend = c("Multiseq", "Dcuts", "Centipede"), lty = c(1, 1, 1), col = 1:3)
  dev.off()
 
  
  setwd("../..")
  


