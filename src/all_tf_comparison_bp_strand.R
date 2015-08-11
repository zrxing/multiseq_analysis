library(multiseq)

nbase = 2^10

normalize = function(x) x/sum(x)

setwd("~/projects/multiseq_analysis")

source("src/rocplot.R")

file.names=list.files("data/")
chip.file.names=file.names[seq(3,247,4)]
dnase.file.names=file.names[seq(2,246,4)]

for(j in 1:length(chip.file.names)){
  factor.name=paste(strsplit(chip.file.names[j],"_")[[1]][1],strsplit(chip.file.names[j],"_")[[1]][2],sep="_")
  try.test=try(dprof<-read.table(file.path("data",dnase.file.names[j])))
  if(class(try.test)=="try-error") next
  cinfo=read.table(file.path("data",chip.file.names[j]),header=TRUE)
  ccount=cinfo[,6]
  
  train.size = round(0.8*dim(dprof)[1])
  set.seed(417)
  train.ind = sample(1:dim(dprof)[1],train.size,replace = FALSE)
  train.ind = (1:dim(dprof)[1])%in%train.ind
  
  dprof.train = dprof[train.ind,]
  dprof.test = dprof[!train.ind,]
  
  ccount.train = ccount[train.ind]
  ccount.test = ccount[!train.ind]
  
  
  nbase = 2^7
  nbase[2] = 2^8
  nbase[3] = 2^9
  nbase[4] = 2^10
  
  base.ind.for = list()
  base.ind.rev = list()
  
  base.ind.for[[1]] = 4033:4160 - 3584
  base.ind.rev[[1]] = 12225:12352 - 10752
  
  base.ind.for[[2]] = 3969:4224 - 3584
  base.ind.rev[[2]] = 12161:12416 - 10752
  
  base.ind.for[[3]] = 3841:4352 - 3584
  base.ind.rev[[3]] = 12033:12544 - 10752
  
  base.ind.for[[4]] = 3585:4608 - 3584
  base.ind.rev[[4]] = 11777:12800 - 10752
  
  
  ccount.post.logmean.for.list = list()
  ccount.post.logmean.rev.list = list()
  
  ccount.post.logmean.list = list()
  
  
  for(l in 1:4){
    ##estimating mean intensity
    nr = 10
    eff.for = matrix(0,nr,nbase[l])
    eff.rev = matrix(0,nr,nbase[l])
    base.for = matrix(0,nr,nbase[l])
    base.rev = matrix(0,nr,nbase[l])
    
    for(i in 1:nr){
      ind = sample(1:dim(dprof.train)[1],3000)
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
        ccount.post.mode.for[i] = exp(ccount.post.val[which(ccount.post.prob.for[i,]==max(ccount.post.prob.for[i,]))])
        ccount.post.mean.for[i] = sum(exp(ccount.post.val)*ccount.post.prob.for[i,])
        ccount.post.logmean.for[i] = sum(ccount.post.val*ccount.post.prob.for[i,])
        ccount.post.mode.rev[i] = exp(ccount.post.val[which(ccount.post.prob.rev[i,]==max(ccount.post.prob.rev[i,]))])
        ccount.post.mean.rev[i] = sum(exp(ccount.post.val)*ccount.post.prob.rev[i,])
        ccount.post.logmean.rev[i] = sum(ccount.post.val*ccount.post.prob.rev[i,])
        
      }
      
      ccount.post.logmean.for.list[[l]] = ccount.post.logmean.for
      ccount.post.logmean.rev.list[[l]] = ccount.post.logmean.rev
    
    ##computing the posterior
    ccount.post.val = ccount.prior.val
    ccount.post.prob = lik.for * lik.rev * (rep(1,length(ccount.test))%o%ccount.prior.prob)
    ccount.post.prob = t(apply(ccount.post.prob,1,normalize))
    
    ccount.post.mean = 0  
    ccount.post.mode = 0
    ccount.post.logmean = 0
    
    for(i in 1:length(ccount.test)){
      ccount.post.mode[i] = exp(ccount.post.val[which(ccount.post.prob[i,]==max(ccount.post.prob[i,]))])
      ccount.post.mean[i] = sum(exp(ccount.post.val)*ccount.post.prob[i,])
      ccount.post.logmean[i] = sum(ccount.post.val*ccount.post.prob[i,])
      
    }
    
    ccount.post.logmean.list[[l]] = ccount.post.logmean   
  }
  
  setwd("results/roc")
  
  unbound.ind = ccount.test <= median(ccount.test)
  bound.ind = ccount.test > quantile(ccount.test, probs = 0.5)
  
  
  pdf(paste("roc", factor.name, "forward_0.5.pdf", sep = "_"))
  roc.res.ms = rocplot(ccount.post.logmean.list[[1]][unbound.ind], ccount.post.logmean.list[[1]][bound.ind])
  plot(roc.res.ms$fpr,roc.res.ms$tpr, type = 'l', xlim = c(0, 1), ylim = c(0, 1), xlab = "FPR", ylab = "TPR", main = "Forward strand, cutoff = 0.5")
  
  roc.res.ms = rocplot(ccount.post.logmean.list[[2]][unbound.ind], ccount.post.logmean.list[[2]][bound.ind])
  lines(roc.res.ms$fpr,roc.res.ms$tpr, col = 2)
  
  roc.res.ms = rocplot(ccount.post.logmean.list[[3]][unbound.ind], ccount.post.logmean.list[[3]][bound.ind])
  lines(roc.res.ms$fpr,roc.res.ms$tpr, col = 3)
  
  roc.res.ms = rocplot(ccount.post.logmean.list[[4]][unbound.ind], ccount.post.logmean.list[[4]][bound.ind])
  lines(roc.res.ms$fpr,roc.res.ms$tpr, col = 4)
  
  legend("bottomright", legend = c("128", "256", "512", "1024"), lty = rep(1, 4), col = 1:4)
  dev.off()
  
  pdf(paste("roc", factor.name, "both_0.5.pdf", sep = "_"))
  roc.res.ms = rocplot(ccount.post.logmean.for.list[[1]][unbound.ind], ccount.post.logmean.for.list[[1]][bound.ind])
  plot(roc.res.ms$fpr,roc.res.ms$tpr, type = 'l', xlim = c(0, 1), ylim = c(0, 1), xlab = "FPR", ylab = "TPR", main = "Both strands, cutoff = 0.5")
  
  roc.res.ms = rocplot(ccount.post.logmean.for.list[[2]][unbound.ind], ccount.post.logmean.for.list[[2]][bound.ind])
  lines(roc.res.ms$fpr,roc.res.ms$tpr, col = 2)
  
  roc.res.ms = rocplot(ccount.post.logmean.for.list[[3]][unbound.ind], ccount.post.logmean.for.list[[3]][bound.ind])
  lines(roc.res.ms$fpr,roc.res.ms$tpr, col = 3)
  
  roc.res.ms = rocplot(ccount.post.logmean.for.list[[4]][unbound.ind], ccount.post.logmean.for.list[[4]][bound.ind])
  lines(roc.res.ms$fpr,roc.res.ms$tpr, col = 4)

  legend("bottomright", legend = c("128", "256", "512", "1024"), lty = rep(1, 4), col = 1:4)
  dev.off()
  
  
  pdf(paste("roc", factor.name, "strand_comp_0.5.pdf", sep = "_"))
  roc.res.ms = rocplot(ccount.post.logmean.list[[2]][unbound.ind], ccount.post.logmean.list[[2]][bound.ind])
  plot(roc.res.ms$fpr,roc.res.ms$tpr, type = 'l', xlim = c(0, 1), ylim = c(0, 1), xlab = "FPR", ylab = "TPR", main = "Strand comparison, cutoff = 0.5")
  
  roc.res.ms = rocplot(ccount.post.logmean.for.list[[2]][unbound.ind], ccount.post.logmean.for.list[[2]][bound.ind])
  lines(roc.res.ms$fpr,roc.res.ms$tpr, col = 2)

  legend("bottomright", legend = c("Forward", "Both"), lty = rep(1, 2), col = 1:2)
  dev.off()
  
  
  pdf(paste("roc", factor.name, "dcut_0.5.pdf", sep = "_"))
  dcut = rowSums(dprof.test[, c(463:562, 1487:1586)])
  roc.res.dcut = rocplot(dcut[unbound.ind], dcut[bound.ind])
  
  plot(roc.res.dcut$fpr,roc.res.dcut$tpr, type='l', xlim = c(0, 1), ylim = c(0, 1), xlab = "FPR", ylab = "TPR", main = "#DNAse cuts, cutoff = 0.5")
  
  
  dcut = rowSums(dprof.test[, c(413:612, 1437:1636)])
  roc.res.dcut = rocplot(dcut[unbound.ind], dcut[bound.ind])
  
  lines(roc.res.dcut$fpr,roc.res.dcut$tpr, col=2)
  
  dcut = rowSums(dprof.test[, c(313:712, 1337: 1736)])
  roc.res.dcut = rocplot(dcut[unbound.ind], dcut[bound.ind])
  
  lines(roc.res.dcut$fpr,roc.res.dcut$tpr, col=3)
  
  dcut = rowSums(dprof.test)
  roc.res.dcut = rocplot(dcut[unbound.ind], dcut[bound.ind])
  
  lines(roc.res.dcut$fpr,roc.res.dcut$tpr, col=4)

  legend("bottomright", legend = c("100", "200", "400", "1024"), lty = rep(1, 4), col = 1:4)
  dev.off()
  
  
  pdf(paste("roc", factor.name, "ms_vs_dcut_0.5.pdf", sep = "_"))
  roc.res.ms = rocplot(ccount.post.logmean.list[[2]][unbound.ind], ccount.post.logmean.list[[2]][bound.ind])
  plot(roc.res.ms$fpr,roc.res.ms$tpr, type = 'l', xlim = c(0, 1), ylim = c(0, 1), xlab = "FPR", ylab = "TPR", main = "Multiseq vs #DNase cuts, cutoff = 0.5")
  
  dcut = rowSums(dprof.test[, c(413:612, 1437:1636)])
  roc.res.dcut = rocplot(dcut[unbound.ind], dcut[bound.ind])
  
  lines(roc.res.dcut$fpr,roc.res.dcut$tpr, col=2)
  
  legend("bottomright", legend = c("Multiseq", "DNase cuts"), lty = rep(1, 2), col = 1:2)
  dev.off()
  
  unbound.ind = ccount.test <= median(ccount.test)
  bound.ind = ccount.test > quantile(ccount.test, probs = 0.95)
  
  
  pdf(paste("roc", factor.name, "forward_0.95.pdf", sep = "_"))
  roc.res.ms = rocplot(ccount.post.logmean.list[[1]][unbound.ind], ccount.post.logmean.list[[1]][bound.ind])
  plot(roc.res.ms$fpr,roc.res.ms$tpr, type = 'l', xlim = c(0, 1), ylim = c(0, 1), xlab = "FPR", ylab = "TPR", main = "Forward strand, cutoff = 0.95")
  
  roc.res.ms = rocplot(ccount.post.logmean.list[[2]][unbound.ind], ccount.post.logmean.list[[2]][bound.ind])
  lines(roc.res.ms$fpr,roc.res.ms$tpr, col = 2)
  
  roc.res.ms = rocplot(ccount.post.logmean.list[[3]][unbound.ind], ccount.post.logmean.list[[3]][bound.ind])
  lines(roc.res.ms$fpr,roc.res.ms$tpr, col = 3)
  
  roc.res.ms = rocplot(ccount.post.logmean.list[[4]][unbound.ind], ccount.post.logmean.list[[4]][bound.ind])
  lines(roc.res.ms$fpr,roc.res.ms$tpr, col = 4)
  
  legend("bottomright", legend = c("128", "256", "512", "1024"), lty = rep(1, 4), col = 1:4)
  dev.off()
  
  pdf(paste("roc", factor.name, "both_0.95.pdf", sep = "_"))
  roc.res.ms = rocplot(ccount.post.logmean.for.list[[1]][unbound.ind], ccount.post.logmean.for.list[[1]][bound.ind])
  plot(roc.res.ms$fpr,roc.res.ms$tpr, type = 'l', xlim = c(0, 1), ylim = c(0, 1), xlab = "FPR", ylab = "TPR", main = "Both strands, cutoff = 0.95")
  
  roc.res.ms = rocplot(ccount.post.logmean.for.list[[2]][unbound.ind], ccount.post.logmean.for.list[[2]][bound.ind])
  lines(roc.res.ms$fpr,roc.res.ms$tpr, col = 2)
  
  roc.res.ms = rocplot(ccount.post.logmean.for.list[[3]][unbound.ind], ccount.post.logmean.for.list[[3]][bound.ind])
  lines(roc.res.ms$fpr,roc.res.ms$tpr, col = 3)
  
  roc.res.ms = rocplot(ccount.post.logmean.for.list[[4]][unbound.ind], ccount.post.logmean.for.list[[4]][bound.ind])
  lines(roc.res.ms$fpr,roc.res.ms$tpr, col = 4)
  
  legend("bottomright", legend = c("128", "256", "512", "1024"), lty = rep(1, 4), col = 1:4)
  dev.off()
  
  
  pdf(paste("roc", factor.name, "strand_comp_0.95.pdf", sep = "_"))
  roc.res.ms = rocplot(ccount.post.logmean.list[[2]][unbound.ind], ccount.post.logmean.list[[2]][bound.ind])
  plot(roc.res.ms$fpr,roc.res.ms$tpr, type = 'l', xlim = c(0, 1), ylim = c(0, 1), xlab = "FPR", ylab = "TPR", main = "Strand comparison, cutoff = 0.95")
  
  roc.res.ms = rocplot(ccount.post.logmean.for.list[[2]][unbound.ind], ccount.post.logmean.for.list[[2]][bound.ind])
  lines(roc.res.ms$fpr,roc.res.ms$tpr, col = 2)
  
  legend("bottomright", legend = c("Forward", "Both"), lty = rep(1, 2), col = 1:2)
  dev.off()
  
  
  pdf(paste("roc", factor.name, "dcut_0.95.pdf", sep = "_"))
  dcut = rowSums(dprof.test[, c(463:562, 1487:1586)])
  roc.res.dcut = rocplot(dcut[unbound.ind], dcut[bound.ind])
  
  plot(roc.res.dcut$fpr,roc.res.dcut$tpr, type='l', xlim = c(0, 1), ylim = c(0, 1), xlab = "FPR", ylab = "TPR", main = "#DNAse cuts, cutoff = 0.95")
  
  
  dcut = rowSums(dprof.test[, c(413:612, 1437:1636)])
  roc.res.dcut = rocplot(dcut[unbound.ind], dcut[bound.ind])
  
  lines(roc.res.dcut$fpr,roc.res.dcut$tpr, col=2)
  
  dcut = rowSums(dprof.test[, c(313:712, 1337: 1736)])
  roc.res.dcut = rocplot(dcut[unbound.ind], dcut[bound.ind])
  
  lines(roc.res.dcut$fpr,roc.res.dcut$tpr, col=3)
  
  dcut = rowSums(dprof.test)
  roc.res.dcut = rocplot(dcut[unbound.ind], dcut[bound.ind])
  
  lines(roc.res.dcut$fpr,roc.res.dcut$tpr, col=4)
  
  legend("bottomright", legend = c("100", "200", "400", "1024"), lty = rep(1, 4), col = 1:4)
  dev.off()
    
  
  pdf(paste("roc", factor.name, "ms_vs_dcut_0.95.pdf", sep = "_"))
  roc.res.ms = rocplot(ccount.post.logmean.list[[2]][unbound.ind], ccount.post.logmean.list[[2]][bound.ind])
  plot(roc.res.ms$fpr,roc.res.ms$tpr, type = 'l', xlim = c(0, 1), ylim = c(0, 1), xlab = "FPR", ylab = "TPR", main = "Multiseq vs #DNase cuts, cutoff = 0.95")
  
  dcut = rowSums(dprof.test[, c(413:612, 1437:1636)])
  roc.res.dcut = rocplot(dcut[unbound.ind], dcut[bound.ind])
  
  lines(roc.res.dcut$fpr,roc.res.dcut$tpr, col=2)
  
  legend("bottomright", legend = c("Multiseq", "DNase cuts"), lty = rep(1, 2), col = 1:2)
  dev.off()
  
  
  setwd("../..")
  
}
  