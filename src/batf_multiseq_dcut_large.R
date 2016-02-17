library(multiseq)
library(pROC)

normalize = function(x) x/sum(x)

setwd("~/projects/multiseq_analysis")

source("src/rocplot.R")


reg.size = 2^14


dprof = read.table("data/Batf_S356_dnase_data_large.txt.gz")
cinfo = read.table("data/Batf_S356_site_and_chip_data.txt.gz",header = TRUE)


ccount=cinfo[,6]
peak.info = read.table("data/macs_peaks/Batf_q1pc_peaks.bed")
peak.ind = 0
for(i in 1:dim(cinfo)[1]){
  peak.ind[i] = sum((as.character(peak.info[, 1]) == as.character(cinfo[i, 1])) & (peak.info[, 2] < cinfo[i, 2]) & (peak.info[, 3] > cinfo[i, 3])) > 0
}

#unbound sites are defined as motif instances that do not lie within a chipseq peak & has less chipseq reads than the overall median
unbound.ind.all = !(peak.ind) & cinfo[, 6] < median(cinfo[, 6])
bound.ind.all = as.logical(peak.ind)

# print("Read data")

train.size = round(0.8*dim(dprof)[1])
set.seed(923)
train.ind = sample(1:dim(dprof)[1],train.size,replace = FALSE)
train.ind = (1:dim(dprof)[1])%in%train.ind

dprof.train = dprof[train.ind,]
dprof.test = dprof[!train.ind,]

ccount.train = ccount[train.ind]
ccount.test = ccount[!train.ind]

nbase = 2^(6:14)

base.ind.for = list()
base.ind.rev = list()

for(l in 1:9){
  base.ind.for[[l]] = (reg.size/2 - nbase[l]/2 + 1):(reg.size/2 + nbase[l]/2)
  base.ind.rev[[l]] = reg.size + (reg.size/2 - nbase[l]/2 + 1):(reg.size/2 + nbase[l]/2)
}

ccount.post.logmean.list = list()


for(l in 1:9){
  ##estimating mean intensity
  nr = 20
  eff.for = matrix(0,nr,nbase[l])
  eff.rev = matrix(0,nr,nbase[l])
  base.for = matrix(0,nr,nbase[l])
  base.rev = matrix(0,nr,nbase[l])
  
  for(i in 1:nr){
    ind = sample(1:dim(dprof.train)[1],min(round(dim(dprof.train)[1]/3),1000))
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

setwd("results/roc/batf_large")

unbound.ind = unbound.ind.all[!train.ind]
bound.ind = bound.ind.all[!train.ind]

dcut = list()
for(l in 1:9){
  dcut[[l]] = rowSums(dprof.test[, c(base.ind.for[[l]], base.ind.rev[[l]])])
}  

pdf(paste("batf", "roc", "binary.pdf", sep = "_"))
roc.res.ms = rocplot(ccount.post.logmean.list[[1]][unbound.ind], ccount.post.logmean.list[[1]][bound.ind])
plot(roc.res.ms$fpr,roc.res.ms$tpr, type = 'l', xlim = c(0, 1), ylim = c(0, 1), xlab = "FPR", ylab = "TPR", main = "multiseq and dcuts, bound vs unbound")

for(l in 2:9){
  roc.res.ms = rocplot(ccount.post.logmean.list[[l]][unbound.ind], ccount.post.logmean.list[[l]][bound.ind])
  lines(roc.res.ms$fpr,roc.res.ms$tpr, col = l)
}

for(l in 1:9){
  roc.res.dcut = rocplot(dcut[[l]][unbound.ind], dcut[[l]][bound.ind])
  lines(roc.res.dcut$fpr,roc.res.dcut$tpr, lty = 2, col = l)
}

legend("bottomright", legend = c("64-MS", "128-MS", "256-MS", "512-MS", "1024-MS", "2048-MS", "4096-MS", "8192-MS", "16384-MS", "64-Dcuts", "128-Dcuts", "256-Dcuts", "512-Dcuts", "1024-Dcuts", "2048-Dcuts", "4096-Dcuts", "8192-Dcuts", "16384-Dcuts"), lty = rep(1:2, each = 9), col = rep(1:9, times = 2))
dev.off()

auc.ms = 0
auc.dcut = 0
for(l in 1:9){
  roc.res.ms = roc(controls = ccount.post.logmean.list[[l]][unbound.ind], cases = ccount.post.logmean.list[[l]][bound.ind])
  auc.ms[l] = roc.res.ms$auc
  roc.res.dcut = roc(controls = dcut[[l]][unbound.ind], cases = dcut[[l]][bound.ind])
  auc.dcut[l] = roc.res.dcut$auc
}

pdf(paste("batf", "auc", "binary.pdf", sep = "_"))
plot(6:14, auc.ms, type = "b", ylim = c(0.8, 1), xlab = "window", ylab = "auc", main = "multiseq vs dcuts, auc")
lines(6:14, auc.dcut, col = 2, type = "b")
legend("bottomright", legend = c("Multiseq", "Dcuts"), lty = c(1, 1), col = 1:2)
dev.off()


corr.ms = NULL
corr.dcut = NULL
corr.rank.ms = NULL
corr.rank.dcut = NULL

for(l in 1:9){
  corr.ms[l] = cor(ccount.post.logmean.list[[l]], log(ccount.test + 1))
  corr.dcut[l] = cor(log(dcut[[l]] + 1), log(ccount.test + 1))
  corr.rank.ms[l] = cor(ccount.post.logmean.list[[l]], log(ccount.test + 1), method = "spearman")
  corr.rank.dcut[l] = cor(log(dcut[[l]] + 1), log(ccount.test + 1), method = "spearman")
}

pdf(paste("batf", "corr_pearson.pdf", sep = "_"))
plot(6:14, corr.ms, type = "b", ylim = c(0, 1), xlab = "window size (log2)", ylab = "correlation", main = "multiseq vs dcuts, correlation")
lines(6:14, corr.dcut, col = 2, type = "b")
legend("bottomright", legend = c("Multiseq", "Dcuts"), lty = c(1, 1), col = 1:2)
dev.off()

pdf(paste("batf", "corr_spearman.pdf", sep = "_"))
plot(6:14, corr.rank.ms, type = "b", ylim = c(0, 1), xlab = "window size (log2)", ylab = "correlation", main = "multiseq vs dcuts, rank correlation")
lines(6:14, corr.rank.dcut, col = 2, type = "b")
legend("bottomright", legend = c("Multiseq", "Dcuts"), lty = c(1, 1), col = 1:2)
dev.off()

save.image("batf_large.RData")