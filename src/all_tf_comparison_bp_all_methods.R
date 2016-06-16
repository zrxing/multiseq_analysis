library(multiseq)
library(CENTIPEDE)
library(pROC)


args = commandArgs(TRUE)

j = as.integer(args[1])


reg.size = 2^10

normalize = function(x) x/sum(x)

setwd("~/projects/multiseq_analysis")

source("src/rocplot.R")



file.names=list.files("data/", pattern = "*.txt.gz")
chip.file.names=file.names[seq(5,310,5)]
dnase.file.names=file.names[seq(2,307,5)]
map.file.names=file.names[seq(3,308,5)]

mscent.data.path = "~/projects/multiseq_analysis/data/mscentipede/"
mscent.data.dest.path = "~/projects/multiseq_analysis/results/roc/mscentipede/"




factor.name.short = strsplit(chip.file.names[j],"_")[[1]][1]
factor.name = paste(strsplit(chip.file.names[j],"_")[[1]][1], strsplit(chip.file.names[j],"_")[[1]][2],sep="_")
if(as.numeric(file.info(file.path("data",dnase.file.names[j]))[1])>(12*1e+6)) stop("Data is too large")
try.test = try(dprof<-read.table(file.path("data",dnase.file.names[j])))
if(class(try.test)=="try-error") stop("DNase data is incomplete")
cinfo = read.table(file.path("data", chip.file.names[j]), header = TRUE)
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
print(factor.name)
print(j)

train.size = round(0.8*dim(dprof)[1])
set.seed(1024)
train.ind = sample(1:dim(dprof)[1],train.size,replace = FALSE)
train.ind = (1:dim(dprof)[1])%in%train.ind

dprof.train = dprof[train.ind,]
dprof.test = dprof[!train.ind,]

ccount.train = ccount[train.ind]
ccount.test = ccount[!train.ind]

nbase = 2^(6:10)

base.ind.for = list()
base.ind.rev = list()

for(l in 1:5){
  base.ind.for[[l]] = (reg.size/2 - nbase[l]/2 + 1):(reg.size/2 + nbase[l]/2)
  base.ind.rev[[l]] = reg.size + (reg.size/2 - nbase[l]/2 + 1):(reg.size/2 + nbase[l]/2)
}

ccount.post.logmean.list = list()

##run multiseq
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


mscentFit = list()
mscentFit.all = list()

for(l in 1:5){
  mscent.model = paste0(mscent.data.dest.path, paste0(factor.name, "_", nbase[l], "_msCentipede_model_parameters.pkl"))
  mscent.binding = paste0(mscent.data.dest.path, paste0(factor.name, "_", nbase[l], "_msCentipede_binding_posterior.txt.gz"))
  mscent.log = paste0(mscent.data.dest.path, paste0(factor.name, "_", nbase[l], "_msCentipede_log.txt"))
  
  mscent.learn = paste0(mscent.data.path, paste0(factor.name, "_data_learn.txt.gz"))
  if(l != 5){
    sys.mscent.learn = paste0("/data/tools/Python2_latest/bin/python2.7 ~/mscentipede/msCentipede/call_binding.py --task learn --mintol 1e-4 --model_file ", mscent.model, " --log_file ", mscent.log, " --window ", nbase[l], " ", mscent.learn, " ", mscent.data.path, "wgEncodeUwDnaseGm12878AlnRep1.sort.bam ", mscent.data.path, "wgEncodeUwDnaseGm12878AlnRep2.sort.bam")
  }else{
    sys.mscent.learn = paste0("/data/tools/Python2_latest/bin/python2.7 ~/mscentipede/msCentipede/call_binding.py --task learn --mintol 1e-4 --batch 7500 --model_file ", mscent.model, " --log_file ", mscent.log, " --window ", nbase[l], " ", mscent.learn, " ", mscent.data.path, "wgEncodeUwDnaseGm12878AlnRep1.sort.bam ", mscent.data.path, "wgEncodeUwDnaseGm12878AlnRep2.sort.bam")
  }
  system(sys.mscent.learn)
  
  mscent.infer = paste0(mscent.data.path, paste0(factor.name, "_data_infer.txt.gz"))
  sys.mscent.infer = paste0("/data/tools/Python2_latest/bin/python2.7 ~/mscentipede/msCentipede/call_binding.py --task infer --mintol 1e-4 --model_file ", mscent.model, " --posterior_file ", mscent.binding, " --log_file ", mscent.log, " --window ", nbase[l], " ", mscent.infer, " ", mscent.data.path, "wgEncodeUwDnaseGm12878AlnRep1.sort.bam ", mscent.data.path, "wgEncodeUwDnaseGm12878AlnRep2.sort.bam")
  system(sys.mscent.infer)
  
  mscent.post = read.table(mscent.binding, header = TRUE)
  
  mscentFit.all[[l]] = mscent.post$LogPosOdds
  mscentFit[[l]] = mscentFit.all[[l]][!train.ind]
}


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

unbound.ind = unbound.ind.all[!train.ind]
bound.ind = bound.ind.all[!train.ind]



auc.ms = 0
auc.dcut = 0
auc.cent = 0
auc.mscent = 0
auc.dcut.all = 0
auc.cent.all = 0
auc.mscent.all = 0
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
  roc.res.mscent = roc(controls = mscentFit[[l]][unbound.ind], cases = mscentFit[[l]][bound.ind])
  auc.mscent[l] = roc.res.mscent$auc
  roc.res.mscent.all = roc(controls = mscentFit.all[[l]][unbound.ind.all], cases = mscentFit.all[[l]][bound.ind.all])
  auc.mscent.all[l] = roc.res.mscent.all$auc
}

setwd("results/roc")

pdf(paste(factor.name, "auc", "binary.pdf", sep = "_"))
plot(6:10, auc.ms, type = "b", ylim = c(0, 1), xlab = "window", ylab = "auc", main = "test motifs, auroc")
lines(6:10, auc.dcut, col = 2, type = "b")
lines(6:10, auc.cent, col = 3, type = "b")
lines(6:10, auc.mscent, col = 4, type = "b")
legend("bottomright", legend = c("Multiseq", "Dcuts", "Centipede", "msCentipede"), lty = c(1, 1, 1, 1), col = 1:4)
dev.off()

pdf(paste(factor.name, "auc", "binary_all.pdf", sep = "_"))
plot(6:10, auc.dcut.all, type = "b", col = 2, xlim = c(5.7, 10), ylim = c(0, 1), xlab = "window", ylab = "auc", main = "all motifs, auroc")
lines(6:10, auc.cent.all, col = 3, type = "b")
lines(6:10, auc.mscent.all, col = 4, type = "b")
text(x = 6:10, y = auc.dcut.all, labels = round(auc.dcut.all, 3), pos = 3, cex = 0.7, col = 2)
text(x = 6:10, y = auc.cent.all, labels = round(auc.cent.all, 3), pos = 1, cex = 0.7, col = 3)
text(x = 6:10, y = auc.mscent.all, labels = round(auc.mscent.all, 3), pos = 2, cex = 0.7, col = 4)
legend("bottomright", legend = c("Dcuts", "Centipede", "msCentipede"), lty = c(1, 1, 1), col = 2:4)
dev.off()

pdf(paste(factor.name, "roc", "binary.pdf", sep = "_"))
roc.res.ms = rocplot(ccount.post.logmean.list[[1]][unbound.ind], ccount.post.logmean.list[[1]][bound.ind])
plot(roc.res.ms$fpr,roc.res.ms$tpr, type = 'l', xlim = c(0, 1), ylim = c(0, 1), xlab = "FPR", ylab = "TPR", main = "test motifs, bound vs unbound roc")

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

for(l in 1:5){
  roc.res.mscent = rocplot(mscentFit[[l]][unbound.ind], mscentFit[[l]][bound.ind])
  lines(roc.res.mscent$fpr,roc.res.mscent$tpr, lty = 3, col = l)
}

legend("bottomright", legend = c("64-MS", "128-MS", "256-MS", "512-MS", "1024-MS", "64-Dcuts", "128-Dcuts", "256-Dcuts", "512-Dcuts", "1024-Dcuts", "64-Centipede", "128-Centipede", "256-Centipede", "512-Centipede", "1024-Centipede", "64-msCentipede", "128-msCentipede", "256-msCentipede", "512-msCentipede", "1024-msCentipede"), lty = rep(c(1, 2, 6, 3), each = 5), col = rep(1:5, times = 4))
dev.off()



corr.ms = NULL
corr.dcut = NULL
corr.rank.ms = NULL
corr.rank.dcut = NULL
corr.rank.cent = NULL
corr.rank.mscent = NULL

corr.rank.dcut.all = NULL
corr.rank.cent.all = NULL
corr.rank.mscent.all = NULL


for(l in 1:5){
  corr.ms[l] = cor(ccount.post.logmean.list[[l]], log(ccount.test + 1))
  corr.dcut[l] = cor(log(dcut[[l]] + 1), log(ccount.test + 1))
  corr.rank.ms[l] = cor(ccount.post.logmean.list[[l]], log(ccount.test + 1), method = "spearman")
  corr.rank.dcut[l] = cor(log(dcut[[l]] + 1), log(ccount.test + 1), method = "spearman")
  corr.rank.cent[l] = cor(centFit.pwm[[l]], log(ccount.test + 1), method = "spearman")
  corr.rank.mscent[l] = cor(mscentFit[[l]], log(ccount.test + 1), method = "spearman")
  
  corr.rank.dcut.all[l] = cor(log(dcut.all[[l]] + 1), log(ccount + 1), method = "spearman")
  corr.rank.cent.all[l] = cor(centFit.pwm.all[[l]], log(ccount + 1), method = "spearman")
  corr.rank.mscent.all[l] = cor(mscentFit.all[[l]], log(ccount + 1), method = "spearman")
}


pdf(paste(factor.name, "corr_pearson.pdf", sep = "_"))
plot(6:10, corr.ms, type = "b", ylim = c(0, 1), xlab = "window size (log2)", ylab = "correlation", main = "test motifs, pearson correlation")
lines(6:10, corr.dcut, col = 2, type = "b")
legend("bottomright", legend = c("Multiseq", "Dcuts"), lty = c(1, 1), col = 1:2)
dev.off()

pdf(paste(factor.name, "corr_spearman.pdf", sep = "_"))
plot(6:10, corr.rank.ms, type = "b", ylim = c(0, 1), xlab = "window size (log2)", ylab = "correlation", main = "test motifs, rank correlation")
lines(6:10, corr.rank.dcut, col = 2, type = "b")
lines(6:10, corr.rank.cent, col = 3, type = "b")
lines(6:10, corr.rank.mscent, col = 4, type = "b")
legend("bottomright", legend = c("Multiseq", "Dcuts", "Centipede", "msCentipede"), lty = c(1, 1, 1, 1), col = 1:4)
dev.off()

pdf(paste(factor.name, "corr_spearman_all.pdf", sep = "_"))
plot(6:10, corr.rank.dcut.all, col = 2, type = "b", ylim = c(0, 1), xlab = "window size (log2)", ylab = "correlation", main = "all motifs, rank correlation")
lines(6:10, corr.rank.cent.all, col = 3, type = "b")
lines(6:10, corr.rank.mscent.all, col = 4, type = "b")
legend("bottomright", legend = c("Dcuts", "Centipede", "msCentipede"), lty = c(1, 1, 1), col = 2:4)
dev.off()

save(auc.ms,
     auc.dcut,
     auc.cent,
     auc.mscent,
     auc.dcut.all,
     auc.cent.all,
     auc.mscent.all,
     corr.ms,
     corr.dcut,
     corr.rank.ms,
     corr.rank.dcut,
     corr.rank.cent,
     corr.rank.mscent,
     corr.rank.dcut.all,
     corr.rank.cent.all,
     corr.rank.mscent.all,
     file = paste(factor.name, "res.Robj", sep = "_"))