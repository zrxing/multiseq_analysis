library(multiseq)

nbase = 2^7

normalize = function(x) x/sum(x)

dprof = read.table("data/Batf_S356_dnase_data.txt.gz")
cinfo = read.table("data/Batf_S356_site_and_chip_data.txt.gz",header = TRUE)
ccount = cinfo[,6]

train.size = round(0.8*dim(dprof)[1])
set.seed(417)
train.ind = sample(1:dim(dprof)[1],train.size,replace = FALSE)
train.ind = (1:dim(dprof)[1])%in%train.ind

dprof.train = dprof[train.ind,]
dprof.test = dprof[!train.ind,]

ccount.train = ccount[train.ind]
ccount.test = ccount[!train.ind]


peak.info = read.table("data/macs_peaks/Batf_q1pc_peaks.bed")
peak.ind = 0
for(i in 1:dim(cinfo)[1]){
  peak.ind[i] = sum((as.character(peak.info[, 1]) == as.character(cinfo[i, 1])) & (peak.info[, 2] < cinfo[i, 2]) & (peak.info[, 3] > cinfo[i, 3])) > 0
}

unbound.ind = !(peak.ind) & cinfo[, 6] < median(cinfo[, 6])
bound.ind = as.logical(peak.ind)

filter.ind = bound.ind == 1 | unbound.ind == 1

unbound.ind = unbound.ind[filter.ind]
bound.ind = bound.ind[filter.ind]

dprof = dprof[filter.ind, ]
cinfo = cinfo[filter.ind, ]
ccount = cinfo[,6]

train.size = round(0.8*dim(dprof)[1])
set.seed(417)
train.ind = sample(1:dim(dprof)[1],train.size,replace = FALSE)
train.ind = (1:dim(dprof)[1])%in%train.ind

dprof.train = dprof[train.ind,]
dprof.test = dprof[!train.ind,]

ccount.train = ccount[train.ind]
ccount.test = ccount[!train.ind]

unbound.ind.train = unbound.ind[train.ind]
unbound.ind.test = unbound.ind[!train.ind]

bound.ind.train = bound.ind[train.ind]
bound.ind.test = bound.ind[!train.ind]

base.ind.for = 4033:4160 - 3584
base.ind.rev = 12225:12352 - 10752

nr = 10
eff.for = matrix(0,nr,nbase)
eff.rev = matrix(0,nr,nbase)
base.for = matrix(0,nr,nbase)
base.rev = matrix(0,nr,nbase)

for(i in 1:nr){
  ind = sample(1:dim(dprof.train)[1],1000)
  est.for = multiseq(as.matrix(dprof.train[ind,base.ind.for]),as.numeric(bound.ind.train[ind]),lm.approx = TRUE)
  est.rev = multiseq(as.matrix(dprof.train[ind,base.ind.rev]),as.numeric(bound.ind.train[ind]),lm.approx = TRUE)
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

ccount.prior.prob = c(1 - mean(bound.ind.train), mean(bound.ind.train))

##computing the likelihood
lambda.for = exp(matrix(c(base.mean.for, base.mean.for + eff.mean.for), nr = 2, byrow = TRUE))
lambda.rev = exp(matrix(c(base.mean.rev, base.mean.rev + eff.mean.rev), nr = 2, byrow = TRUE))
lik.for = matrix(0, nr = length(ccount.test), nc = 2)
lik.rev = matrix(0, nr = length(ccount.test), nc = 2)
for(i in 1:length(ccount.test)){
  loglik.ini.for = rowSums(t(apply(lambda.for,1,dpois,x = as.numeric(dprof.test[i,base.ind.for]),log = TRUE)))
  loglik.ini.for = loglik.ini.for - max(loglik.ini.for)
  loglik.ini.rev = rowSums(t(apply(lambda.rev,1,dpois,x = as.numeric(dprof.test[i,base.ind.rev]),log = TRUE)))
  loglik.ini.rev = loglik.ini.rev - max(loglik.ini.rev)
  lik.for[i,] = exp(loglik.ini.for)
  lik.rev[i,] = exp(loglik.ini.rev)
}

##computing the posterior
# ccount.post.val = ccount.prior.val
# ccount.post.prob.for = lik.for * (rep(1,length(ccount.test))%o%ccount.prior.prob)
# ccount.post.prob.rev = lik.rev * (rep(1,length(ccount.test))%o%ccount.prior.prob)
# ccount.post.prob.for = t(apply(ccount.post.prob.for,1,normalize))
# ccount.post.prob.rev = t(apply(ccount.post.prob.rev,1,normalize))
# 
# ccount.post.mean.for = 0  
# ccount.post.mode.for = 0
# ccount.post.logmean.for = 0
# ccount.post.mean.rev = 0
# ccount.post.mode.rev = 0
# ccount.post.logmean.rev = 0
# for(i in 1:length(ccount.test)){
#   ccount.post.mode.for[i] = exp(ccount.post.val[which(ccount.post.prob.for[i,]==max(ccount.post.prob.for[i,]))])
#   ccount.post.mean.for[i] = sum(exp(ccount.post.val)*ccount.post.prob.for[i,])
#   ccount.post.logmean.for[i] = sum(ccount.post.val*ccount.post.prob.for[i,])
#   ccount.post.mode.rev[i] = exp(ccount.post.val[which(ccount.post.prob.rev[i,]==max(ccount.post.prob.rev[i,]))])
#   ccount.post.mean.rev[i] = sum(exp(ccount.post.val)*ccount.post.prob.rev[i,])
#   ccount.post.logmean.rev[i] = sum(ccount.post.val*ccount.post.prob.rev[i,])
#   
# }


ccount.post.prob = lik.for * lik.rev * (rep(1,length(ccount.test))%o%ccount.prior.prob)
ccount.post.prob = t(apply(ccount.post.prob,1,normalize))



roc.res.ms = roc(bound.ind.test ~ ccount.post.prob[,1])
plot(roc.res.ms)

pscore = cinfo[,5]
pscore.test = pscore[!train.ind]


centFit <- fitCentipede(Xlist = list(DNase=as.matrix(dprof.test[, c(413:612, 1437:1636)])), Y=matrix(rep(1,dim(dprof.test)[1], nc = 1)), sweeps = 300)
centFit.pwm <- fitCentipede(Xlist = list(DNase=as.matrix(dprof.test)[, c(413:612, 1437:1636)]), Y=cbind(rep(1,dim(dprof.test)[1]),pscore.test))
#roc.res.cent = rocplot(centFit$PostPr[unbound.ind], centFit$PostPr[bound.ind])
#roc.res.cent.pwm = rocplot(centFit.pwm$PostPr[unbound.ind], centFit.pwm$PostPr[bound.ind])

#lines(roc.res.cent$fpr, roc.res.cent$tpr, col = 2)
#lines(roc.res.cent.pwm$fpr, roc.res.cent.pwm$tpr, col = 3)

roc.res.cent = roc(bound.ind.test ~ as.vector(centFit$PostPr))
roc.res.cent.pwm = roc(bound.ind.test ~ as.vector(centFit.pwm$PostPr))

lines(roc.res.cent, col = 2)
lines(roc.res.cent.pwm, col = 3)

dcut = rowSums(dprof.test[, c(413:612, 1437:1636)])
#roc.res.dcut = rocplot(dcut[unbound.ind], dcut[bound.ind])

#lines(roc.res.dcut$fpr,roc.res.dcut$tpr, col = 4)

roc.res.dcut = roc(bound.ind.test ~ dcut)

lines(roc.res.dcut, col = 4)


plot(as.vector(centFit$PostPr), ccount.post.prob)



##look at differences between centipede and multiseq
plot(centFit$MultiNomLogRatio,as.numeric(apply(as.matrix(dprof.test[, base.ind.for]), 1, dmultinom, size = NULL, prob = lambda.for[2, ], log = TRUE) - apply(as.matrix(dprof.test[, base.ind.for]), 1, dmultinom, size = NULL, prob = lambda.for[1, ], log = TRUE)))
abline(0, 1, col = 2)
plot(centFit$NegBinLogRatio, as.numeric(dpois(rowSums(as.matrix(dprof.test[, base.ind.for])), lambda = sum(lambda.for[2, ]), log = TRUE) - dpois(rowSums(as.matrix(dprof.test[, base.ind.for])), lambda = sum(lambda.for[1, ]), log = TRUE)))
abline(0, 1, col = 2)
