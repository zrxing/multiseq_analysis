library(multiseq)
library(CENTIPEDE)

sep.strand = FALSE

nbase = 256

base.ind.for.cur = base.ind.for[[2]]
base.ind.rev.cur = base.ind.rev[[2]]


nr = 10
eff.for = matrix(0,nr,nbase)
eff.rev = matrix(0,nr,nbase)
base.for = matrix(0,nr,nbase)
base.rev = matrix(0,nr,nbase)

for(i in 1:nr){
  ind = sample(1:dim(dprof.train)[1],3000)
  est.for = multiseq(as.matrix(dprof.train[ind,base.ind.for.cur]),log(ccount.train[ind]+1),lm.approx = TRUE)
  est.rev = multiseq(as.matrix(dprof.train[ind,base.ind.rev.cur]),log(ccount.train[ind]+1),lm.approx = TRUE)
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

prior.breaks.ini = function(x0,n){
  x = c(-x0,x0)
  for(i in 3:n){
    x[i] = 2 * log(i - 1) - x[i - 1]
  }
  return(x)
}
prior.breaks = c(prior.breaks.ini(0.445,11),seq(2.5,8.5,0.1))
ccount.prior.info = hist(log(ccount.train+1),breaks = prior.breaks, plot = FALSE)
ccount.prior.val = ccount.prior.info$mids
ccount.prior.prob = normalize(smooth.spline(ccount.prior.info$mids,ccount.prior.info$counts)$y)

##computing the likelihood
lambda.for = exp(rep(1,length(ccount.prior.val))%o%base.mean.for + ccount.prior.val%o%eff.mean.for)
lambda.rev = exp(rep(1,length(ccount.prior.val))%o%base.mean.rev + ccount.prior.val%o%eff.mean.rev)
lik.for = matrix(0,length(ccount.test),length(ccount.prior.val))
lik.rev = matrix(0,length(ccount.test),length(ccount.prior.val))
for(i in 1:length(ccount.test)){
#   loglik.ini.for = rowSums(t(apply(lambda.for,1,dpois,x = as.numeric(dprof.test[i,base.ind.for.cur]),log = TRUE)))
#   loglik.ini.for = loglik.ini.for - max(loglik.ini.for)
#   loglik.ini.rev = rowSums(t(apply(lambda.rev,1,dpois,x = as.numeric(dprof.test[i,base.ind.rev.cur]),log = TRUE)))
#   loglik.ini.rev = loglik.ini.rev - max(loglik.ini.rev)
  loglik.ini.for = dpois(sum(as.numeric(dprof.test[i, base.ind.for.cur])), rowSums(lambda.for), log = TRUE) + apply(lambda.for, 1, dmultinom, x = as.numeric(dprof.test[i, base.ind.for.cur]), size = NULL, log = TRUE)
  loglik.ini.for = loglik.ini.for - max(loglik.ini.for)
  loglik.ini.rev = dpois(sum(as.numeric(dprof.test[i, base.ind.rev.cur])), rowSums(lambda.rev), log = TRUE) + apply(lambda.rev, 1, dmultinom, x = as.numeric(dprof.test[i, base.ind.rev.cur]), size = NULL, log = TRUE)
  loglik.ini.rev = loglik.ini.rev - max(loglik.ini.rev)
  lik.for[i,] = exp(loglik.ini.for)
  lik.rev[i,] = exp(loglik.ini.rev)
}

##computing the posterior
if(sep.strand){
  ccount.post.val = ccount.prior.val
  ccount.post.prob.for = lik.for * (rep(1,length(ccount.test))%o%ccount.prior.prob)
  #ccount.post.prob.rev = lik.rev * (rep(1,length(ccount.test))%o%ccount.prior.prob)
  ccount.post.prob.for = t(apply(ccount.post.prob.for,1,normalize))
  #ccount.post.prob.rev = t(apply(ccount.post.prob.rev,1,normalize))
  
  ccount.post.mean.for = 0  
  ccount.post.mode.for = 0
  ccount.post.logmean.for = 0
  #ccount.post.mean.rev = 0
  #ccount.post.mode.rev = 0
  #ccount.post.logmean.rev = 0
  for(i in 1:length(ccount.test)){
    ccount.post.mode.for[i] = exp(ccount.post.val[which(ccount.post.prob.for[i,]==max(ccount.post.prob.for[i,]))])
    ccount.post.mean.for[i] = sum(exp(ccount.post.val)*ccount.post.prob.for[i,])
    ccount.post.logmean.for[i] = sum(ccount.post.val*ccount.post.prob.for[i,])
    #  ccount.post.mode.rev[i] = exp(ccount.post.val[which(ccount.post.prob.rev[i,]==max(ccount.post.prob.rev[i,]))])
    #  ccount.post.mean.rev[i] = sum(exp(ccount.post.val)*ccount.post.prob.rev[i,])
    #  ccount.post.logmean.rev[i] = sum(ccount.post.val*ccount.post.prob.rev[i,])
    
  }
  plot(ccount.post.logmean.for,log(ccount.test+1), xlim = c(0, 9), ylim = c(0, 9))
  abline(0, 1, col = 2)
  
  cor(ccount.post.logmean.for,log(ccount.test+1))
}else{
  ccount.post.val = ccount.prior.val
  ccount.post.prob = lik.for * lik.rev * (rep(1,length(ccount.test))%o%ccount.prior.prob)
  ccount.post.prob = t(apply(ccount.post.prob,1,normalize))
  
  ccount.post.mean = 0  
  ccount.post.mode = 0
  ccount.post.logmean = 0
  ccount.post.prob.bound = 0
  ind.bound = ccount.prior.val > quantile(log(ccount+1),0.5)
  ind.unbound = ccount.prior.val <= quantile(log(ccount+1),0.5)
  
  
  for(i in 1:length(ccount.test)){
    ccount.post.mode[i] = exp(ccount.post.val[which(ccount.post.prob[i,]==max(ccount.post.prob[i,]))])
    ccount.post.mean[i] = sum(exp(ccount.post.val)*ccount.post.prob[i,])
    ccount.post.logmean[i] = sum(ccount.post.val*ccount.post.prob[i,])
    ccount.post.prob.bound[i] = exp(log(sum(ccount.post.prob[i, ind.bound]))-log(sum(ccount.post.prob[i, ind.bound | ind.unbound])))
  } 
  plot(ccount.post.logmean,log(ccount.test+1), xlim = c(0, 9), ylim = c(0, 9))
  abline(0, 1, col = 2)
  
  cor(ccount.post.logmean,log(ccount.test+1))
}

plot(ccount.post.prob.bound[bound.ind], log(ccount.test[bound.ind]+1), xlim = c(0, 1), ylim = c(0, 9))
points(ccount.post.prob.bound[unbound.ind], log(ccount.test[unbound.ind]+1), col = 2)


plot(centFit$PostPr[bound.ind], log(ccount.test[bound.ind]+1), xlim = c(0, 1), ylim = c(0, 9))
points(centFit$PostPr[unbound.ind], log(ccount.test[unbound.ind]+1), col = 2)

roc.res.ms = rocplot(ccount.post.logmean[unbound.ind], ccount.post.logmean[bound.ind])
plot(roc.res.ms$fpr,roc.res.ms$tpr, type = 'l', xlim = c(0, 1))
centFit <- fitCentipede(Xlist = list(DNase=as.matrix(dprof[, c(413:612, 1437:1636)])), Y=matrix(rep(1,dim(dprof)[1], nc = 1)))
roc.res.cent = rocplot(centFit$PostPr[!train.ind][unbound.ind], centFit$PostPr[!train.ind][bound.ind])
lines(roc.res.cent$fpr,roc.res.cent$tpr, col = 2)
dcut = rowSums(dprof.test[, c(413:612, 1437:1636)])
roc.res.dcut = rocplot(dcut[unbound.ind], dcut[bound.ind])
lines(roc.res.dcut$fpr,roc.res.dcut$tpr, col = 4)



centFit <- fitCentipede(Xlist = list(DNase=as.matrix(dprof.test[, c(413:612, 1437:1636)])),DampLambda = 0, DampNegBin = 0)

misclass = which(unbound.ind == TRUE & ccount.post.prob.bound > 0.95)

#5880

centFit$PriorLogRatio[misclass]
centFit$MultiNomLogRatio[misclass]
centFit$NegBinLogRatio[misclass]


as.numeric(apply(as.matrix(dprof.test[misclass, base.ind.for.cur]), 1, dmultinom, size = NULL, prob = lambda.for[37, ], log = TRUE) - apply(as.matrix(dprof.test[misclass, base.ind.for.cur]), 1, dmultinom, size = NULL, prob = lambda.for[5, ], log = TRUE))
as.numeric(dpois(rowSums(as.matrix(dprof.test[misclass, base.ind.for.cur])), lambda = sum(lambda.for[37, ]), log = TRUE) - dpois(rowSums(as.matrix(dprof.test[misclass, base.ind.for.cur])), lambda = sum(lambda.for[5, ]), log = TRUE))

#index 37 corresponds to median of top 50% of chipseq counts, index 5 corresponds to bottom 50%
plot(centFit$MultiNomLogRatio[misclass],as.numeric(apply(as.matrix(dprof.test[misclass, base.ind.for.cur]), 1, dmultinom, size = NULL, prob = lambda.for[37, ], log = TRUE) - apply(as.matrix(dprof.test[misclass, base.ind.for.cur]), 1, dmultinom, size = NULL, prob = lambda.for[5, ], log = TRUE)))
abline(0, 1, col = 2)
plot(centFit$NegBinLogRatio[misclass], as.numeric(dpois(rowSums(as.matrix(dprof.test[misclass, base.ind.for.cur])), lambda = sum(lambda.for[37, ]), log = TRUE) - dpois(rowSums(as.matrix(dprof.test[misclass, base.ind.for.cur])), lambda = sum(lambda.for[5, ]), log = TRUE)))
abline(0, 1, col = 2)

centFit.all = fitCentipede(Xlist = list(DNase=as.matrix(dprof[, c(413:612, 1437:1636)])),DampLambda = 0, DampNegBin = 0)
overall.ratio.cent = centFit.all$NegBinLogRatio
overall.ratio.ms = as.numeric(dpois(rowSums(as.matrix(dprof[, base.ind.for.cur])), lambda = sum(lambda.for[37, ]), log = TRUE) - dpois(rowSums(as.matrix(dprof[, base.ind.for.cur])), lambda = sum(lambda.for[5, ]), log = TRUE))

plot(overall.ratio.cent, overall.ratio.ms)
abline(0, 1, col = 2)
lmfit = lm(overall.ratio.ms[overall.ratio.cent > 0] ~ overall.ratio.cent[overall.ratio.cent > 0])

#Look at multinomial portion across all sites
multi.ratio.cent = centFit.all$MultiNomLogRatio
multi.ratio.ms = as.numeric(apply(as.matrix(dprof[, base.ind.for.cur]), 1, dmultinom, size = NULL, prob = lambda.for[37, ], log = TRUE) - apply(as.matrix(dprof[, base.ind.for.cur]), 1, dmultinom, size = NULL, prob = lambda.for[5, ], log = TRUE))
plot(multi.ratio.cent, multi.ratio.ms)
abline(0, 1, col = 2)

#Look at the ratio of the negbin vs the multinom components
plot(overall.ratio.cent, multi.ratio.cent)
abline(0, 1, col = 2)
abline(0, -1, col = 2)

plot(overall.ratio.ms, multi.ratio.ms, xlim = c(0, 1000), ylim = c(-200, 100))
abline(0, 1, col = 2)
abline(0, -1, col = 2)
