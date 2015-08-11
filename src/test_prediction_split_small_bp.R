library(multiseq)
load("results/test_prediction_bp_for_n_rev.RData")

sep.strand = FALSE

nbase = 128

base.ind.for.cur = base.ind.for[[1]]
base.ind.rev.cur = base.ind.rev[[1]]


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
  loglik.ini.for = rowSums(t(apply(lambda.for,1,dpois,x = as.numeric(dprof.test[i,base.ind.for.cur]),log = TRUE)))
  loglik.ini.for = loglik.ini.for - max(loglik.ini.for)
  loglik.ini.rev = rowSums(t(apply(lambda.rev,1,dpois,x = as.numeric(dprof.test[i,base.ind.rev.cur]),log = TRUE)))
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
  
  for(i in 1:length(ccount.test)){
    ccount.post.mode[i] = exp(ccount.post.val[which(ccount.post.prob[i,]==max(ccount.post.prob[i,]))])
    ccount.post.mean[i] = sum(exp(ccount.post.val)*ccount.post.prob[i,])
    ccount.post.logmean[i] = sum(ccount.post.val*ccount.post.prob[i,])
    
  } 
  plot(ccount.post.logmean,log(ccount.test+1), xlim = c(0, 9), ylim = c(0, 9))
  abline(0, 1, col = 2)
  
  cor(ccount.post.logmean,log(ccount.test+1))
}




#split into 10 equally spaced bins based on quantiles, and run smash on each bin separately
##prior
lccount.train = log(ccount.train + 1)
bin.cuts = quantile(lccount.train, probs = 1:10/10)
ccount.prior.val.quant.sep = as.vector(filter(c(0, as.vector(bin.cuts)), c(0.5, 0.5), method = "convolution", sides = 1)[-1])
ccount.prior.prob.quant.sep = c(1/10, rep(1/10, 9))
nq = length(bin.cuts)

lambda.for.quant.sep = matrix(0, nr = nq, nc = nbase)
lambda.rev.quant.sep = matrix(0, nr = nq, nc = nbase)

dprof.train.temp = dprof.train[lccount.train <= bin.cuts[1], ]
lambda.for.quant.sep[1, ] = exp(multiseq(as.matrix(dprof.train.temp[, base.ind.for.cur]), lm.approx = FALSE)$baseline.mean)
lambda.rev.quant.sep[1, ] = exp(multiseq(as.matrix(dprof.train.temp[, base.ind.rev.cur]), lm.approx = FALSE)$baseline.mean)
for(i in 2:nq){
  dprof.train.temp = dprof.train[(lccount.train > bin.cuts[i-1]) & (lccount.train <= bin.cuts[i]), ]
  lambda.for.quant.sep[i, ] = exp(multiseq(as.matrix(dprof.train.temp[, base.ind.for.cur]), lm.approx = FALSE)$baseline.mean)
  lambda.rev.quant.sep[i, ] = exp(multiseq(as.matrix(dprof.train.temp[, base.ind.rev.cur]), lm.approx = FALSE)$baseline.mean)
}

##computing the likelihood
lik.for.quant.sep = matrix(0, length(ccount.test), nq)
lik.rev.quant.sep = matrix(0, length(ccount.test), nq)
for(i in 1:length(ccount.test)){
  loglik.ini.for = rowSums(t(apply(lambda.for.quant.sep, 1, dpois, x = as.numeric(dprof.test[i, base.ind.for.cur]), log = TRUE)))
  loglik.ini.for = loglik.ini.for - max(loglik.ini.for)
  loglik.ini.rev = rowSums(t(apply(lambda.rev.quant.sep,1,dpois,x = as.numeric(dprof.test[i,base.ind.rev.cur]),log = TRUE)))
  loglik.ini.rev = loglik.ini.rev - max(loglik.ini.rev)
  lik.for.quant.sep[i,] = exp(loglik.ini.for)
  lik.rev.quant.sep[i,] = exp(loglik.ini.rev)
}

##computing the posterior
if(sep.strand){
  ccount.post.val.quant.sep = ccount.prior.val.quant.sep
  ccount.post.prob.for.quant.sep = lik.for.quant.sep * (rep(1,length(ccount.test))%o%ccount.prior.prob.quant.sep)
  #ccount.post.prob.rev.quant.sep = lik.rev.quant.sep * (rep(1,length(ccount.test))%o%ccount.prior.prob.quant.sep)
  ccount.post.prob.for.quant.sep = t(apply(ccount.post.prob.for.quant.sep,1,normalize))
  #ccount.post.prob.rev.quant.sep = t(apply(ccount.post.prob.rev.quant.sep,1,normalize))
  
  ccount.post.mean.for.quant.sep = 0  
  ccount.post.mode.for.quant.sep = 0
  ccount.post.logmean.for.quant.sep = 0
  #ccount.post.mean.rev.quant.sep = 0
  #ccount.post.mode.rev.quant.sep = 0
  #ccount.post.logmean.rev.quant.sep = 0
  for(i in 1:length(ccount.test)){
    ccount.post.mode.for.quant.sep[i] = exp(ccount.post.val.quant.sep[which(ccount.post.prob.for.quant.sep[i,]==max(ccount.post.prob.for.quant.sep[i,]))])
    ccount.post.mean.for.quant.sep[i] = sum(exp(ccount.post.val.quant.sep)*ccount.post.prob.for.quant.sep[i,])
    ccount.post.logmean.for.quant.sep[i] = sum(ccount.post.val.quant.sep*ccount.post.prob.for.quant.sep[i,])
    #  ccount.post.mode.rev.quant.sep[i] = exp(ccount.post.val.quant.sep[which(ccount.post.prob.rev.quant.sep[i,]==max(ccount.post.prob.rev.quant.sep[i,]))])
    #  ccount.post.mean.rev.quant.sep[i] = sum(exp(ccount.post.val.quant.sep)*ccount.post.prob.rev.quant.sep[i,])
    #  ccount.post.logmean.rev.quant.sep[i] = sum(ccount.post.val.quant.sep*ccount.post.prob.rev.quant.sep[i,])
    
  }
  plot(ccount.post.logmean.for.quant.sep,log(ccount.test+1), xlim = c(0, 9), ylim = c(0, 9))
  abline(0, 1, col = 2)
  
  cor(ccount.post.logmean.for.quant.sep,log(ccount.test+1))
}else{
  ccount.post.val.quant.sep = ccount.prior.val.quant.sep
  ccount.post.prob.quant.sep = lik.for.quant.sep * lik.rev.quant.sep * (rep(1,length(ccount.test))%o%ccount.prior.prob.quant.sep)
  ccount.post.prob.quant.sep = t(apply(ccount.post.prob.quant.sep,1,normalize))
  
  ccount.post.mean.quant.sep = 0  
  ccount.post.mode.quant.sep = 0
  ccount.post.logmean.quant.sep = 0
  
  for(i in 1:length(ccount.test)){
    ccount.post.mode.quant.sep[i] = exp(ccount.post.val.quant.sep[which(ccount.post.prob.quant.sep[i,]==max(ccount.post.prob.quant.sep[i,]))])
    ccount.post.mean.quant.sep[i] = sum(exp(ccount.post.val.quant.sep)*ccount.post.prob.quant.sep[i,])
    ccount.post.logmean.quant.sep[i] = sum(ccount.post.val.quant.sep*ccount.post.prob.quant.sep[i,])
    
  } 
  plot(ccount.post.logmean.quant.sep,log(ccount.test+1), xlim = c(0, 9), ylim = c(0, 9))
  abline(0, 1, col = 2)
  
  cor(ccount.post.logmean.quant.sep,log(ccount.test+1))
}



#split into 9 bins based on log chipseq count values, and run smash on each bin separately
lccount.train = log(ccount.train + 1)
bin.cuts = 1:9
ccount.prior.val.quant.sep.bin = as.vector(filter(c(0, as.vector(bin.cuts)), c(0.5, 0.5), method = "convolution", sides = 1)[-1])
ccount.prior.prob.quant.sep.bin = hist(lccount.train, breaks = 0:9 ,plot = FALSE)$density
nq = length(bin.cuts)

lambda.for.quant.sep.bin = matrix(0, nr = nq, nc = nbase)
lambda.rev.quant.sep.bin = matrix(0, nr = nq, nc = nbase)
dprof.train.temp = dprof.train[lccount.train <= bin.cuts[1], ]
lambda.for.quant.sep.bin[1, ] = exp(multiseq(as.matrix(dprof.train.temp[, base.ind.for.cur]), lm.approx = TRUE)$baseline.mean)
lambda.rev.quant.sep.bin[1, ] = exp(multiseq(as.matrix(dprof.train.temp[, base.ind.rev.cur]), lm.approx = TRUE)$baseline.mean)
for(i in 2:nq){
  dprof.train.temp = dprof.train[(lccount.train > bin.cuts[i-1]) & (lccount.train <= bin.cuts[i]), ]
  lambda.for.quant.sep.bin[i, ] = exp(multiseq(as.matrix(dprof.train.temp[, base.ind.for.cur]), lm.approx = TRUE)$baseline.mean)
  lambda.rev.quant.sep.bin[i, ] = exp(multiseq(as.matrix(dprof.train.temp[, base.ind.rev.cur]), lm.approx = TRUE)$baseline.mean)
}

##computing the likelihood
lik.for.quant.sep.bin = matrix(0, length(ccount.test), nq)
lik.rev.quant.sep.bin = matrix(0, length(ccount.test), nq)
for(i in 1:length(ccount.test)){
  loglik.ini.for = rowSums(t(apply(lambda.for.quant.sep.bin, 1, dpois, x = as.numeric(dprof.test[i, base.ind.for.cur]), log = TRUE)))
  loglik.ini.for = loglik.ini.for - max(loglik.ini.for)
  loglik.ini.rev = rowSums(t(apply(lambda.rev.quant.sep.bin,1,dpois,x = as.numeric(dprof.test[i,base.ind.rev.cur]),log = TRUE)))
  loglik.ini.rev = loglik.ini.rev - max(loglik.ini.rev)
  lik.for.quant.sep.bin[i,] = exp(loglik.ini.for)
  lik.rev.quant.sep.bin[i,] = exp(loglik.ini.rev)
}

##computing the posterior
if(sep.strand){
  ccount.post.val.quant.sep.bin = ccount.prior.val.quant.sep.bin
  ccount.post.prob.for.quant.sep.bin = lik.for.quant.sep.bin * (rep(1,length(ccount.test))%o%ccount.prior.prob.quant.sep.bin)
  #ccount.post.prob.rev.quant.sep.bin = lik.rev.quant.sep.bin * (rep(1,length(ccount.test))%o%ccount.prior.prob.quant.sep.bin)
  ccount.post.prob.for.quant.sep.bin = t(apply(ccount.post.prob.for.quant.sep.bin,1,normalize))
  #ccount.post.prob.rev.quant.sep.bin = t(apply(ccount.post.prob.rev.quant.sep.bin,1,normalize))
  
  ccount.post.mean.for.quant.sep.bin = 0  
  ccount.post.mode.for.quant.sep.bin = 0
  ccount.post.logmean.for.quant.sep.bin = 0
  #ccount.post.mean.rev.quant.sep.bin = 0
  #ccount.post.mode.rev.quant.sep.bin = 0
  #ccount.post.logmean.rev.quant.sep.bin = 0
  for(i in 1:length(ccount.test)){
    ccount.post.mode.for.quant.sep.bin[i] = exp(ccount.post.val.quant.sep.bin[which(ccount.post.prob.for.quant.sep.bin[i,]==max(ccount.post.prob.for.quant.sep.bin[i,]))])
    ccount.post.mean.for.quant.sep.bin[i] = sum(exp(ccount.post.val.quant.sep.bin)*ccount.post.prob.for.quant.sep.bin[i,])
    ccount.post.logmean.for.quant.sep.bin[i] = sum(ccount.post.val.quant.sep.bin*ccount.post.prob.for.quant.sep.bin[i,])
    #  ccount.post.mode.rev.quant.sep.bin[i] = exp(ccount.post.val.quant.sep.bin[which(ccount.post.prob.rev.quant.sep.bin[i,]==max(ccount.post.prob.rev.quant.sep.bin[i,]))])
    #  ccount.post.mean.rev.quant.sep.bin[i] = sum(exp(ccount.post.val.quant.sep.bin)*ccount.post.prob.rev.quant.sep.bin[i,])
    #  ccount.post.logmean.rev.quant.sep.bin[i] = sum(ccount.post.val.quant.sep.bin*ccount.post.prob.rev.quant.sep.bin[i,])
    
  }
  plot(ccount.post.logmean.for.quant.sep.bin,log(ccount.test+1), xlim = c(0, 9), ylim = c(0, 9))
  abline(0, 1, col = 2)
  
  cor(ccount.post.logmean.for.quant.sep.bin,log(ccount.test+1))
}else{
  ccount.post.val.quant.sep.bin = ccount.prior.val.quant.sep.bin
  ccount.post.prob.quant.sep.bin = lik.for.quant.sep.bin * lik.rev.quant.sep.bin * (rep(1,length(ccount.test))%o%ccount.prior.prob.quant.sep.bin)
  ccount.post.prob.quant.sep.bin = t(apply(ccount.post.prob.quant.sep.bin,1,normalize))
  
  ccount.post.mean.quant.sep.bin = 0  
  ccount.post.mode.quant.sep.bin = 0
  ccount.post.logmean.quant.sep.bin = 0
  
  for(i in 1:length(ccount.test)){
    ccount.post.mode.quant.sep.bin[i] = exp(ccount.post.val.quant.sep.bin[which(ccount.post.prob.quant.sep.bin[i,]==max(ccount.post.prob.quant.sep.bin[i,]))])
    ccount.post.mean.quant.sep.bin[i] = sum(exp(ccount.post.val.quant.sep.bin)*ccount.post.prob.quant.sep.bin[i,])
    ccount.post.logmean.quant.sep.bin[i] = sum(ccount.post.val.quant.sep.bin*ccount.post.prob.quant.sep.bin[i,])
    
  } 
  plot(ccount.post.logmean.quant.sep.bin,log(ccount.test+1), xlim = c(0, 9), ylim = c(0, 9))
  abline(0, 1, col = 2)
  
  cor(ccount.post.logmean.quant.sep.bin,log(ccount.test+1))
}


#split into 10 equally spaced bins based on quantiles, and run multiseq on adjacent bins
#single baseline from original analysis; add effect from each bin to baseline
##prior
lccount.train = log(ccount.train + 1)
bin.cuts = quantile(lccount.train, probs = 1:10/10)
ccount.prior.val.quant.common.base = as.vector(filter(c(0, as.vector(bin.cuts)), c(0.5, 0.5), method = "convolution", sides = 1)[-1])
ccount.prior.prob.quant.common.base = rep(1/10, 10)
nq = length(bin.cuts)

lambda.for.quant.common.base = matrix(0, nr = nq, nc = nbase)
lambda.for.quant.common.base[1, ] = base.mean.for
lambda.rev.quant.common.base = matrix(0, nr = nq, nc = nbase)
lambda.rev.quant.common.base[1, ] = base.mean.rev
dprof.train.temp1 = dprof.train[lccount.train <= bin.cuts[1], ]
dprof.train.temp2 = dprof.train[(lccount.train > bin.cuts[1]) & (lccount.train <= bin.cuts[2]), ]
dprof.train.temp = rbind(dprof.train.temp1, dprof.train.temp2)
g.temp = c(rep(0, dim(dprof.train.temp1)[1]), rep(1, dim(dprof.train.temp2)[1]))
lambda.for.quant.common.base[2, ] = lambda.for.quant.common.base[1, ] + multiseq(as.matrix(dprof.train.temp[, base.ind.for.cur]), g.temp, lm.approx = TRUE)$effect.mean
lambda.rev.quant.common.base[2, ] = lambda.rev.quant.common.base[1, ] + multiseq(as.matrix(dprof.train.temp[, base.ind.rev.cur]), g.temp, lm.approx = TRUE)$effect.mean
for(i in 3:nq){
  dprof.train.temp1 = dprof.train[(lccount.train > bin.cuts[i-2]) & (lccount.train <= bin.cuts[i-1]), ]
  dprof.train.temp2 = dprof.train[(lccount.train > bin.cuts[i-1]) & (lccount.train <= bin.cuts[i]), ]
  dprof.train.temp = rbind(dprof.train.temp1, dprof.train.temp2)
  g.temp = c(rep(0, dim(dprof.train.temp1)[1]), rep(1, dim(dprof.train.temp2)[1]))
  lambda.for.quant.common.base[i, ] = lambda.for.quant.common.base[i-1, ] + multiseq(as.matrix(dprof.train.temp[, base.ind.for.cur]), g.temp, lm.approx = TRUE)$effect.mean
  lambda.rev.quant.common.base[i, ] = lambda.rev.quant.common.base[i-1, ] + multiseq(as.matrix(dprof.train.temp[, base.ind.rev.cur]), g.temp, lm.approx = TRUE)$effect.mean
}
lambda.for.quant.common.base = exp(lambda.for.quant.common.base)
lambda.rev.quant.common.base = exp(lambda.rev.quant.common.base)

##computing the likelihood
lik.for.quant.common.base = matrix(0, length(ccount.test), nq)
lik.rev.quant.common.base = matrix(0, length(ccount.test), nq)
for(i in 1:length(ccount.test)){
  loglik.ini.for = rowSums(t(apply(lambda.for.quant.common.base, 1, dpois, x = as.numeric(dprof.test[i, base.ind.for.cur]), log = TRUE)))
  loglik.ini.for = loglik.ini.for - max(loglik.ini.for)
  loglik.ini.rev = rowSums(t(apply(lambda.rev.quant.common.base,1,dpois,x = as.numeric(dprof.test[i,base.ind.rev.cur]),log = TRUE)))
  loglik.ini.rev = loglik.ini.rev - max(loglik.ini.rev)
  lik.for.quant.common.base[i,] = exp(loglik.ini.for)
  lik.rev.quant.common.base[i,] = exp(loglik.ini.rev)
}

##computing the posterior
if(sep.strand){
  ccount.post.val.quant.common.base = ccount.prior.val.quant.common.base
  ccount.post.prob.for.quant.common.base = lik.for.quant.common.base * (rep(1,length(ccount.test))%o%ccount.prior.prob.quant.common.base)
  #ccount.post.prob.rev.quant.common.base = lik.rev.quant.common.base * (rep(1,length(ccount.test))%o%ccount.prior.prob.quant.common.base)
  ccount.post.prob.for.quant.common.base = t(apply(ccount.post.prob.for.quant.common.base,1,normalize))
  #ccount.post.prob.rev.quant.common.base = t(apply(ccount.post.prob.rev.quant.common.base,1,normalize))
  
  ccount.post.mean.for.quant.common.base = 0  
  ccount.post.mode.for.quant.common.base = 0
  ccount.post.logmean.for.quant.common.base = 0
  #ccount.post.mean.rev.quant.common.base = 0
  #ccount.post.mode.rev.quant.common.base = 0
  #ccount.post.logmean.rev.quant.common.base = 0
  for(i in 1:length(ccount.test)){
    ccount.post.mode.for.quant.common.base[i] = exp(ccount.post.val.quant.common.base[which(ccount.post.prob.for.quant.common.base[i,]==max(ccount.post.prob.for.quant.common.base[i,]))])
    ccount.post.mean.for.quant.common.base[i] = sum(exp(ccount.post.val.quant.common.base)*ccount.post.prob.for.quant.common.base[i,])
    ccount.post.logmean.for.quant.common.base[i] = sum(ccount.post.val.quant.common.base*ccount.post.prob.for.quant.common.base[i,])
    #  ccount.post.mode.rev.quant.common.base[i] = exp(ccount.post.val.quant.common.base[which(ccount.post.prob.rev.quant.common.base[i,]==max(ccount.post.prob.rev.quant.common.base[i,]))])
    #  ccount.post.mean.rev.quant.common.base[i] = sum(exp(ccount.post.val.quant.common.base)*ccount.post.prob.rev.quant.common.base[i,])
    #  ccount.post.logmean.rev.quant.common.base[i] = sum(ccount.post.val.quant.common.base*ccount.post.prob.rev.quant.common.base[i,])
    
  }
  plot(ccount.post.logmean.for.quant.common.base,log(ccount.test+1), xlim = c(0, 9), ylim = c(0, 9))
  abline(0, 1, col = 2)
  
  cor(ccount.post.logmean.for.quant.common.base,log(ccount.test+1))
}else{
  ccount.post.val.quant.common.base = ccount.prior.val.quant.common.base
  ccount.post.prob.quant.common.base = lik.for.quant.common.base * lik.rev.quant.common.base * (rep(1,length(ccount.test))%o%ccount.prior.prob.quant.common.base)
  ccount.post.prob.quant.common.base = t(apply(ccount.post.prob.quant.common.base,1,normalize))
  
  ccount.post.mean.quant.common.base = 0  
  ccount.post.mode.quant.common.base = 0
  ccount.post.logmean.quant.common.base = 0
  
  for(i in 1:length(ccount.test)){
    ccount.post.mode.quant.common.base[i] = exp(ccount.post.val.quant.common.base[which(ccount.post.prob.quant.common.base[i,]==max(ccount.post.prob.quant.common.base[i,]))])
    ccount.post.mean.quant.common.base[i] = sum(exp(ccount.post.val.quant.common.base)*ccount.post.prob.quant.common.base[i,])
    ccount.post.logmean.quant.common.base[i] = sum(ccount.post.val.quant.common.base*ccount.post.prob.quant.common.base[i,])
    
  } 
  plot(ccount.post.logmean.quant.common.base,log(ccount.test+1), xlim = c(0, 9), ylim = c(0, 9))
  abline(0, 1, col = 2)
  
  cor(ccount.post.logmean.quant.common.base,log(ccount.test+1))
}






#single baseline from first two bins; add effect from each bin to baseline
##prior
lccount.train = log(ccount.train + 1)
bin.cuts = quantile(lccount.train, probs = 1:10/10)
ccount.prior.val.quant.single.base = as.vector(filter(c(0, as.vector(bin.cuts)), c(0.5, 0.5), method = "convolution", sides = 1)[-1])
ccount.prior.prob.quant.single.base = rep(1/10, 10)
nq = length(bin.cuts)

lambda.for.quant.single.base = matrix(0, nr = nq, nc = nbase)
lambda.rev.quant.single.base = matrix(0, nr = nq, nc = nbase)
dprof.train.temp1 = dprof.train[lccount.train <= bin.cuts[1], ]
dprof.train.temp2 = dprof.train[(lccount.train > bin.cuts[1]) & (lccount.train <= bin.cuts[2]), ]
dprof.train.temp = rbind(dprof.train.temp1, dprof.train.temp2)
g.temp = c(rep(0, dim(dprof.train.temp1)[1]), rep(1, dim(dprof.train.temp2)[1]))
multiseq.out.temp = multiseq(as.matrix(dprof.train.temp[, base.ind.for.cur]), g.temp, lm.approx = TRUE)
lambda.for.quant.single.base[1, ] = multiseq.out.temp$baseline.mean
lambda.for.quant.single.base[2, ] = lambda.for.quant.single.base[1, ] + multiseq.out.temp$effect.mean
multiseq.out.temp = multiseq(as.matrix(dprof.train.temp[, base.ind.rev.cur]), g.temp, lm.approx = TRUE)
lambda.rev.quant.single.base[1, ] = multiseq.out.temp$baseline.mean
lambda.rev.quant.single.base[2, ] = lambda.rev.quant.single.base[1, ] + multiseq.out.temp$effect.mean
for(i in 3:nq){
  dprof.train.temp1 = dprof.train[(lccount.train > bin.cuts[i-2]) & (lccount.train <= bin.cuts[i-1]), ]
  dprof.train.temp2 = dprof.train[(lccount.train > bin.cuts[i-1]) & (lccount.train <= bin.cuts[i]), ]
  dprof.train.temp = rbind(dprof.train.temp1, dprof.train.temp2)
  g.temp = c(rep(0, dim(dprof.train.temp1)[1]), rep(1, dim(dprof.train.temp2)[1]))
  lambda.for.quant.single.base[i, ] = lambda.for.quant.single.base[i-1, ] + multiseq(as.matrix(dprof.train.temp[, base.ind.for.cur]), g.temp, lm.approx = TRUE)$effect.mean
  lambda.rev.quant.single.base[i, ] = lambda.rev.quant.single.base[i-1, ] + multiseq(as.matrix(dprof.train.temp[, base.ind.rev.cur]), g.temp, lm.approx = TRUE)$effect.mean
}
lambda.for.quant.single.base = exp(lambda.for.quant.single.base)
lambda.rev.quant.single.base = exp(lambda.rev.quant.single.base)

##computing the likelihood
lik.for.quant.single.base = matrix(0, length(ccount.test), nq)
lik.rev.quant.single.base = matrix(0, length(ccount.test), nq)
for(i in 1:length(ccount.test)){
  loglik.ini.for = rowSums(t(apply(lambda.for.quant.single.base, 1, dpois, x = as.numeric(dprof.test[i, base.ind.for.cur]), log = TRUE)))
  loglik.ini.for = loglik.ini.for - max(loglik.ini.for)
  loglik.ini.rev = rowSums(t(apply(lambda.rev.quant.single.base,1,dpois,x = as.numeric(dprof.test[i,base.ind.rev.cur]),log = TRUE)))
  loglik.ini.rev = loglik.ini.rev - max(loglik.ini.rev)
  lik.for.quant.single.base[i,] = exp(loglik.ini.for)
  lik.rev.quant.single.base[i,] = exp(loglik.ini.rev)
}

##computing the posterior
if(sep.strand){
  ccount.post.val.quant.single.base = ccount.prior.val.quant.single.base
  ccount.post.prob.for.quant.single.base = lik.for.quant.single.base * (rep(1,length(ccount.test))%o%ccount.prior.prob.quant.single.base)
  #ccount.post.prob.rev.quant.single.base = lik.rev.quant.single.base * (rep(1,length(ccount.test))%o%ccount.prior.prob.quant.single.base)
  ccount.post.prob.for.quant.single.base = t(apply(ccount.post.prob.for.quant.single.base,1,normalize))
  #ccount.post.prob.rev.quant.single.base = t(apply(ccount.post.prob.rev.quant.single.base,1,normalize))
  
  ccount.post.mean.for.quant.single.base = 0  
  ccount.post.mode.for.quant.single.base = 0
  ccount.post.logmean.for.quant.single.base = 0
  #ccount.post.mean.rev.quant.single.base = 0
  #ccount.post.mode.rev.quant.single.base = 0
  #ccount.post.logmean.rev.quant.single.base = 0
  for(i in 1:length(ccount.test)){
    ccount.post.mode.for.quant.single.base[i] = exp(ccount.post.val.quant.single.base[which(ccount.post.prob.for.quant.single.base[i,]==max(ccount.post.prob.for.quant.single.base[i,]))])
    ccount.post.mean.for.quant.single.base[i] = sum(exp(ccount.post.val.quant.single.base)*ccount.post.prob.for.quant.single.base[i,])
    ccount.post.logmean.for.quant.single.base[i] = sum(ccount.post.val.quant.single.base*ccount.post.prob.for.quant.single.base[i,])
    #  ccount.post.mode.rev.quant.single.base[i] = exp(ccount.post.val.quant.single.base[which(ccount.post.prob.rev.quant.single.base[i,]==max(ccount.post.prob.rev.quant.single.base[i,]))])
    #  ccount.post.mean.rev.quant.single.base[i] = sum(exp(ccount.post.val.quant.single.base)*ccount.post.prob.rev.quant.single.base[i,])
    #  ccount.post.logmean.rev.quant.single.base[i] = sum(ccount.post.val.quant.single.base*ccount.post.prob.rev.quant.single.base[i,])
    
  }
  plot(ccount.post.logmean.for.quant.single.base,log(ccount.test+1), xlim = c(0, 9), ylim = c(0, 9))
  abline(0, 1, col = 2)
  
  cor(ccount.post.logmean.for.quant.single.base,log(ccount.test+1))
}else{
  ccount.post.val.quant.single.base = ccount.prior.val.quant.single.base
  ccount.post.prob.quant.single.base = lik.for.quant.single.base * lik.rev.quant.single.base * (rep(1,length(ccount.test))%o%ccount.prior.prob.quant.single.base)
  ccount.post.prob.quant.single.base = t(apply(ccount.post.prob.quant.single.base,1,normalize))
  
  ccount.post.mean.quant.single.base = 0  
  ccount.post.mode.quant.single.base = 0
  ccount.post.logmean.quant.single.base = 0
  
  for(i in 1:length(ccount.test)){
    ccount.post.mode.quant.single.base[i] = exp(ccount.post.val.quant.single.base[which(ccount.post.prob.quant.single.base[i,]==max(ccount.post.prob.quant.single.base[i,]))])
    ccount.post.mean.quant.single.base[i] = sum(exp(ccount.post.val.quant.single.base)*ccount.post.prob.quant.single.base[i,])
    ccount.post.logmean.quant.single.base[i] = sum(ccount.post.val.quant.single.base*ccount.post.prob.quant.single.base[i,])
    
  } 
  plot(ccount.post.logmean.quant.single.base,log(ccount.test+1), xlim = c(0, 9), ylim = c(0, 9))
  abline(0, 1, col = 2)
  
  cor(ccount.post.logmean.quant.single.base,log(ccount.test+1))
}






#baseline + effect for each bin separately
##prior
lccount.train = log(ccount.train + 1)
bin.cuts = quantile(lccount.train, probs = 2:20/20)
ccount.prior.val.quant.sep.base = as.vector(filter(c(0, as.vector(bin.cuts)), c(0.8, 0.2), method = "convolution", sides = 1)[-1])
ccount.prior.prob.quant.sep.base = c(2/20, rep(1/20, 18))
nq = length(bin.cuts)

lambda.for.quant.sep.base = matrix(0, nr = nq, nc = nbase)
lambda.rev.quant.sep.base = matrix(0, nr = nq, nc = nbase)
dprof.train.temp1 = dprof.train[lccount.train <= bin.cuts[1], ]
dprof.train.temp2 = dprof.train[(lccount.train > bin.cuts[1]) & (lccount.train <= bin.cuts[2]), ]
dprof.train.temp = rbind(dprof.train.temp1, dprof.train.temp2)
g.temp = c(rep(0, dim(dprof.train.temp1)[1]), rep(1, dim(dprof.train.temp2)[1]))
multiseq.out.temp = multiseq(as.matrix(dprof.train.temp[, base.ind.for.cur]), g.temp, lm.approx = TRUE)
lambda.for.quant.sep.base[1, ] = multiseq.out.temp$baseline.mean
lambda.for.quant.sep.base[2, ] = lambda.for.quant.sep.base[1, ] + multiseq.out.temp$effect.mean
multiseq.out.temp = multiseq(as.matrix(dprof.train.temp[, base.ind.rev.cur]), g.temp, lm.approx = TRUE)
lambda.rev.quant.sep.base[1, ] = multiseq.out.temp$baseline.mean
lambda.rev.quant.sep.base[2, ] = lambda.rev.quant.sep.base[1, ] + multiseq.out.temp$effect.mean
for(i in 3:nq){
  dprof.train.temp1 = dprof.train[(lccount.train > bin.cuts[i-2]) & (lccount.train <= bin.cuts[i-1]), ]
  dprof.train.temp2 = dprof.train[(lccount.train > bin.cuts[i-1]) & (lccount.train <= bin.cuts[i]), ]
  dprof.train.temp = rbind(dprof.train.temp1, dprof.train.temp2)
  g.temp = c(rep(0, dim(dprof.train.temp1)[1]), rep(1, dim(dprof.train.temp2)[1]))
  multiseq.out.temp = multiseq(as.matrix(dprof.train.temp[, base.ind.for.cur]), g.temp, lm.approx = TRUE)
  lambda.for.quant.sep.base[i, ] = multiseq.out.temp$baseline.mean + multiseq.out.temp$effect.mean 
  multiseq.out.temp = multiseq(as.matrix(dprof.train.temp[, base.ind.rev.cur]), g.temp, lm.approx = TRUE)
  lambda.rev.quant.sep.base[i, ] = multiseq.out.temp$baseline.mean + multiseq.out.temp$effect.mean 
  # lambda.for.quant.sep.base[i, ] = (multiseq.out.temp$baseline.mean + lambda.for.quant.sep.base[i-1, ])/2 + multiseq.out.temp$effect.mean 
}
lambda.for.quant.sep.base = exp(lambda.for.quant.sep.base)
lambda.rev.quant.sep.base = exp(lambda.rev.quant.sep.base)


##computing the likelihood
lik.for.quant.sep.base = matrix(0, length(ccount.test), nq)
lik.rev.quant.sep.base = matrix(0, length(ccount.test), nq)
for(i in 1:length(ccount.test)){
  loglik.ini.for = rowSums(t(apply(lambda.for.quant.sep.base, 1, dpois, x = as.numeric(dprof.test[i, base.ind.for.cur]), log = TRUE)))
  loglik.ini.for = loglik.ini.for - max(loglik.ini.for)
  loglik.ini.rev = rowSums(t(apply(lambda.rev.quant.sep.base,1,dpois,x = as.numeric(dprof.test[i,base.ind.rev.cur]),log = TRUE)))
  loglik.ini.rev = loglik.ini.rev - max(loglik.ini.rev)
  lik.for.quant.sep.base[i,] = exp(loglik.ini.for)
  lik.rev.quant.sep.base[i,] = exp(loglik.ini.rev)
}

##computing the posterior
if(sep.strand){
  ccount.post.val.quant.sep.base = ccount.prior.val.quant.sep.base
  ccount.post.prob.for.quant.sep.base = lik.for.quant.sep.base * (rep(1,length(ccount.test))%o%ccount.prior.prob.quant.sep.base)
  #ccount.post.prob.rev.quant.sep.base = lik.rev.quant.sep.base * (rep(1,length(ccount.test))%o%ccount.prior.prob.quant.sep.base)
  ccount.post.prob.for.quant.sep.base = t(apply(ccount.post.prob.for.quant.sep.base,1,normalize))
  #ccount.post.prob.rev.quant.sep.base = t(apply(ccount.post.prob.rev.quant.sep.base,1,normalize))
  
  ccount.post.mean.for.quant.sep.base = 0  
  ccount.post.mode.for.quant.sep.base = 0
  ccount.post.logmean.for.quant.sep.base = 0
  #ccount.post.mean.rev.quant.sep.base = 0
  #ccount.post.mode.rev.quant.sep.base = 0
  #ccount.post.logmean.rev.quant.sep.base = 0
  for(i in 1:length(ccount.test)){
    ccount.post.mode.for.quant.sep.base[i] = exp(ccount.post.val.quant.sep.base[which(ccount.post.prob.for.quant.sep.base[i,]==max(ccount.post.prob.for.quant.sep.base[i,]))])
    ccount.post.mean.for.quant.sep.base[i] = sum(exp(ccount.post.val.quant.sep.base)*ccount.post.prob.for.quant.sep.base[i,])
    ccount.post.logmean.for.quant.sep.base[i] = sum(ccount.post.val.quant.sep.base*ccount.post.prob.for.quant.sep.base[i,])
    #  ccount.post.mode.rev.quant.sep.base[i] = exp(ccount.post.val.quant.sep.base[which(ccount.post.prob.rev.quant.sep.base[i,]==max(ccount.post.prob.rev.quant.sep.base[i,]))])
    #  ccount.post.mean.rev.quant.sep.base[i] = sum(exp(ccount.post.val.quant.sep.base)*ccount.post.prob.rev.quant.sep.base[i,])
    #  ccount.post.logmean.rev.quant.sep.base[i] = sum(ccount.post.val.quant.sep.base*ccount.post.prob.rev.quant.sep.base[i,])
    
  }
  plot(ccount.post.logmean.for.quant.sep.base,log(ccount.test+1), xlim = c(0, 9), ylim = c(0, 9))
  abline(0, 1, col = 2)
  
  cor(ccount.post.logmean.for.quant.sep.base,log(ccount.test+1))
}else{
  ccount.post.val.quant.sep.base = ccount.prior.val.quant.sep.base
  ccount.post.prob.quant.sep.base = lik.for.quant.sep.base * lik.rev.quant.sep.base * (rep(1,length(ccount.test))%o%ccount.prior.prob.quant.sep.base)
  ccount.post.prob.quant.sep.base = t(apply(ccount.post.prob.quant.sep.base,1,normalize))
  
  ccount.post.mean.quant.sep.base = 0  
  ccount.post.mode.quant.sep.base = 0
  ccount.post.logmean.quant.sep.base = 0
  
  for(i in 1:length(ccount.test)){
    ccount.post.mode.quant.sep.base[i] = exp(ccount.post.val.quant.sep.base[which(ccount.post.prob.quant.sep.base[i,]==max(ccount.post.prob.quant.sep.base[i,]))])
    ccount.post.mean.quant.sep.base[i] = sum(exp(ccount.post.val.quant.sep.base)*ccount.post.prob.quant.sep.base[i,])
    ccount.post.logmean.quant.sep.base[i] = sum(ccount.post.val.quant.sep.base*ccount.post.prob.quant.sep.base[i,])
    
  } 
  plot(ccount.post.logmean.quant.sep.base,log(ccount.test+1), xlim = c(0, 9), ylim = c(0, 9))
  abline(0, 1, col = 2)
  
  cor(ccount.post.logmean.quant.sep.base,log(ccount.test+1))
}

