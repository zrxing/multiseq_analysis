rsdprof=rowSums(dprof.train)
plot(log(ccount.train+1),log(rsdprof+1))
abline(0,1,col=2)






nr = 20
eff.for.high = matrix(0,nr,nbase)
eff.rev.high = matrix(0,nr,nbase)
base.for.high = matrix(0,nr,nbase)
base.rev.high = matrix(0,nr,nbase)

eff.for.low = matrix(0,nr,nbase)
eff.rev.low = matrix(0,nr,nbase)
base.for.low = matrix(0,nr,nbase)
base.rev.low = matrix(0,nr,nbase)

ind.high = ccount.train > 25
ind.low = ccount.train <= 25

dprof.train.high = dprof.train[ind.high, ]
dprof.train.low = dprof.train[ind.low, ]
ccount.train.high = ccount.train[ind.high]
ccount.train.low = ccount.train[ind.low]


for(i in 1:nr){
  ind = sample(1:dim(dprof.train.high)[1],1000)
  est.for = multiseq(as.matrix(dprof.train.high[ind,1:nbase]),log(ccount.train.high[ind]+1),lm.approx = TRUE)
  #est.rev = multiseq(as.matrix(dprof.train.high[ind,(nbase + 1):(2 * nbase)]),log(ccount.train.high[ind]+1),lm.approx = TRUE)
  base.for.high[i,] = est.for$baseline.mean
  #base.rev[i,] = est.rev$baseline.mean
  eff.for.high[i,] = est.for$effect.mean
  #eff.rev[i,] = est.rev$effect.mean
  print(i)
}

base.mean.for.high = colMeans(base.for.high)
#base.mean.rev.high = colMeans(base.rev.high)
eff.mean.for.high = colMeans(eff.for.high)
#eff.mean.rev.high = colMeans(eff.rev.high)

for(i in 1:nr){
  ind = sample(1:dim(dprof.train.low)[1],1000)
  est.for = multiseq(as.matrix(dprof.train.low[ind,1:nbase]),log(ccount.train.low[ind]+1),lm.approx = TRUE)
  #est.rev = multiseq(as.matrix(dprof.train.low[ind,(nbase + 1):(2 * nbase)]),log(ccount.train.low[ind]+1),lm.approx = TRUE)
  base.for.low[i,] = est.for$baseline.mean
  #base.rev[i,] = est.rev$baseline.mean
  eff.for.low[i,] = est.for$effect.meanx
  #eff.rev[i,] = est.rev$effect.mean
  print(i)
}

base.mean.for.low = colMeans(base.for.low)
#base.mean.rev.low = colMeans(base.rev.low)
eff.mean.for.low = colMeans(eff.for.low)
#eff.mean.rev.low = colMeans(eff.rev.low)







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
lambda.for.high = exp(rep(1,length(ccount.prior.val))%o%base.mean.for.high + ccount.prior.val%o%eff.mean.for.high)
lambda.for.low = exp(rep(1,length(ccount.prior.val))%o%base.mean.for.low + ccount.prior.val%o%eff.mean.for.low)
lambda.for = rbind(lambda.for.low[ccount.prior.val<=log(25),], lambda.for.high[ccount.prior.val>log(25),])

#lambda.rev = exp(rep(1,length(ccount.prior.val))%o%base.mean.rev + ccount.prior.val%o%eff.mean.rev)
lik.for = matrix(0,length(ccount.test),length(ccount.prior.val))
#lik.rev = matrix(0,length(ccount.test),length(ccount.prior.val))
for(i in 1:length(ccount.test)){
  loglik.ini.for = rowSums(t(apply(lambda.for,1,dpois,x = as.numeric(dprof.test[i,1:nbase]),log = TRUE)))
  loglik.ini.for = loglik.ini.for - max(loglik.ini.for)
#  loglik.ini.rev = rowSums(t(apply(lambda.rev,1,dpois,x = as.numeric(dprof.test[i,(nbase + 1):(2 * nbase)]),log = TRUE)))
#  loglik.ini.rev = loglik.ini.rev - max(loglik.ini.rev)
  lik.for[i,] = exp(loglik.ini.for)
#  lik.rev[i,] = exp(loglik.ini.rev)
}

##computing the posterior
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

