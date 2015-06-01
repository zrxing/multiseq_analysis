library(multiseq)

nbase = 2^10

normalize = function(x) x/sum(x)

dprof = read.table("data/Ctcf_S2_dnase_data.txt.gz")
cinfo = read.table("data/Ctcf_S2_site_and_chip_data.txt.gz",header = TRUE)
ccount = cinfo[,6]

train.size = round(0.8*dim(dprof)[1])
set.seed(417)
train.ind = sample(1:dim(dprof)[1],train.size,replace = FALSE)
train.ind = (1:dim(dprof)[1])%in%train.ind

dprof.train = dprof[train.ind,]
dprof.test = dprof[!train.ind,]

ccount.train = ccount[train.ind]
ccount.test = ccount[!train.ind]


##estimating mean intensity
nr = 10
eff.for = matrix(0,nr,nbase)
eff.rev = matrix(0,nr,nbase)
base.for = matrix(0,nr,nbase)
base.rev = matrix(0,nr,nbase)

for(i in 1:nr){
  ind = sample(1:dim(dprof.train)[1],1000)
  est.for = multiseq(as.matrix(dprof.train[ind,1:nbase]),log(ccount.train[ind]+1),lm.approx = TRUE)
  est.rev = multiseq(as.matrix(dprof.train[ind,(nbase + 1):(2 * nbase)]),log(ccount.train[ind]+1),lm.approx = TRUE)
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
  loglik.ini.for = rowSums(t(apply(lambda.for,1,dpois,x = as.numeric(dprof.test[i,1:nbase]),log = TRUE)))
  loglik.ini.for = loglik.ini.for - max(loglik.ini.for)
  loglik.ini.rev = rowSums(t(apply(lambda.rev,1,dpois,x = as.numeric(dprof.test[i,(nbase + 1):(2 * nbase)]),log = TRUE)))
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

plot(ccount.post.mean.for, ccount.test, xlim=c(0,2000), ylim=c(0,1000))
abline(0, 1, col = 2)
plot(ccount.post.mean.for,ccount.test)
abline(0, 1, col = 2)

plot(ccount.post.logmean.for,log(ccount.test+1), xlim = c(0, 9), ylim = c(0, 9))
abline(0, 1, col = 2)

cor(ccount.post.logmean.for,log(ccount.test+1))
