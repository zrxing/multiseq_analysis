library(multiseq)

nbase = 2^13

normalize = function(x) x/sum(x)

dprof = read.table("data/Ctcf_S2_dnase_data.txt.gz")
cinfo = read.table("data/Ctcf_S2_site_and_chip_data.txt.gz",header = TRUE)
ccount = cinfo[,6]

#ind.thresh=(rowSums(dprof)>0)&(rowSums(dprof)<2000)
#dprof=dprof[ind.thresh,]
#ccount=ccount[ind.thresh]


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
  ind = sample(1:dim(dprof.train)[1],5000)
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

prior.breaks = log(quantile(ccount.train,c(0,0.5,0.7,0.9,0.99,1))+1)
prior.breaks[length(prior.breaks)] = ceiling(prior.breaks[length(prior.breaks)])
names(prior.breaks) = NULL
ccount.prior.info = hist(log(ccount.train+1),breaks = prior.breaks, plot = FALSE)
#ccount.prior.val = ccount.prior.info$mids
#ccount.prior.val = prior.breaks[1:(length(prior.breaks)-1)]
lccount.train = log(ccount.train+1)
ccount.prior.val = 0
for(i in 1:(length(prior.breaks)-1)){
  ccount.prior.val[i] = mean(lccount.train[lccount.train>=prior.breaks[i]&lccount.train<prior.breaks[i+1]])
}
#ccount.prior.prob = normalize(smooth.spline(ccount.prior.info$mids,ccount.prior.info$counts)$y)
ccount.prior.prob = normalize(ccount.prior.info$counts)

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
ccount.post.class.for = 0
ccount.post.mean.rev = 0
ccount.post.class.rev = 0
for(i in 1:length(ccount.test)){
  ccount.post.class.for[i] = which(ccount.post.prob.for[i,]==max(ccount.post.prob.for[i,]))
  ccount.post.mean.for[i] = sum((ccount.post.val)*ccount.post.prob.for[i,])
  ccount.post.class.rev[i] = which(ccount.post.prob.rev[i,]==max(ccount.post.prob.rev[i,]))
  ccount.post.mean.rev[i] = sum((ccount.post.val)*ccount.post.prob.rev[i,])
}


ccount.test.m = log(ccount.test+1)
ccount.test.class = 0
for(i in 1:(length(prior.breaks)-1)){
  ccount.test.class = ccount.test.class + i*(ccount.test.m>=prior.breaks[i]&ccount.test.m<prior.breaks[i+1])
}

ccount.post.meanclass.for = 0
for(i in 1:(length(prior.breaks)-1)){
  ccount.post.meanclass.for = ccount.post.meanclass.for + i*(ccount.post.mean.for>=prior.breaks[i]&ccount.post.mean.for<prior.breaks[i+1])
}



res.multiclass=sum(ccount.post.class.for==ccount.test.class)/length(ccount.test.class)
sum(ccount.post.meanclass.for==ccount.test.class)/length(ccount.test.class)
