library(multiseq)

nbase = 2^10

normalize = function(x) x/sum(x)

dprof = read.table("data/Ebf_S79_dnase_data.txt.gz")
cinfo = read.table("data/Ebf_S79_site_and_chip_data.txt.gz",header = TRUE)
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
	  ind = sample(1:dim(dprof.train)[1],3000)
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

save(base.mean.for, base.mean.rev, eff.mean.for, eff.mean.rev, ccount.train,  dprof.test, ccount.test, file = "results/res_ms_ebf.Robj")
