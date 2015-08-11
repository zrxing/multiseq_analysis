library(multiseq)

nbase = 2^13
nsample = 3000

dir.name = "~/projects/multiseq_analysis/"
setwd(dir.name)

normalize = function(x) x/sum(x)

dprof = read.table("data/Ctcf_S2_dnase_data_large.txt.gz", fill = TRUE)
dprof[is.na(dprof)] = 0
ind.for = (1/2 * nbase + 1):(3/2 * nbase)
ind.rev = 2 * nbase + (1/2 * nbase + 1):(3/2 * nbase)
dprof = dprof[, c(ind.for, ind.rev)]
cinfo = read.table("data/Ctcf_S2_site_and_chip_data.txt.gz",header = TRUE)
ccount = cinfo[,6]

train.size = round(0.8*dim(dprof)[1])
set.seed(527)
train.ind = sample(1:dim(dprof)[1],train.size,replace = FALSE)
train.ind = (1:dim(dprof)[1])%in%train.ind

dprof.train = dprof[train.ind,]
dprof.test = dprof[!train.ind,]

ccount.train = ccount[train.ind]
ccount.test = ccount[!train.ind]

rm(dprof)

##estimating mean intensity
nr = 20
eff.for = matrix(0,nr,nbase)
eff.rev = matrix(0,nr,nbase)
base.for = matrix(0,nr,nbase)
base.rev = matrix(0,nr,nbase)

for(i in 1:nr){
  ind = sample(1:dim(dprof.train)[1],min(round(dim(dprof.train)[1]/3), nsample))
  est.for = multiseq(as.matrix(dprof.train[ind, 1:nbase]), log(ccount.train[ind]+1), lm.approx = TRUE)
  est.rev = multiseq(as.matrix(dprof.train[ind, (nbase + 1):(2 * nbase)]), log(ccount.train[ind]+1), lm.approx = TRUE)
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

save.image(file.path(dir.name,"results/ctcf_large","ctcf_large.RData"))