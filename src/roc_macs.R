library(CENTIPEDE)
library(pROC)

peak.info = read.table("data/macs_peaks/Ctcf_q1pc_peaks.bed")
peak.ind = 0
nbase = 2^10
for(i in 1:dim(cinfo)[1]){
  peak.ind[i] = sum((as.character(peak.info[, 1]) == as.character(cinfo[i, 1])) & (peak.info[, 2] < cinfo[i, 2]) & (peak.info[, 3] > cinfo[i, 3])) > 0
}


unbound.ind = !(peak.ind[!train.ind]) & ccount.test < median(ccount)
bound.ind = as.logical(peak.ind[!train.ind])

boxplot(log(ccount.test + 1)[unbound.ind], log(ccount.test + 1)[bound.ind])


#roc.res.ms = rocplot(ccount.post.logmean[unbound.ind], ccount.post.logmean[bound.ind])
#plot(roc.res.ms$fpr,roc.res.ms$tpr, type = 'l', xlim = c(0, 1))

#roc.res.ms = roc(bound.ind ~ ccount.post.logmean)
roc.res.ms = roc(bound.ind ~ ccount.post.logmean.list[[1]])

plot(roc.res.ms)

pscore = cinfo[,5]
pscore.test = pscore[!train.ind]


centFit <- fitCentipede(Xlist = list(DNase=as.matrix(dprof.test[, c(413:612, 1437:1636)])), Y=matrix(rep(1,dim(dprof.test)[1], nc = 1)), sweeps = 300)
centFit.pwm <- fitCentipede(Xlist = list(DNase=as.matrix(dprof.test)[, c(413:612, 1437:1636)]), Y=cbind(rep(1,dim(dprof.test)[1]),pscore.test))
#roc.res.cent = rocplot(centFit$PostPr[unbound.ind], centFit$PostPr[bound.ind])
#roc.res.cent.pwm = rocplot(centFit.pwm$PostPr[unbound.ind], centFit.pwm$PostPr[bound.ind])

#lines(roc.res.cent$fpr, roc.res.cent$tpr, col = 2)
#lines(roc.res.cent.pwm$fpr, roc.res.cent.pwm$tpr, col = 3)

roc.res.cent = roc(bound.ind ~ as.vector(centFit$PostPr))
roc.res.cent.pwm = roc(bound.ind ~ as.vector(centFit.pwm$PostPr))

lines(roc.res.cent, col = 2)
lines(roc.res.cent.pwm, col = 3)

dcut = rowSums(dprof.test[, c(413:612, 1437:1636)])
#roc.res.dcut = rocplot(dcut[unbound.ind], dcut[bound.ind])

#lines(roc.res.dcut$fpr,roc.res.dcut$tpr, col = 4)

roc.res.dcut = roc(bound.ind ~ dcut)

lines(roc.res.dcut, col = 4)

##find differences between centipede and ms
diff.ind = centFit.pwm$PostPr > 0.95 & ccount.post.logmean < 2 & bound.ind == 1
diff.ind = which(diff.ind == 1)

ccount.test[diff.ind[1]]

centFit.pwm$PriorLogRatio[diff.ind[1]]
centFit.pwm$MultiNomLogRatio[diff.ind[1]]
centFit.pwm$NegBinLogRatio[diff.ind[1]]

as.numeric(apply(as.matrix(dprof.test[diff.ind[1], base.ind.for[[4]]]), 1, dmultinom, size = NULL, prob = lambda.for[37, ], log = TRUE) - apply(as.matrix(dprof.test[diff.ind[1], base.ind.for[[4]]]), 1, dmultinom, size = NULL, prob = lambda.for[5, ], log = TRUE))
as.numeric(dpois(rowSums(as.matrix(dprof.test[diff.ind[1], base.ind.for[[4]]])), lambda = sum(lambda.for[37, ]), log = TRUE) - dpois(rowSums(as.matrix(dprof.test[diff.ind[1], base.ind.for[[4]]])), lambda = sum(lambda.for[5, ]), log = TRUE))

#all sites, no multiseq

unbound.ind = !(peak.ind) & ccount < median(ccount)
bound.ind = as.logical(peak.ind)

boxplot(log(ccount + 1)[unbound.ind], log(ccount + 1)[bound.ind])



pscore = cinfo[,5]

#c(413:612, 1437:1636)

centFit <- fitCentipede(Xlist = list(DNase=as.matrix(dprof[, ])), Y=matrix(rep(1,dim(dprof)[1], nc = 1)), sweeps = 300)
centFit.pwm <- fitCentipede(Xlist = list(DNase=as.matrix(dprof)[, ]), Y=cbind(rep(1,dim(dprof)[1]),pscore),sweeps = 500)
#roc.res.cent = rocplot(centFit$PostPr[unbound.ind], centFit$PostPr[bound.ind])
#roc.res.cent.pwm = rocplot(centFit.pwm$PostPr[unbound.ind], centFit.pwm$PostPr[bound.ind])

#lines(roc.res.cent$fpr, roc.res.cent$tpr, col = 2)
#lines(roc.res.cent.pwm$fpr, roc.res.cent.pwm$tpr, col = 3)

roc.res.cent = roc(bound.ind ~ as.vector(centFit$PostPr))
roc.res.cent.pwm = roc(bound.ind ~ as.vector(centFit.pwm$PostPr))

plot(roc.res.cent, col = 2)
lines(roc.res.cent.pwm, col = 3)

dcut = rowSums(dprof[, ])
#roc.res.dcut = rocplot(dcut[unbound.ind], dcut[bound.ind])

#lines(roc.res.dcut$fpr,roc.res.dcut$tpr, col = 4)

roc.res.dcut = roc(bound.ind ~ dcut)

lines(roc.res.dcut, col = 4)

