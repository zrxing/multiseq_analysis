library(CENTIPEDE)
source("src/rocplot.R")


unbound.ind = ccount.test <= median(ccount.test)
bound.ind = ccount.test > quantile(ccount.test, probs = 0.5)

unbound.pred = ccount.post.logmean.for[unbound.ind] <= median(ccount.post.logmean.for)
bound.pred = ccount.post.logmean.for[bound.ind] > quantile(ccount.post.logmean.for, probs = 0.95)

mean(unbound.pred)
mean(bound.pred)


# roc.res.ms = rocplot(ccount.post.logmean.for[unbound.ind], ccount.post.logmean.for[bound.ind])
# 
# plot(roc.res.ms$fpr,roc.res.ms$tpr, type = 'l', xlim = c(0, 1))
# 
# roc.res.ms = rocplot(ccount.post.logmean.for.quant.sep[unbound.ind], ccount.post.logmean.for.quant.sep[bound.ind])
# 
# lines(roc.res.ms$fpr,roc.res.ms$tpr, col = 2)
# 
# roc.res.ms = rocplot(ccount.post.logmean.for.quant.sep.bin[unbound.ind], ccount.post.logmean.for.quant.sep.bin[bound.ind])
# 
# lines(roc.res.ms$fpr,roc.res.ms$tpr, col = 3)
# 
# roc.res.ms = rocplot(ccount.post.logmean.for.quant.common.base[unbound.ind], ccount.post.logmean.for.quant.common.base[bound.ind])
# 
# lines(roc.res.ms$fpr,roc.res.ms$tpr, col = 4)
# 
# roc.res.ms = rocplot(ccount.post.logmean.for.quant.single.base[unbound.ind], ccount.post.logmean.for.quant.single.base[bound.ind])
# 
# lines(roc.res.ms$fpr,roc.res.ms$tpr, col = 5)
# 
# roc.res.ms = rocplot(ccount.post.logmean.for.quant.sep.base[unbound.ind], ccount.post.logmean.for.quant.sep.base[bound.ind])
# 
# lines(roc.res.ms$fpr,roc.res.ms$tpr, col = 6)
# 
# 
# pscore = cinfo[,5]
# pscore.test = pscore[!train.ind]
# 
# 
# centFit <- fitCentipede(Xlist = list(DNase=as.matrix(dprof.test[, c(413:612, 1437:1636)])), Y=matrix(rep(1,dim(dprof.test)[1], nc = 1)))
# centFit.pwm <- fitCentipede(Xlist = list(DNase=as.matrix(dprof.test)[, c(413:612, 1437:1636)]), Y=cbind(rep(1,dim(dprof.test)[1]),pscore.test))
# roc.res.cent = rocplot(centFit$PostPr[unbound.ind], centFit$PostPr[bound.ind])
# roc.res.cent.pwm = rocplot(centFit.pwm$PostPr[unbound.ind], centFit.pwm$PostPr[bound.ind])
# 
# lines(roc.res.cent$fpr,roc.res.cent$tpr, col = 2)
# lines(roc.res.cent.pwm$fpr,roc.res.cent.pwm$tpr, col = 3)
# 
# dcut = rowSums(dprof.test[, 413:612])
# roc.res.dcut = rocplot(dcut[unbound.ind], dcut[bound.ind])
# 
# lines(roc.res.dcut$fpr,roc.res.dcut$tpr, col = 4)
# 




roc.res.ms = rocplot(ccount.post.logmean.for[unbound.ind], ccount.post.logmean.for[bound.ind])

plot(roc.res.ms$fpr,roc.res.ms$tpr, type = 'l', xlim = c(0, 1))

roc.res.ms = rocplot(ccount.post.logmean.for.quant.sep[unbound.ind], ccount.post.logmean.for.quant.sep[bound.ind])

lines(roc.res.ms$fpr,roc.res.ms$tpr, col = 2)

roc.res.ms = rocplot(ccount.post.logmean.for.quant.sep.bin[unbound.ind], ccount.post.logmean.for.quant.sep.bin[bound.ind])

lines(roc.res.ms$fpr,roc.res.ms$tpr, col = 3)

roc.res.ms = rocplot(ccount.post.logmean.for.quant.common.base[unbound.ind], ccount.post.logmean.for.quant.common.base[bound.ind])

lines(roc.res.ms$fpr,roc.res.ms$tpr, col = 4)

roc.res.ms = rocplot(ccount.post.logmean.for.quant.single.base[unbound.ind], ccount.post.logmean.for.quant.single.base[bound.ind])

lines(roc.res.ms$fpr,roc.res.ms$tpr, col = 5)

roc.res.ms = rocplot(ccount.post.logmean.for.quant.sep.base[unbound.ind], ccount.post.logmean.for.quant.sep.base[bound.ind])

lines(roc.res.ms$fpr,roc.res.ms$tpr, col = 6)


pscore = cinfo[,5]
pscore.test = pscore[!train.ind]


centFit <- fitCentipede(Xlist = list(DNase=as.matrix(dprof.test[, c(413:612, 1437:1636)])), Y=matrix(rep(1,dim(dprof.test)[1], nc = 1)))
centFit.pwm <- fitCentipede(Xlist = list(DNase=as.matrix(dprof.test)[, c(413:612, 1437:1636)]), Y=cbind(rep(1,dim(dprof.test)[1]),pscore.test))
roc.res.cent = rocplot(centFit$PostPr[unbound.ind], centFit$PostPr[bound.ind])
roc.res.cent.pwm = rocplot(centFit.pwm$PostPr[unbound.ind], centFit.pwm$PostPr[bound.ind])

lines(roc.res.cent$fpr,roc.res.cent$tpr, col = 2)
lines(roc.res.cent.pwm$fpr,roc.res.cent.pwm$tpr, col = 3)

dcut = rowSums(dprof.test[, c(413:612, 1437:1636)])
roc.res.dcut = rocplot(dcut[unbound.ind], dcut[bound.ind])

lines(roc.res.dcut$fpr,roc.res.dcut$tpr, col = 4)






#centipede plots
plot(centFit.pwm$PostPr[bound.ind], log(ccount.test[bound.ind]+1),ylim=c(0,9))
points(centFit.pwm$PostPr[unbound.ind], log(ccount.test[unbound.ind]+1),col=2)

plot(log(dcut+1),log(ccount.test+1))
abline(h=log(quantile(ccount.test,0.95)+1))
abline(h=log(quantile(ccount.test,0.5)+1))
