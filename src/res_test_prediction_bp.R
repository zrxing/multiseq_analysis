source("src/rocplot.R")

par(mfrow = c(2,2))

plot(ccount.post.logmean.for.list[[1]],log(ccount.test+1), xlim = c(0, 9), ylim = c(0, 9))
abline(0, 1, col = 2)

plot(ccount.post.logmean.for.list[[2]],log(ccount.test+1), xlim = c(0, 9), ylim = c(0, 9))
abline(0, 1, col = 2)

plot(ccount.post.logmean.for.list[[3]],log(ccount.test+1), xlim = c(0, 9), ylim = c(0, 9))
abline(0, 1, col = 2)

plot(ccount.post.logmean.for.list[[4]],log(ccount.test+1), xlim = c(0, 9), ylim = c(0, 9))
abline(0, 1, col = 2)



cor(ccount.post.logmean.for.list[[1]],log(ccount.test+1))
cor(ccount.post.logmean.for.list[[2]],log(ccount.test+1))
cor(ccount.post.logmean.for.list[[3]],log(ccount.test+1))
cor(ccount.post.logmean.for.list[[4]],log(ccount.test+1))



roc.res.ms = rocplot(ccount.post.logmean.for.list[[1]][unbound.ind], ccount.post.logmean.for.list[[1]][bound.ind])
plot(roc.res.ms$fpr,roc.res.ms$tpr, type = 'l', xlim = c(0, 1))

roc.res.ms = rocplot(ccount.post.logmean.for.list[[2]][unbound.ind], ccount.post.logmean.for.list[[2]][bound.ind])
lines(roc.res.ms$fpr,roc.res.ms$tpr, col = 2)

roc.res.ms = rocplot(ccount.post.logmean.for.list[[3]][unbound.ind], ccount.post.logmean.for.list[[3]][bound.ind])
lines(roc.res.ms$fpr,roc.res.ms$tpr, col = 3)

roc.res.ms = rocplot(ccount.post.logmean.for.list[[4]][unbound.ind], ccount.post.logmean.for.list[[4]][bound.ind])
lines(roc.res.ms$fpr,roc.res.ms$tpr, col = 4)
