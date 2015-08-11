load("results/prediction_log_val_large.RData")


unbound.ind = ccount.test <= median(ccount.test)
bound.ind = ccount.test > quantile(ccount.test, probs = 0.95)


dcut = rowSums(dprof.test[, 3997:4196])
roc.res.dcut = rocplot(dcut[unbound.ind], dcut[bound.ind])

plot(roc.res.dcut$fpr,roc.res.dcut$tpr, type='l')


dcut = rowSums(dprof.test[, 3897:4296])
roc.res.dcut = rocplot(dcut[unbound.ind], dcut[bound.ind])

lines(roc.res.dcut$fpr,roc.res.dcut$tpr, col=2)

dcut = rowSums(dprof.test[, 3597:4596])
roc.res.dcut = rocplot(dcut[unbound.ind], dcut[bound.ind])

lines(roc.res.dcut$fpr,roc.res.dcut$tpr, col=3)

dcut = rowSums(dprof.test[, 3097:5096])
roc.res.dcut = rocplot(dcut[unbound.ind], dcut[bound.ind])

lines(roc.res.dcut$fpr,roc.res.dcut$tpr, col=4)

dcut = rowSums(dprof.test[, 2097:6096])
roc.res.dcut = rocplot(dcut[unbound.ind], dcut[bound.ind])

lines(roc.res.dcut$fpr,roc.res.dcut$tpr, col=5)

dcut = rowSums(dprof.test)
roc.res.dcut = rocplot(dcut[unbound.ind], dcut[bound.ind])

lines(roc.res.dcut$fpr,roc.res.dcut$tpr, col=6)


load("results/test_region_size_ms.Robj")

plot(roc.res.ms.13.50$fpr,roc.res.ms.13.50$tpr, type = 'l', xlim = c(0, 1))
lines(roc.res.ms.10.50$fpr,roc.res.ms.10.50$tpr, col=2)

plot(roc.res.ms.13.95$fpr,roc.res.ms.13.95$tpr, type = 'l', xlim = c(0, 1))
lines(roc.res.ms.10.95$fpr,roc.res.ms.10.95$tpr, col=2)
