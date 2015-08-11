eff.mean.for.low.mod = (base.mean.for.high + 3.25 * eff.mean.for.high - base.mean.for.low) / 3.25
lambda.for.high = exp(rep(1,length(ccount.prior.val))%o%base.mean.for.high + ccount.prior.val%o%eff.mean.for.high)
lambda.for.low = exp(rep(1,length(ccount.prior.val))%o%base.mean.for.low + ccount.prior.val%o%eff.mean.for.low.mod)
lambda.for = rbind(lambda.for.low[ccount.prior.val<=log(25),], lambda.for.high[ccount.prior.val>log(25),])
plot(lambda.for[18,],type='l')
lines(lambda.for[19,],co=2)
lines(lambda.for[20,],co=3)
lines(lambda.for[21,],co=4)
lines(lambda.for[17,],co=5)
lines(lambda.for[16,],co=6)
lines(lambda.for[15,],co=7)
lines(lambda.for[22,],co=8)
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
