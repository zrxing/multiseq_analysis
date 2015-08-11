library(CENTIPEDE)
load("results/prediction_log_val.RData")

centFit.pwm <- fitCentipede(Xlist = list(DNase=as.matrix(dprof.test[, c(413:612, 1437:1636)])), Y=cbind(rep(1,dim(dprof.test)[1]),pscore.test))

plot(as.numeric(dprof.test[bound.ind, ][1417, 1025:2048]))
centFit.pwm$PostPr[bound.ind][1417]
centFit.pwm$LogRatios[bound.ind][1417]
centFit.pwm$PriorLogRatio[bound.ind][1417]
centFit.pwm$NegBinLogRatio[bound.ind][1417]
centFit.pwm$MultiNomLogRatio[bound.ind][1417]

plotProfile(centFit.pwm$LambdaParList[[1]])


plot(ashsmooth.pois(as.numeric(dprof.test[bound.ind, ][1417, 1:1024])),type='l')
