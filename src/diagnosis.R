nbase=2^13

dprof.int.mean=matrix(0,8,nbase)

dprof.int.mean[1,]=colMeans(as.matrix(dprof[log(ccount+1)<1,1:nbase]))
dprof.int.mean[8,]=colMeans(as.matrix(dprof[log(ccount+1)>7,1:nbase]))


for(i in 2:7){
  dprof.int.mean[i,]=colMeans(as.matrix(dprof[log(ccount+1)>(i-1)&log(ccount+1)<=i,1:nbase]))
}


plot(log(ma(dprof.int.mean[7,3585:4608],10)),type='l',ylim=c(-6,0))
for(i in 6:1){
  lines(log(ma(dprof.int.mean[i,3585:4608],10)),col=9-i)
}



plot(log(ma(dprof.int.mean[7,1:8192],100)),type='l',ylim=c(-6,0))
for(i in 6:1){
  lines(log(ma(dprof.int.mean[i,1:8192],100)),col=9-i)
}


plot(log(dprof.int.mean[7,3585:4608]),type='l',ylim=c(-6,0))
for(i in 6:1){
  lines(log(dprof.int.mean[i,3585:4608]),col=9-i)
}
