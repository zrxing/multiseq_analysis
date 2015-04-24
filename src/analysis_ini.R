library(multiseq)


dprof=read.table("data/Ctcf_S2_dnase_data.txt.gz")
cinfo=read.table("data/Ctcf_S2_site_and_chip_data.txt.gz",header=TRUE)
ccount=cinfo[,6]

ccount.nz=ccount[ccount>10]
hist((log(ccount.nz)),freq=FALSE,breaks=100)
lines(dexp(0:4000,1/165),col=2)

qqplot(ccount.nz,rpois(10000,165))





eff.for=matrix(0,20,1024)
eff.rev=matrix(0,20,1024)

for(i in 1:20){
  ind=sample(1:dim(dprof)[1],5000)
  eff.for[i,]=multiseq(as.matrix(dprof[ind,1:1024]),log(ccount[ind]+1),lm.approx=TRUE)$effect.mean
  eff.rev[i,]=multiseq(as.matrix(dprof[ind,1025:2048]),log(ccount[ind]+1),lm.approx=TRUE)$effect.mean
  print(i)
}

eff.nl.for=matrix(0,20,1024)
eff.nl.rev=matrix(0,20,1024)

for(i in 1:20){
  ind=sample(1:dim(dprof)[1],5000)
  eff.nl.for[i,]=multiseq(as.matrix(dprof[ind,1:1024]),ccount[ind],lm.approx=TRUE)$effect.mean
  eff.nl.rev[i,]=multiseq(as.matrix(dprof[ind,1025:2048]),ccount[ind],lm.approx=TRUE)$effect.mean
  print(i)
}


pdf("ctcf_multiseq_est_log.pdf",width=12,height=7)
plot(colMeans(eff.for),type='l',ylim=c(-0.02,1),xlab="base pair",ylab="log(effect)",main="Effect estimated from multiseq, based on log(ChIP-seq counts)")
lines(colMeans(eff.rev),col=2)
dev.off()

plot(colMeans(eff.nl.for),type='l',ylim=c(0,0.01),xlab="base pair",ylab="log(effect)",main="Effect estimated from multiseq, based on ChIP-seq counts")
lines(colMeans(eff.nl.rev),col=2)




dprof.high.s1=as.matrix(dprof[ccount>200,1:1024])
dprof.med.s1=as.matrix(dprof[ccount>10&ccount<=200,1:1024])
dprof.low.s1=as.matrix(dprof[ccount<=10,1:1024])

dprof.high.s2=as.matrix(dprof[ccount>200,1025:2048])
dprof.med.s2=as.matrix(dprof[ccount>10&ccount<=200,1025:2048])
dprof.low.s2=as.matrix(dprof[ccount<=10,1025:2048])


#dprof.high.s1[dprof.high.s1>10]=NA
#dprof.med.s1[dprof.high.s1>10]=NA
#dprof.low.s1[dprof.high.s1>10]=NA
#
#dprof.high.s2[dprof.high.s2>10]=NA
#dprof.med.s2[dprof.high.s2>10]=NA
#dprof.low.s2[dprof.high.s2>10]=NA


pdf("plots/ctcf_group_average.pdf")
par(mfrow=c(2,1))
plot(colMeans(dprof.high.s1,na.rm=TRUE),type='l',xlab="base pair",ylab="reads",main="Average DNase profile for 1st strand, grouped by ChIP-seq counts")
lines(colMeans(dprof.med.s1,na.rm=TRUE),col=2)
lines(colMeans(dprof.low.s1,na.rm=TRUE),col=3)

plot(apply(dprof.high.s2,2,mean,na.rm=TRUE),type='l',xlab="base pair",ylab="reads",main="Average DNase profile for 2nd strand, grouped by ChIP-seq counts")
lines(colMeans(dprof.med.s2,na.rm=TRUE),col=2)
lines(colMeans(dprof.low.s2,na.rm=TRUE),col=3)
dev.off()

pdf("plots/ctcf_group_high.pdf",width=12,height=7)
plot(colMeans(dprof.high.s1),type='l',ylim=c(0,1.2),xlab="base pair",ylab="reads",main="Average DNase profile for both strands across sites with high ChIP-seq counts")
lines(colMeans(dprof.high.s2),col=2)
abline(v=504,lty=2)
abline(v=533,lty=2)
dev.off()

which(colMeans(dprof.high.s2[,1:520])==sort(colMeans(dprof.high.s2[,1:520]),decreasing=TRUE)[1])
#504
which(colMeans(dprof.high.s2[,521:1024])==sort(colMeans(dprof.high.s2[,521:1024]),decreasing=TRUE)[1])+520
#533

