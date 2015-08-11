library(multiseq)

nbase=2^14
nsample=5000

setwd("~/projects/multiseq_analysis")

file.names=list.files("data/")
chip.file.names=file.names[seq(3,247,4)]
dnase.file.names=file.names[seq(2,246,4)]

for(j in 1:length(chip.file.names)){
  factor.name=paste(strsplit(chip.file.names[j],"_")[[1]][1],strsplit(chip.file.names[j],"_")[[1]][2],sep="_")
  try.test=try(dprof<-read.table(file.path("data",dnase.file.names[j])))
  if(class(try.test)=="try-error") next
  cinfo=read.table(file.path("data",chip.file.names[j]),header=TRUE)
  ccount=cinfo[,6]

  eff.for=matrix(0,10,nbase)
  eff.rev=matrix(0,10,nbase)

  for(i in 1:10){
    ind=sample(1:dim(dprof)[1],min(round(dim(dprof)[1]/3),nsample))
    eff.for[i,]=multiseq(as.matrix(dprof[ind,1:nbase]),log(ccount[ind]+1),lm.approx=TRUE)$effect.mean
    eff.rev[i,]=multiseq(as.matrix(dprof[ind,(nbase+1):(2*nbase)]),log(ccount[ind]+1),lm.approx=TRUE)$effect.mean
  }

#  eff.nl.for=matrix(0,10,nbase)
#  eff.nl.rev=matrix(0,10,nbase)
#
#  for(i in 1:10){
#    ind=sample(1:dim(dprof)[1],min(round(dim(dprof)[1]/3),nsample))
#    eff.nl.for[i,]=multiseq(as.matrix(dprof[ind,1:nbase]),ccount[ind],lm.approx=TRUE)$effect.mean
#    eff.nl.rev[i,]=multiseq(as.matrix(dprof[ind,(nbase+1):(2*nbase)]),ccount[ind],lm.approx=TRUE)$effect.mean
#  }

  cmef=colMeans(eff.for)
  cmer=colMeans(eff.rev)

  if(is.na(max(c(cmef,cmer)))|!is.finite(max(c(cmef,cmer)))) next

  pdf(paste(paste0("results/plots/",factor.name),"multiseq","est","log.pdf",sep="_"),width=12,height=7)
  plot(cmef,type='l',ylim=c(min(c(cmef,cmer))-0.05,max(c(cmef,cmer))+0.05),xlab="base pair",ylab="log(effect)",main="Effect estimated from multiseq, based on log(ChIP-seq counts)")
  lines(cmer,col=2,lty=2)
  dev.off()


  bound=quantile(ccount,c(0.9,0.99))

  dprof.high.s1=as.matrix(dprof[ccount>bound[2],1:nbase])
  dprof.med.s1=as.matrix(dprof[ccount>bound[1]&ccount<=bound[2],1:nbase])
  dprof.low.s1=as.matrix(dprof[ccount<=bound[1],1:nbase])

  dprof.high.s2=as.matrix(dprof[ccount>bound[2],(nbase+1):(2*nbase)])
  dprof.med.s2=as.matrix(dprof[ccount>bound[1]&ccount<=bound[2],(nbase+1):(2*nbase)])
  dprof.low.s2=as.matrix(dprof[ccount<=bound[1],(nbase+1):(2*nbase)])

  cmh1=colMeans(dprof.high.s1)
  cmh2=colMeans(dprof.high.s2)
  cmm1=colMeans(dprof.med.s1)
  cmm2=colMeans(dprof.med.s2)
  cml1=colMeans(dprof.low.s1)
  cml2=colMeans(dprof.low.s2)


  pdf(paste(paste0("results/plots/",factor.name),"group","average.pdf",sep="_"))
  par(mfrow=c(2,1))
  plot(cmh1,ylim=c(0,max(cmh1)+0.1),type='l',xlab="base pair",ylab="reads",main="Average DNase profile for 1st strand, grouped by ChIP-seq counts")
  lines(cmm1,col=2)
  lines(cml1,col=3)

  plot(cmh2,ylim=c(0,max(cmh2)+0.1),type='l',xlab="base pair",ylab="reads",main="Average DNase profile for 2nd strand, grouped by ChIP-seq counts")
  lines(cmm2,col=2)
  lines(cml2,col=3)
  dev.off()

  pdf(paste(paste0("results/plots/",factor.name),"group","high.pdf",sep="_"),width=12,height=7)
  plot(cmh1,type='l',ylim=c(0,max(c(cmh1,cmh2))+0.1),xlab="base pair",ylab="reads",main="Average DNase profile for both strands across sites with high ChIP-seq counts")
  lines(cmh2,col=2,lty=2)
  dev.off()

  print(j)
}