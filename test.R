setwd("~/projects/multiseq_analysis")

dprof=read.table("data/Ctcf_S2_dnase_data.txt.gz")
ccount=read.table("data/Ctcf_S2_site_and_chip_data.txt.gz",header=TRUE)

library(multiseq)

eff.nl.lm=matrix(0,20,1024)
eff.nl.glm=matrix(0,20,1024)
for(i in 1:20){
ind=sample(1:dim(dprof)[1],10000)

#ind=1:1000

res.nl.lm=multiseq(as.matrix(dprof[ind,1:1024]),ccount[ind,6],lm.approx=TRUE)
res.nl.glm=multiseq(as.matrix(dprof[ind,1:1024]),ccount[ind,6],lm.approx=FALSE)
eff.nl.lm[i,]=res.nl.lm$effect.mean
eff.nl.glm[i,]=res.nl.glm$effect.mean
#res.l.lm=multiseq(as.matrix(dprof[ind,1:1024]),log(ccount[ind,6]),lm.approx=TRUE)
#res.l.glm=multiseq(as.matrix(dprof[ind,1:1024]),log(ccount[ind,6]),lm.approx=FALSE)
}


pdf("res_nl_lm.pdf")
plot(res.nl.lm$effect.mean,type='l')
dev.off()

pdf("res_nl_glm.pdf")
plot(res.nl.glm$effect.mean,type='l')
dev.off()