library(R.utils)

# dprof=read.table("data/Batf_S356_dnase_data.txt.gz")
cinfo=read.table("data/Batf_S356_site_and_chip_data.txt.gz",header=TRUE)
mapsites=read.table("data/Batf_S356_mappability_flags_large.txt.gz")
# dprof=dprof[mapsites[,1]==1&cinfo[,5]>=13,]
cinfo=cinfo[mapsites[,1]==1&cinfo[,5]>=13,-6]


write.table(cinfo,"data/Batf_S356_learn_data.txt",col.names=names(cinfo),quote=FALSE,row.names=FALSE,sep="\t")
gzip("data/Batf_S356_learn_data.txt")

