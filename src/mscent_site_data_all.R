library(R.utils)

setwd("~/projects/multiseq_analysis")


file.names=list.files("data/", pattern = "*.txt.gz")
chip.file.names=file.names[seq(5,310,5)]
map.file.names=file.names[seq(3,308,5)]


for(j in 1:length(chip.file.names)){
  factor.name = paste(strsplit(chip.file.names[j],"_")[[1]][1], strsplit(chip.file.names[j],"_")[[1]][2],sep="_")
  cinfo = read.table(file.path("data", chip.file.names[j]), header = TRUE)
  mapsites = read.table(file.path("data", map.file.names[j]))
  cinfo.learn = cinfo[mapsites[, 1] == 1 & cinfo[, 5] >= 13, -6]
  cinfo.infer = cinfo[mapsites[, 1] == 1 & cinfo[, 5] >= 10, -6]
  
  write.table(cinfo.learn, file.path("data/mscentipede", paste0(factor.name, "_data_learn.txt")), col.names = names(cinfo.learn), quote = FALSE, row.names = FALSE, sep = "\t")
  gzip(file.path("data/mscentipede", paste0(factor.name, "_data_learn.txt")))
  write.table(cinfo.infer, file.path("data/mscentipede", paste0(factor.name, "_data_infer.txt")), col.names = names(cinfo.infer), quote = FALSE, row.names = FALSE, sep = "\t")
  gzip(file.path("data/mscentipede", paste0(factor.name, "_data_infer.txt"))) 
}