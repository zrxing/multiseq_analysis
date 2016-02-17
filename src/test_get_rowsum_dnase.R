inputFile <- "data/Ctcf_S2_dnase_data.txt.gz"
con  <- file(inputFile, open = "r")

rs.10 = NULL
rs.11 = NULL
rs.12 = NULL
rs.13 = NULL
rs.14 = NULL

base.ind.for = list()
base.ind.rev = list()

base.ind.for[[1]] = 7681:8704
base.ind.for[[2]] = 7169:9216
base.ind.for[[3]] = 6145:10240
base.ind.for[[4]] = 4097:12288
base.ind.for[[5]] = 1:16384

base.ind.rev[[1]] = 7681:8704 + 16384
base.ind.rev[[2]] = 7169:9216 + 16384
base.ind.rev[[3]] = 6145:10240 + 16384
base.ind.rev[[4]] = 4097:12288 + 16384
base.ind.rev[[5]] = 1:16384 + 16384



while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0){
  myVector = (strsplit(oneLine, " "))
  myVector = as.numeric(myVector[[1]])
  rsum = sum(myVector[c(base.ind.for[[1]], base.ind.rev[[1]])])
  rs.10 <- c(rs.10, rsum)
  rsum = sum(myVector[c(base.ind.for[[2]], base.ind.rev[[2]])])
  rs.11 <- c(rs.11, rsum)
  rsum = sum(myVector[c(base.ind.for[[3]], base.ind.rev[[3]])])
  rs.12 <- c(rs.12, rsum)
  rsum = sum(myVector[c(base.ind.for[[4]], base.ind.rev[[4]])])
  rs.13 <- c(rs.13, rsum)
  rsum = sum(myVector[c(base.ind.for[[5]], base.ind.rev[[5]])])
  rs.14 <- c(rs.14, rsum)
} 

close(con)

rsum = cbind(rs.10, rs.11, rs.12, rs.13, rs.14)

write(t(rsum), file = "test.txt", ncolumns = 5)
