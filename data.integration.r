
source('lib.r')
raw = readRDS('raw.filter.rds')

runConos =function(datlp2,appname,n.cores=8){
  print(names(datlp2))
  con <- Conos$new(datlp2,n.cores=n.cores)
  con$buildGraph()
  con$findCommunities()

  con$embedGraph()

  f1=paste(appname,'_conos_clustering.pdf',sep='')
  p1=con$plotGraph()
  ggsave(f1,p1)
  con$embedGraph(method = 'UMAP',spread=7)
  return(con)
}


p2lis2=lapply(raw2,function(x) basicP2proc(x,n.cores = 10,min.transcripts.per.cell =0))

#saveRDS(p2lis2,'p2lis2.rds')

c1=runConos(p2lis2,appname)
saveRDS(c1,'dat.conos.rds')


datraw=raw2
genelists <- lapply(datraw, function(x) rownames(x))
str(genelists)
commongenes <- Reduce(intersect,genelists)
length(commongenes)
matrices_raw <- mapply(function(mm, name) {
         mm[commongenes,]
        },
    datraw,
    names(datraw))


bigM2 <- Reduce(cbind, matrices_raw)

p2=basicP2proc(bigM2,min.cells.per.gene = 0,n.cores = 10)



