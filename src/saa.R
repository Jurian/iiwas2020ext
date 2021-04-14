
source('src/R-functions.R')

load.saa.ground.truth <- function() {
  truth <- fread('data/saa-truth.tsv.gz',  header = T, sep = '\t', blank.lines.skip = T, quote = '#')
  truth <- truth[!is.na(truth$ILN_GROUP)]
  truth[,uri:=as.factor(paste0(DATASET, RESOURCE_URI))]
  truth[,cluster:=as.integer(as.factor(paste0(as.integer(as.factor(CLUSTER_ID)), ILN_GROUP)))]
  setkey(truth, uri)
  return(truth[,list(uri,cluster)])
}


saa.results <- perform.tests(load.saa.ground.truth(), dataset = 'saa_embedding', cutoff = seq(0.50, 0.99, 0.01), k = 2, n_trees = 1000, calc.fsc = T, normalize = T)
fwrite(saa.results, "results/saa.results.k2.tsv", sep = '\t', col.names = T, row.names = F)
