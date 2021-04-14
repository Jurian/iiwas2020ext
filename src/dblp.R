
source('src/R-functions.R')

load.dblp.ground.truth <- function() {
  truth <- fread('data/dblp-truth.tsv.gz', header = F, sep = '\t', blank.lines.skip = T, quote = '#')
  colnames(truth) <- c('uri','cluster')
  truth[,uri:= paste0('dblp:',substring(uri, 33))]
  truth[,uri:=as.factor(uri)]
  truth[,cluster:=as.integer(as.factor(cluster))]
  setkey(truth, uri)
  return(truth)
}


dblp.results <- perform.tests(load.dblp.ground.truth(), dataset = 'dblp_embedding', cutoff = seq(0.55, 0.99, 0.01), k = 2, n_trees = 1000, calc.fsc = T, normalize = T)
fwrite(dblp.results, "results/dblp.results.k2.tsv", sep = '\t', col.names = T, row.names = F)

