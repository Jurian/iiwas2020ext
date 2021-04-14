
install.dependencies <- function() {
  list.of.packages.cran <- c('data.table', 'R.utils', 'lpSolve', 'Matrix', 'pbapply', 'Rcpp', 'RcppGSL', 'RcppZiggurat', 'RcppArmadillo', 'RcppAnnoy', 'tictoc', 'aricode')
  new.packages.cran <- list.of.packages.cran[!(list.of.packages.cran %in% installed.packages()[,'Package'])]
  if(length(new.packages.cran)) install.packages(new.packages.cran)
}

embedding.load.glove <- function(dataset) {
  output <- list()
  
  data <- fread(file = paste0('data/', dataset, '.tsv.gz'), header = F, sep = '\t')
  data[,V1:=as.factor(V1)]
  output$dict <- data[,V1]
  output$vect <- t(as.matrix(data[,2:ncol(data)]))
  
  return(output)
}

get.nn.indexes <- function(tree, k = 1, eps = 0) {
  nn <- tree$querySelf(k = k+1, eps = eps)
  # The closest point to a query point is the point itself, remove it
  nn$nn.idx <- nn$nn.idx[,2:(k+1), drop = F]
  nn$nn.dists <- nn$nn.dists[,2:(k+1), drop = F]
  return(nn)
}

performance <- function(labels1, labels2, calc.fsc = T, calc.ami = F, calc.ari = F, calc.nmi = F) {
  
  if(!any(c(calc.fsc, calc.ami, calc.ari, calc.nmi))) {
    stop('No performance measure set')
  }
  
  out <- list()
  
  if(calc.fsc) {
    
    comembership <- comembershipTable(labels1, labels2)
    precision <- comembership[1,1] / (comembership[1,1] + comembership[1,2])
    recall <- comembership[1,1] / (comembership[1,1] + comembership[2,1])
    
    out$prc <- precision
    out$rec <- recall
    out$f10 <- 2 * ((precision * recall) / (precision + recall))
    out$f05 <- (1 + 0.5^2) * ((precision * recall) / ((0.5^2 * precision) + recall))
  }
  
  if(calc.ari) out$ari <- aricode::ARI(labels1, labels2)
  if(calc.ami) out$ami <- aricode::AMI(labels1, labels2)
  if(calc.nmi) out$nmi <- aricode::NMI(labels1, labels2)
  
  return(out)
  
}

editClustering <- function(compIndex, vect, theta, epsilon, max.in.size) {
  
  # A component of size 1 means an entity that is matched with nothing else
  if(length(compIndex) == 1) { return(compIndex) }
  # We cannot handle problems that are too large
  if(length(compIndex) > max.in.size) { return(voteClustering(compIndex, vect, theta, epsilon)) }
  # We do not have to fix a component of size 2 in any case
  if(length(compIndex) == 2) { 
    return(rep(compIndex[1],2)) 
  }
  
  #print(paste('cluster edit on component of size', length(compIndex)))
  
  n <- length(compIndex) # The number of vertexes
  e <- n*(n-1)/2 # The number of edges
  b <- round((1/6) * (n - 2) * (n - 1) * n) # The number of triangles
  c <- 3 * b # The number of constraints
  
  edge.triangles <- getallTriangles(n)
  
  # Create a matrix of coordinates of all coefficients we need to set to some non-zero value (first 2 columns)
  # plus a column of values
  f.con <- matrix(
    data = c(
      rep(1:b, times = 3, each = 3) + rep((0:2) * b, each = c), 
      rep(c(edge.triangles), times = 3),
      c(rep(c(1,1,-1),b), rep(c(1,-1,1),b), rep(c(-1,1,1),b))
    ), 
    ncol = 3
  )
  
  # Set inequality/equality signs
  f.dir <- rep('<=', c)
  
  # Set right hand side coefficients
  f.rhs <- rep(1, c)
  
  # Set coefficients of the objective function
  f.obj <- adjustedWeight(compIndex, vect, theta, epsilon)
  
  # Perform the integer linear programming operation
  ILP <- lp(
    direction = 'max', 
    objective.in = f.obj, 
    dense.const = f.con, 
    const.dir = f.dir, 
    const.rhs = f.rhs, 
    all.bin = TRUE
  )

  return(extractSolution(compIndex, ILP$solution))
}

load.annoy.index <- function(annoy.file, dim) {
  a <- new(AnnoyEuclidean, dim)
  a$load(annoy.file)
  return(a)
}

build.annoy.index <- function(annoy.file, embedding, n_trees) {
  total.nr.entities <- length(embedding$dict)
  dim <- nrow(embedding$vect)
  
  a <- new(AnnoyEuclidean, dim)
  
  for (i in 1:total.nr.entities) a$addItem(i - 1, embedding$vect[,i])
  
  a$build(n_trees)
  a$save(annoy.file)
  
  return(a)
}


perform.tests <- function(
  ground.truth, 
  dataset, 
  cutoff, 
  k, 
  epsilon = 1e-6, 
  max.size.cluster.edit = 50, 
  n_trees = 250,
  normalize = F,
  calc.fsc = T,
  calc.ami = F,
  calc.ari = F,
  calc.nmi = F) {
  
  install.dependencies()
  
  require(data.table)
  require(R.utils)
  require(lpSolve)
  require(Matrix)
  require(pbapply)
  require(Rcpp)
  require(RcppZiggurat)
  require(RcppArmadillo)
  require(RcppAnnoy)
  require(parallel)
  require(tictoc)
  if(calc.ami | calc.ari | calc.nmi) require(aricode)
  
  pbo = pboptions(type="txt")
  Rcpp::sourceCpp('src/CPP-functions.cpp')
  
  
  tryCatch ({
    
    embedding <- embedding.load.glove(dataset)
    total.nr.truth.entities <- nrow(ground.truth)
    total.nr.entities <- length(embedding$dict)
    dim <- nrow(embedding$vect)
    
    if(normalize) {
      embedding$vect <- normalize(embedding$vect)
    }
    
    annoy.file <- paste0('data/', dataset, '.annoy')
    
    if(file.exists(annoy.file)) {
      print('Loading precomputed annoy index')
      annoy <- load.annoy.index(annoy.file, dim)
    } else {
      print('Computing new annoy index')
      annoy <- build.annoy.index(annoy.file, embedding, n_trees)
    }
    
    
    tic('Computing nearest neighbors')
    # NB key2.index is zero indexed
    key2.index <- c(t(sapply(0:(total.nr.entities - 1), function(i){annoy$getNNsByItem(i, k+1)}))[,2:(k+1), drop = F])
    toc()
    annoy$unload()
    
    # NB key1.index is zero indexed
    key1.index <- rep(0:(total.nr.entities - 1), k)
    
    # Create the candidate pairs and their features
    candidate.pairs <- data.table(
      cosine.similarity = cosineSimilarityV(embedding$vect, key1.index, key2.index),
      key1.index = key1.index,
      key2.index = key2.index
    )
    
    # Remove duplicated rows
    candidate.pairs <- candidate.pairs[!duplicated(candidate.pairs)]
    
    results <- lapply(cutoff, function(theta){
      
      print(paste0('------- Processing theta = ', theta, ' -------'))
      
      valid.candidate.pairs <- candidate.pairs[cosine.similarity >= theta]
      
      tic('Finding connected components')
      # NB components is zero indexed
      components <- connectedComponents(total.nr.entities, valid.candidate.pairs$key1.index, valid.candidate.pairs$key2.index)
      toc()
      
      results <- data.table(uri = embedding$dict[unlist(components)+1], component.size = unlist(lapply(components, function(c){rep(length(c),length(c))})))
      
      print('Computing clusterings')

      output <- list()
      output$cutoff <- theta
      output$n.component <- length(components)
      output$n.too.large.comp <- sum(sapply(components, function(c){length(c) > max.size.cluster.edit}))
      output$max.size.comp <- max(sapply(components, length))
      
      
      tic('transitive closure')
      clustering <- pblapply(components, function(component){rep(component[1],length(component))}, cl = detectCores() - 1)
      output$clst.closure.scr <- clusteringCost(components, clustering, embedding$vect, theta, epsilon)
      results[,closure.result           := unlist(clustering)]
      toc()
      
      tic('vote clustering')
      clustering <- pblapply(components, voteClustering, embedding$vect, theta, epsilon, cl = detectCores() - 1)
      output$clst.vote.scr <- clusteringCost(components, clustering, embedding$vect, theta, epsilon)
      results[,corr.clst.vote.result    := unlist(clustering)]
      toc()
      
      tic('center clustering')
      clustering <- pblapply(components, centerClustering, embedding$vect, theta, epsilon, cl = detectCores() - 1)
      output$clst.center.scr <- clusteringCost(components, clustering, embedding$vect, theta, epsilon)
      results[,center.clst.result       := unlist(clustering)]
      toc()
      
      tic('center merge clustering')
      clustering <- pblapply(components, centerMergeClustering, embedding$vect, theta, epsilon, cl = detectCores() - 1)
      output$clst.center.merge.scr <- clusteringCost(components, clustering, embedding$vect, theta, epsilon)
      results[,center.merge.clst.result := unlist(clustering)]
      toc()
      
      tic('cluster edit + vote')
      clustering <- pblapply(components, editClustering, embedding$vect, theta, epsilon, max.in.size = max.size.cluster.edit, cl = detectCores() - 1)
      output$clst.edit.vote.scr <- clusteringCost(components, clustering, embedding$vect, theta, epsilon)
      results[,cluster.edit.vote.result := unlist(clustering)]
      toc()
      
      results <- merge(results, ground.truth, by = 'uri')

      print("Computing performances...")
      tic("Vote")
      output <- c(output, clst.vote = performance(results$corr.clst.vote.result, results$cluster, calc.fsc, calc.ami,  calc.ari, calc.nmi))
      toc()
      tic("Cluster edit + Vote")
      output <- c(output, clst.edit.vote = performance(results$cluster.edit.vote.result, results$cluster,  calc.fsc, calc.ami,  calc.ari, calc.nmi))
      toc()
      tic("Center cluster")
      output <- c(output, clst.center = performance(results$center.clst.result, results$cluster, calc.fsc, calc.ami,  calc.ari, calc.nmi))
      toc()
      tic("Center merge cluster")
      output <- c(output, clst.center.merge = performance(results$center.merge.clst.result, results$cluster, calc.fsc, calc.ami,  calc.ari, calc.nmi))
      toc()
      tic("Transitive closure")
      output <- c(output, clst.closure = performance(results$closure.result, results$cluster, calc.fsc, calc.ami,  calc.ari, calc.nmi))
      toc()
      
      return(output)
    })
    
    return(data.table(t(sapply(results, unlist))))
  },
  finally = {
    annoy$unload()
    #tree$delete_tree()
  })
}

