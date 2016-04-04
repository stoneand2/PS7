# Andy Stone
# Problem Set 7
# April 7, 2016

sg.int<-function(g, ..., lower, upper){ 
  require("SparseGrid")
  
  # Rounds down to integers
  lower <- floor(lower)
  # Rounds up to integers
  upper <- ceiling(upper)
  # Checks to make sure no lower value is greater than any upper value
  if (any(lower>upper)) stop("lower must be smaller than upper")
  
  # Matrix of all permutations of pairings of lower[1], upper[1]-1 by 1, and same for 2nd index
  gridss <- as.matrix(expand.grid(seq(lower[1], upper[1]-1, by=1), seq(lower[2], upper[2]-1, by=1)))
  # Creates nodes, weights used for integration
  sp.grid <- createIntegrationGrid('KPU', dimension=2, k=5)
  nodes <- gridss[1,] + sp.grid$nodes
  weights <- sp.grid$weights
  
  for (i in 2:nrow(gridss)){
    nodes <- rbind(nodes, gridss[i,] + sp.grid$nodes)  
    weights<-c(weights,sp.grid$weights)
  }
  gx.sp <- apply(nodes, 1, g,...)
  val.sp <- gx.sp %*%weights
  val.sp
}

mixDist <- function(x){
  (.8*dnorm(x, mean=1, sd=1)+.2*dnorm(x, mean=-1, sd=.4))
}


# Function allowing greater number of dimensions
sg.int<-function(g, ..., lower, upper, dimensions, parallel.cores=1){ 
  require("SparseGrid")
  # Allowing for parallel processing
  require("doMC")
  registerDoMC(cores=parallel.cores)
  
  # Rounds down to integers
  lower <- floor(lower)
  # Rounds up to integers
  upper <- ceiling(upper)
  # Check to make sure no lower value is greater than any upper value
  if (any(lower>upper)) stop("lower must be smaller than upper")
  # Check to make sure number of dimensions isn't larger than number of elements in lower, upper
  if (dimensions > length(lower) | dimensions > length(upper)) stop("Not enough dimensions specified in vector.")
  
  # Matrix of all permutations of pairings of lower[1], upper[1]-1 by 1, and same for 2nd index
  gridFunction <- function(i){
    seq(lower[i], upper[i]-1, by=1)
  }
  gridss <- as.matrix(expand.grid(lapply(1:dimensions, gridFunction)))
  
  # gridss <- as.matrix(expand.grid(seq(lower[1], upper[1]-1, by=1), seq(lower[2], upper[2]-1, by=1)))
  # Creates nodes, weights used for integration
  sp.grid <- createIntegrationGrid('KPU', dimension=dimensions, k=5)
  nodes <- gridss[1,] + sp.grid$nodes
  weights <- sp.grid$weights
  
  for (i in 2:nrow(gridss)){
    nodes <- rbind(nodes, gridss[i,] + sp.grid$nodes)  
    weights<-c(weights,sp.grid$weights)
  }
  gx.sp <- apply(nodes, 1, g,...)
  val.sp <- gx.sp %*%weights
  val.sp
}







