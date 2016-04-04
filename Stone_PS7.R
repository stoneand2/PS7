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
  # Allowing for parallel processing. doParallel works for Linux, Mac OSX, Windows users (unlike doMC) 
  require("doParallel")
  registerDoParallel(cores=parallel.cores)
  
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


# Measuring gains in speed when running in parallel
library(microbenchmark)
# Three dimensions, 1 and 4 cores
microbenchmark(sg.int(mixDist, lower=c(1,1,2,3,4,5), upper=c(6,6,6,7,8,9), dimensions=3, parallel.cores=1), 
               sg.int(mixDist, lower=c(1,1,2,3,4,5), upper=c(6,6,6,7,8,9), dimensions=3, parallel.cores=4),
               times=100)
# Four dimensions, 1 and 8 cores
microbenchmark(sg.int(mixDist, lower=c(1,1,2,3,4,5), upper=c(6,6,6,7,8,9), dimensions=4, parallel.cores=1), 
               sg.int(mixDist, lower=c(1,1,2,3,4,5), upper=c(6,6,6,7,8,9), dimensions=4, parallel.cores=8),
               times=10)

# Integrating using adaptIntegrate, comparing speed/accuracy with sparse grid algorithm
library(cubature);library(mvtnorm)

# Function to integrate over (just x-dimensional multivariate normal distribution)
myNorm <- function(x){
  dmvnorm(x, mean=rep(0, dimensions), sigma=diag(rep(1, dimensions)))
}
# Setting number of dimentions
dimensions <- 2

### SPEED ###

# Actual answer (uses the actual distribution function)
ans <- as.numeric(pmvnorm(upper=rep(.5, dimensions), mean=rep(0, dimensions), sigma=diag(rep(1, dimensions))))
# Integrate from -100 (where density is practically zero) to 0.5
# Both functions are pretty accurate here in two dimensions, but adaptIntegrate is better
# Error is -7.7156e-08
adaptIntegrate(myNorm, lowerLimit=rep(-100, dimensions), upperLimit=rep(.5, dimensions))$integral - ans
# Error is -0.012055
sg.int(myNorm, lower=rep(-1, dimensions), upper=rep(.5, dimensions), dimensions=2, parallel.cores=1) - ans

### ACCURACy ###

# Speed comparison: adaptiveIntegrate is slower in this case
microbenchmark(adaptIntegrate(myNorm, lowerLimit=rep(-100, dimensions), upperLimit=rep(.5, dimensions)), 
               sg.int(myNorm, lower=rep(-1, dimensions), upper=rep(.5, dimensions), dimensions=2, parallel.cores=1),
               times=100)
# With many cores specified for our function, the differene is about the same
microbenchmark(adaptIntegrate(myNorm, lowerLimit=rep(-100, dimensions), upperLimit=rep(.5, dimensions)), 
               sg.int(myNorm, lower=rep(-1, dimensions), upper=rep(.5, dimensions), dimensions=2, parallel.cores=8),
               times=100)


# Three dimensions, lower=c(1,1,2), upper=c(6,6,6)
# With three dimensions, the algorithms give quite different answers
adaptIntegrate(mixDist, lowerLimit=c(1,1,2), upperLimit=c(6,6,6), fDim=3)$integral
sg.int(mixDist, lower=c(1,1,2), upper=c(6,6,6), dimensions=3, parallel.cores=1)

# With higher dimensions, adaptiveIntegrate gets a little slower in comparison to our function
microbenchmark(adaptIntegrate(mixDist, lowerLimit=c(1,1,2), upperLimit=c(6,6,6), fDim=3)$integral, 
               sg.int(mixDist, lower=c(1,1,2), upper=c(6,6,6), dimensions=3, parallel.cores=8),
               times=10)









