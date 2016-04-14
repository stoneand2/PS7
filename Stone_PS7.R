# Andy Stone
# Problem Set 7
# April 7, 2016

# Setting working directory
setwd("~/github/PS7")

# Libraries to be utilized
library(cubature); library(mvtnorm); library(SparseGrid); library(doParallel); library(testthat)
library(plyr); library(microbenchmark)

##################################
####### The actual function ######
##################################

# Allowing for parallel processing. doParallel works for Linux, Mac OSX, Windows users (unlike doMC)
# Uses max number of cores on the user's computer (detected with detectCores())
registerDoParallel(cores=makeCluster(detectCores()))

# The sparse grid integration function, adapted to allow for a greater number of dimensions and 
# running in parallel (defaults to not running in parallel, if the user wants to run in parallel, he 
# needs to specify the number of cores above)
sg.int<-function(g, ..., lower, upper, dimensions, parallel=F){ 
  require("SparseGrid")

  # Check to make sure lower, upper are numeric
  if (class(lower) != "numeric" | class(upper) != "numeric") stop("Lower, upper need to be numeric")
  
  # Rounds values passed into lower vector down to nearest integer
  lower <- floor(lower)
  # Rounds values passed into upper vector up to nearest integer
  upper <- ceiling(upper)
  # Check to make sure no lower value is greater than any upper value (as we integrate left to right)
  if (any(lower>upper)) stop("lower must be smaller than upper")
  # Check to make sure number of dimensions isn't larger than number of elements in lower, upper
  if (dimensions > length(lower) | dimensions > length(upper)) stop("Not enough dimensions specified in vector.")
  
  # Creating a matrix of all permutations of pairings of lower[i], upper[i]-1 by 1 for every ith dimension
  # Function that will create sequence for each 
  gridFunction <- function(i){
    seq(lower[i], upper[i]-1, by=1)
  }
  # Actually making the matrix, expand.grid gives us all permutations, lapply runs function for each i
  gridss <- as.matrix(expand.grid(llply(1:dimensions, gridFunction, .parallel=parallel)))
  
  # Creating the nodes, weights used for integration
  sp.grid <- createIntegrationGrid('KPU', dimension=dimensions, k=5)
  nodes <- gridss[1,] + sp.grid$nodes
  weights <- sp.grid$weights
  # Adding more nodes, weights for each dimension greater than 1
  for (i in 2:nrow(gridss)){
    nodes <- rbind(nodes, gridss[i,] + sp.grid$nodes)  
    weights<-c(weights,sp.grid$weights)
  }
  # Integrating (finding values of the function across values of nodes), allowing for parallel by
  # switching the function to be aaply
  gx.sp <- aaply(nodes, 1, g,..., .parallel=parallel)
  # Weighted sum (using matrix multiplication) to find the integral value
  val.sp <- gx.sp %*%weights
  # Returning the value
  val.sp
}

# Function to integrate over (just n-dimensional multivariate normal distribution)
myNorm <- function(x){
  dmvnorm(x, mean=rep(0, dimensions), sigma=diag(rep(1, dimensions)))
}

###############################################################
###### Measuring gains in speed when running in parallel ###### 
###############################################################
# Setting dimensions to be 2
dimensions <- 2

# Three dimensions, not running in parallel appears to be faster
microbenchmark(sg.int(myNorm, lower=c(1,1,2,3,4,5), upper=c(6,6,6,7,8,9), dimensions=2, parallel=T), 
               sg.int(myNorm, lower=c(1,1,2,3,4,5), upper=c(6,6,6,7,8,9), dimensions=2, parallel=F),
               times=10)

dimensions <- 4
# Four dimensions, again, not running in parallel seems to be faster
microbenchmark(sg.int(myNorm, lower=c(1,1,2,3,4,5), upper=c(6,6,6,7,8,9), dimensions=4, parallel=T), 
               sg.int(myNorm, lower=c(1,1,2,3,4,5), upper=c(6,6,6,7,8,9), dimensions=4, parallel=F),
               times=10)



#########################################################################################
# Integrating using adaptIntegrate, comparing speed/accuracy with sparse grid algorithm #
#########################################################################################
### ACCURACY ###

# Actual answer (uses the actual distribution function)
dimensions <- 2
ans <- as.numeric(pmvnorm(upper=rep(.5, dimensions), mean=rep(0, dimensions), sigma=diag(rep(1, dimensions))))
# Integrate from -100 (where density is practically zero) to 0.5
# Both functions are pretty accurate here in two dimensions, but adaptIntegrate is better
# Error is -7.7156e-08
adaptIntegrate(myNorm, lowerLimit=rep(-100, dimensions), upperLimit=rep(.5, dimensions))$integral - ans
# Error is -0.012055 (would be better if lower is -100. But that takes a very long time)
sg.int(myNorm, lower=rep(-1, dimensions), upper=rep(.5, dimensions), dimensions=2, parallel=T) - ans

### SPEED ###

# Speed comparison: adaptiveIntegrate is slower in this case
microbenchmark(adaptIntegrate(myNorm, lowerLimit=rep(-1, dimensions), upperLimit=rep(.5, dimensions)), 
               sg.int(myNorm, lower=rep(-1, dimensions), upper=rep(.5, dimensions), dimensions=2, parallel=T),
               times=100)
# With many cores specified for our function, the differene is about the same
microbenchmark(adaptIntegrate(myNorm, lowerLimit=rep(-1, dimensions), upperLimit=rep(.5, dimensions)), 
               sg.int(myNorm, lower=rep(-1, dimensions), upper=rep(.5, dimensions), dimensions=2, parallel=T),
               times=100)

# Setting number of dimensions to be 3 with this multivariate normal distribution takes a VERY LONG time

# Three dimensions, lower=c(1,1,2), upper=c(6,6,6)
# With three dimensions, the algorithms give quite different answers
dimensions <- 3
adaptIntegrate(myNorm, lowerLimit=c(1,1,2), upperLimit=c(6,6,6), fDim=3)$integral
sg.int(myNorm, lower=c(1,1,2), upper=c(6,6,6), dimensions=3, parallel.cores=1)

# With higher dimensions, adaptiveIntegrate gets a little slower in comparison to our function
microbenchmark(adaptIntegrate(myNorm, lowerLimit=c(1,1,2), upperLimit=c(6,6,6), fDim=3)$integral, 
               sg.int(myNorm, lower=c(1,1,2), upper=c(6,6,6), dimensions=3, parallel.cores=8),
               times=10)


######################
#### UNIT TESTING ####
######################

# Function to compare output of function to correct answer
# Integrating a multivariate normal distribution (so we can know the true answer)
# Tolearance can be set by user
dimensions <- 2
within_tolerance <- function(tolerance){
test_that(paste('Function within',as.character(tolerance),'of true answer'),{
          expect_equal(as.numeric(pmvnorm(upper=rep(.5, dimensions), mean=rep(0, dimensions), sigma=diag(rep(1, dimensions)))),
            as.numeric(sg.int(myNorm, lower=rep(-1, dimensions), upper=rep(.5, dimensions), dimensions=2, parallel=F)),
            tolerance=tolerance)
}
)
}
# within_tolerance(2) # Passes
# within_tolerance(.002) # Fails

# Testing to make sure input is in numeric vector form
numeric_inputs <- function(lower, upper){
test_that('Inputs of bounds are numeric',{
  expect_is(lower, "numeric")
  expect_is(upper, "numeric")
  })
}
# numeric_inputs(lower=c(-1,-1), upper=c(1,1)) # Passes
# numeric_inputs(lower="-1", upper="1") # Fails

# Testing to make sure dimensions specified is numeric
numeric_dimensions <- function(dimensions){
  test_that('Dimension is numeric',{
    expect_is(dimensions, "numeric")
  })
}
# numeric_dimensions(2) # Passes
# numeric_dimensions("2") # Fails


################################### 
##### Monte Carlo integration ##### 
################################### 

# Function for Monte Carlo integration
integrateMonteCarlo <- function(func, lower, upper, n, dimensions, parallel=F){
  # Creating the distribution of points to integrate under
  random.points <- matrix(runif(n=(n*dimensions), min=lower, max=upper), ncol=dimensions)
  # Evaluating random points at function to get the x values of our random points
  x.values <- unname(aaply(random.points, 1, func, .parallel=parallel))
  # "Integrating" by taking count of number of points under curve (mean), multiplying by length of integral
  # and taking it to the power of how many dimensions we have
  return(mean(x.values)*(upper-lower)^dimensions)
}

# Example
dimensions <- 3
integrateMonteCarlo(myNorm, -2, 2, n=100000, dimensions=3, parallel=F)



