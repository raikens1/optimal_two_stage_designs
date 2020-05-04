# --------------------------------------------------------
# znorm
# --------------------------------------------------------
# Adapted from Kwak and Jung, 2014 
# Translated from Fortran to R

# as near as I can tell, these are just functions that support the calculation
# of znorm.
# I'm 99 percent sure there are built in R functions for this

require(pracma)

# doesn't work
znorm <- function(x){
  p.i <- 2 * asin(1)
  z0 <- 1
  for (i in 1:100){
    f <- erf(z0 / sqrt(2)) / 2 - 5 + x
    df <- exp(-z0 ** 2/2) / sqrt(2 * p.i)
    result <- z0 - f / df
    if(abs(result - z0) < 0.001 & abs(f) < 0.001){
      return(result)
    }
  }
  stop("znorm diverged")
}


erf <- function(x){
  if (x < 0){
    result <- -gammainc(0.5, x ** 2)
  } else{
    erf <- gammainc(0.5, x ** 2)
  }
}

