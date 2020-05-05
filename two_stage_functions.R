# --------------------------------------------------------
# Supporting functions for two-stage design set-up
# --------------------------------------------------------
# Adapted from Kwak and Jung, 2014 
# Translated from Fortran to R

type1 <- function(c){
  value <- integrate(fun1, -10, c)$value
  return(value - alpha)
}

# annoyingly, this returns values very close to - but not exactly - what the
# original code produces - within about (+/-0.001).  Maybe it's just that R is
# using a more accurate approcimation for the normal cdf?
fun1 <- function(z){
  # depends on c1 and rho0 from global environment
  phi <- 2 * asin(1)
  result <- 1 / sqrt(2 * phi) * exp(-z ** 2 / 2)
  result <- result * pnorm((c1 - rho0 * z))/sqrt(1 - rho0 ** 2)
  return(result)
}

# much better fidelity to fun2 from original code (=< 0.0001 diffference)
fun2 <- function(z){
  # depends on cb1 and rho1 from global environment
  phi <- 2 * asin(1)
  result <- 1 / (2 * phi) ** .5 * exp(-z ** 2 / 2)
  result <- result * pnorm((cb1 - rho1 * z)/ sqrt(1 - rho1 ** 2))
  return(result)
}

#' Find c*
#'
#' Perform a binary search to identify a*, the optimal stopping time.  Stops
#' when the range of possible values for c* is less than 0.001 and difference
#' between alpha and type 1 error rate is < 0.001. Fails if this criterion is
#' not reached within 100 iterations.
#'
#' @param a_range
#'
#' @return
#' @export
#'
#' @examples
find_c_star <- function(){
  c21 <- -3
  f1 <- type1(c21) 
  c22 <- 0
  f2 <- type1(c22)
  
  if(f1 * f2 > 0) {stop("Error finding c*: not in starting range")}
  
  for (i in 1:100){
    c23 <- (c21 + c22)/2
    f3 <- type1(c23)
    if(f1 * f3 < 0){
      c22 <- c23
      f2 <- f3 
    } else{
      c21 <- c23
      f1 <- f3
    }
    if (abs(c21-c22) < 0.001 & abs(f3 < 0.001)){
      return(c23)
    }
  }
  stop("Failed to converge on a value for c* after 100 iterations")
}