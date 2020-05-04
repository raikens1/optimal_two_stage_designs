# --------------------------------------------------------
# Supporting functions
# --------------------------------------------------------
# Adapted from Kwak and Jung, 2014 
# Translated from Fortran to R

#' Fun
#' 
#' TODO: this needs a better name and description
#'
#' @param a the accrual period
#'
#' @return result TODO: what is this?
#'
#' @examples
fun <- function(a) {
  # expected number of events at final analysis under H1
  v1 <- 1 - (exp(-hz1 * b) + exp(-hz1 * (a + b))) / (a * hz1) 
  
  # expected number of events under H0
  v0 <- v1 * hz0 / hz1 
  hr <- hz0 / hz1 # TODO why reassign hr?  does this change the global value?
  w <- v1 - v0
  
  hz.1 <- (hz0 + hz1) / 2 # TODO what is this?
  v <- 1 - (exp(-hz.1 * b) + exp(-hz.1 * (a + b))) / (a * hz.1)

  result <- rate * a - (sqrt(v0) * za + sqrt(v) * zb) ** 2 / w ** 2
  
  return(result)
}

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
    result <- -gammp(0.5, x ** 2)
  } else{
    erf <- gammp(0.5, x ** 2)
  }
}


gammp <- function(a, x){
  if(x < 0 | a <= 0 ){
    stop("bad arguments in gammp")
  }
  if(x < a + 1){
    gser(gamser, a, x, gLn) # TODO: ???? what are gamser and gln?
    result <- gamser
  } else {
    gcf(gammcf, a, x, gLn) # TODO: ??? what are gammcf and gLn?
    result < 1 - gammcf
  }
  return(result)
}


