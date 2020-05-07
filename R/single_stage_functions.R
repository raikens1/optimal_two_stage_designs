# --------------------------------------------------------
# Supporting functions for single stage design set-up
# --------------------------------------------------------
# Adapted from Kwak and Jung, 2014 
# Translated from Fortran to R

#' Equation 3 from paper
#'
#' Gives difference between right and left hand sides of equation 3 in paper (I
#' think?).  Applied by find_a_star to solve for the value of a giving minimal
#' sample size in a one-stage study
#' 
#' NOTE: depends on b, hz1, hz0, za, and zb from global environment.
#'
#' @param a an accrual period
#'
#' @return difference between right and left hand sides of equation 3 in paper (I
#' think?)
#'
#' @examples
fun <- function(a) {
  v1 <- 1 - exp(-hz1 * b)/ (a * hz1) + exp(-hz1 * (a + b))/ (a * hz1)
  
  v0 <- v1 * hz0 / hz1 
  hr <- hz0 / hz1
  w <- v1 - v0
  
  # expected number of events in combined treat and control groups at final analysis
  hz.mean <- (hz0 + hz1) / 2 
  v <- 1 - exp(-hz.mean * b)/ (a * hz.mean) + exp(-hz.mean * (a + b))/ (a * hz.mean)

  result <- rate * a - (sqrt(v0) * za + sqrt(v) * zb) ** 2 / w ** 2
  
  return(result)
}

#' Select a_limits
#'
#' Select a range of values to search for the optimal accrual period a* for a
#' single stage design with the smallest sample size.  The range must be
#' selected so that fun(a_min) and fun(a_max) have opposite signs.  Stops after
#' 10 iterations
#'
#' @param a_min
#' @param a_max
#'
#' @return tuple, range of a values to search
select_a_limits <- function(a_min, a_max){
  f_min <- fun(a_min)
  f_max <- fun(a_max)
  
  # if f_min and f_max have the same sign,
  # move a_min and a_max farther apart until they have different signs
  if (f_min * f_max > 0) {
    for (i in 1:10){
      a_min <- a_min / 2
      a_max <- a_max * 2
      f_min <- fun(a_min)
      f_max <- fun(a_max)
      message(paste(f_min, ", ", f_max))
      
      if (f_min * f_max < 0) {
        break
      }
    }
    stop("failed to identify search range for a*")
  }
  
  return(c(a_min, a_max))
}

#' Find a*
#'
#' Given a search range, a_limits, perform a binary search to identify a*, the
#' optimal accrual period for a single-stage design requiring the smallest
#' sample size.  Stops when the range of possible values for a* is sufficiently
#' small that the optimal sample size is the same for all values in that range.
#' Fails if this criterion is not reached within 100 iterations.
#'
#' @param a_limits
#'
#' @return
#' @export
#'
#' @examples
find_a_star <- function(a_limits){
  a_min <- a_limits[1]
  a_max <- a_limits[2]
  
  f_min <- fun(a_min)
  f_max <- fun(a_max)
  
  # binary search
  # move a_min and a_max closer while keeping the signs of f_min and f_max opposite
  for (i in 1:100){
    a_test <- (a_min + a_max)/2
    f_test <- fun(a_test)
    message(paste("a_test:", a_test))
    message(paste("f_test", f_test))
    
    if(f_min * f_test < 0) { # f_test has same sign as f_max
      a_max <- a_test
    } else{ # f_test has same sign as f_min
      a_min <- a_test
    }
    
    # stop if the range a_min to a_max is small enough that sample size is the
    # same across the whole interval
    if (abs(rate * (a_min - a_max)) < 1) {
      message(abs(rate * (a_min - a_max)))
      break
    }
  }
  
  
  if (abs(rate * (a_min - a_max) > 1)) {
    stop("failed to converge to a value of a* in 100 iterations")
  }
  
  return(a_max)
}