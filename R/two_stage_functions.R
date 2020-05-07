# --------------------------------------------------------
# Supporting functions for two-stage design set-up
# --------------------------------------------------------
# Adapted from Kwak and Jung, 2014 
# Translated from Fortran to R

#' Calculate difference between type 1 error rate and alpha
#'
#' Calculates right hand side of equation 9 from Kwak and Jung, 2014.  Called by
#' find_c_star to search for appropriate critical value (c).
#'
#' @param c critical value
#' @param alpha type 1 error rate
#' @param c1 early stopping value
#' @param rho0 the ratio of sigma values: sigma_1/sigma in the paper, s11/s1 in
#'   code
#'
#' @return alpha - P(type 1 error)
type1 <- function(c, alpha, c1, rho0){
  value <- integrate(f = fun1, lower = -10, upper = c,
                     c1 = c1, rho0 = rho0)$value
  return(value - alpha)
}

#' Equation 9
#'
#' Integrand to calculate type 1 error rate from a given rho and c1
#'
#' @param z quantitiy to integrate over (from -infty to c)
#' @param c1 early stopping value
#' @param rho0 the ratio of sigma values: sigma_1/sigma in the paper, s11/s1 in
#'   code
#'
#' @return float evaluation of integrand
fun1 <- function(z, c1, rho0){
  result <- dnorm(z) * pnorm((c1 - rho0 * z)/sqrt(1 - rho0 ** 2))
  return(result)
}


#' Function 2
#' 
#' Integrand to calculate power from a given c1_bar and rho1 (text, pg 10)
#'
#' @param z quantitiy to integrate over (from -infty to c_bar)
#' @param c1_bar variable calculated from sigmas, w1 and n1 (see text, pg 10)
#' @param rho1 variable calculated from sigmas (see text, pg 9)
#'
#' @return float evaluation of integrand
fun2 <- function(z, c1_bar, rho1){
  phi <- 2 * asin(1)
  result <- 1 / (2 * phi) ** .5 * exp(-z ** 2 / 2)
  result <- result * pnorm((c1_bar - rho1 * z)/ sqrt(1 - rho1 ** 2))
  return(result)
}

#' Find c*
#'
#' Identify c*, the critical value giving the desired type 1 error rate, alpha.
#' Does so by solving equation 9 in the paper (fun1) with a binary search. Stops
#' when the range of possible values for c* is less than 0.001 and difference
#' between alpha and type 1 error rate is < 0.001. Fails if this criterion is
#' not reached within 100 iterations.
#'
#' @param alpha float - type 1 error rate
#' @param c1 float - early stopping value
#' @param rho0 the ratio of sigma values: sigma_1/sigma in the paper, s11/s1 in
#'   code
#'
#' @return critical value giving type 1 error rate alpha
find_c_star <- function(alpha, c1, rho0){
  c_min <- -3
  f_min <- type1(c_min, alpha, c1, rho0) 
  c_max <- 0
  f_max <- type1(c_max, alpha, c1, rho0)
  
  if(f_min * f_max > 0) {stop("Error finding c*: not in starting range")}
  
  for (i in 1:100){
    c_test <- (c_min + c_max)/2
    f_test <- type1(c_test, alpha, c1, rho0)
    if(f_min * f_test < 0){
      c_max <- c_test
      f_max <- f_test 
    } else{
      c_min <- c_test
      f_min <- f_test
    }
    if (abs(c_min-c_max) < 0.001 & abs(f_test < 0.001)){
      return(c_test)
    }
  }
  stop("Failed to converge on a value for c* after 100 iterations")
}


#' Calculate power (and critical value)
#'
#' @param hz0 - hazard in control group
#' @param hz1 - hazard in treatment group
#' @param tau - when to conduct interrim analysis
#' @param a - accrual period
#' @param b - followup period
#' @param c1 - stopping value
#' @param alpha - type 1 error rate
#' @param rate - expected accrual rate
#'
#' @return vector with power and critical value
calculate_power <- function(hz0, hz1, tau, a, b, c1, alpha, rate){
  # calculate sigmas and w
  hr <- hz0/hz1
  v1 <- 1 - (1 - exp(-hz0 * tau)) / (tau * hz0)
  v <- 1 - (1 - exp(-hz0 * a)) * exp(-hz0*b) / (a * hz0)
  
  s11 <- 1 - (1 - exp(-hz1 * tau)) / (tau * hz1)
  s1 <- 1 - (1 - exp(-hz1 * a)) * exp(-hz1 * b) / (a * hz1)
  s01 <- hr * s11 
  s0 <- hr * s1
  
  w1 <- s11 - s01
  w <- s1 - s0
  
  # calculate rho0, rho1
  hz.mean <- (hz0 + hz1) / 2
  s11 <- 1 - (1 - exp(-hz.mean * tau)) / (tau * hz.mean) # awful variable naming...
  s1 <- 1 - (1 - exp(-hz.mean * a)) * exp(-hz.mean*b) / (a * hz.mean) 
  rho0 <- sqrt(v1/v)
  rho1 <- sqrt(s11/s1)
  
  # find critical value
  c <- find_c_star(alpha, c1, rho0)
  
  c1_bar <- sqrt(s01 / s11) * (c1 - w1 * sqrt(rate * tau)/ sqrt(s01))
  c_bar <- sqrt(s0 / s1) * (c - w * sqrt(rate * a) / sqrt(s0))
  
  # calculate power
  power <- integrate(f = fun2, lower = -10, upper = c_bar, 
                   c1_bar = c1_bar, rho1 = rho1)$value
  
  return(c(power, c))
}