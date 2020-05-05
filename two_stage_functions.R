# --------------------------------------------------------
# Supporting functions for two-stage design set-up
# --------------------------------------------------------
# Adapted from Kwak and Jung, 2014 
# Translated from Fortran to R

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
#' @return
#' @export
#'
#' @examples
fun1 <- function(z, c1, rho0){
  result <- dnorm(z) * pnorm((c1 - rho0 * z)/sqrt(1 - rho0 ** 2))
  return(result)
}


#' Function 2
#' 
#' Integrand to calculate power from a given cb1 and rho1
#'
#' @param z quantitiy to integrate over (from -infty to cb)
#' @param cb1 
#' @param rho1 
#'
#' @return
#' @export
#'
#' @examples
fun2 <- function(z, cb1, rho1){
  phi <- 2 * asin(1)
  result <- 1 / (2 * phi) ** .5 * exp(-z ** 2 / 2)
  result <- result * pnorm((cb1 - rho1 * z)/ sqrt(1 - rho1 ** 2))
  return(result)
}

#' Find c*
#'
#' Identify c*, the optimal stopping time to give the desired type 1 error rate,
#' alpha. Does so by solving equation 9 in the paper (fun1) with a binary
#' search. Stops when the range of possible values for c* is less than 0.001 and
#' difference between alpha and type 1 error rate is < 0.001. Fails if this
#' criterion is not reached within 100 iterations.
#'
#' @param alpha float - type 1 error rate
#' @param c1 float - early stopping value
#' @param rho0 the ratio of sigma values: sigma_1/sigma in the paper, s11/s1 in
#'   code
#'
#' @return
#' @export
#'
#' @examples
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
#' @param a - accrual rate
#' @param b - followup period
#' @param c1 - stopping value
#' @param alpha - type 1 error rate
#' @param rate - 
#'
#' @return
#' @export
#'
#' @examples
calculate_power <- function(hz0, hz1, tau, a, b, c1, alpha, rate){
  message("Power calculation with parameters:")
  message(paste("tau: ", tau, "\ta: ", a, "\tc1: ", c1))
  hr <- hz0/hz1
  v1 <- 1 - (1 - exp(-hz0 * tau)) / (tau * hz0)
  v <- 1 - (1 - exp(-hz0 * a)) * exp(-hz0*b) / (a * hz0)
  
  s11 <- 1 - (1 - exp(-hz1 * tau)) / (tau * hz1)
  s1 <- 1 - (1 - exp(-hz1 * a)) * exp(-hz1 * b) / (a * hz1)
  s01 <- hr * s11 
  s0 <- hr * s1
  
  w1 <- s11 - s01
  w <- s1 - s0
  
  hz.mean <- (hz0 + hz1) / 2
  s11 <- 1 - (1 - exp(-hz.mean * tau)) / (tau * hz.mean)
  s1 <- 1 - (1 - exp(-hz.mean * a)) * exp(-hz.mean*b) / (a * hz.mean)
  rho0 <- sqrt(v1/v)
  rho1 <- sqrt(s11/s1)
  
  c <- find_c_star(alpha, c1, rho0)
  
  cb1 <- sqrt(s01 / s11) * (c1 - w1 * sqrt(rate * tau)/ sqrt(s01))
  cb <- sqrt(s0 / s1) * (c - w * sqrt(rate * a) / sqrt(s0))
  
  power <- integrate(f = fun2, lower = -10, upper = cb, 
                   cb1 = cb1, rho1 = rho1)$value
  
  return(list(power = power, c_star = c))
}