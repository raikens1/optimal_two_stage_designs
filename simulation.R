# --------------------------------------------------------
# Minimax & Optimal 2-stage designs for 1-sampLe Log-rank test
# minimizing both expected N and study period
# Same as size.for, but replacing lambda1 with
# lambda=(lambda0+lambda1)/2 in calculating sigma_1^2
# --------------------------------------------------------
# Adapted from Kwak and Jung, 2014 
# Translated from Fortran to R

source("functions.R")
source("znorm.R")

# --------------------------------------------------------
# INPUT PARAMETERS
# rate = expected accrual rate
# b = followup period
# alpha = probability of type 1 error
# pwr0 = Power
# med0, med1 = baseline and treatment group hazard?


# DEFAULTS
# These values are the defaults based on the comments 
alpha <- 0.1
pwr0 <- 0.9
rate <- 30
b <- 1
med0 <- 0.1
med1 <- 0.2

beta <- 1 - pwr0

# SIMULATION PARAMETERS
# --------------------------------------------
idnum <- 1234
nsample <- 10000

message("Single-Stage and Two-Stage Designs")
message("Input Parameters:")
message(paste("alpha: ", alpha))
message(paste("1-beta: ", pwr0))
message(paste("accrual rate: ", rate))
message(paste("followup period: ", b))

za <- -qnorm(alpha)
zb <- -qnorm(beta)

hr <- med1/med0

hz0 <- log(2)/med0 # TODO why this?
hz1 <- log(2)/med1

message(paste("hazard rate: ", hz0, ",", hz1))
message(paste("med: ", med0, ",", med1))
message(paste("hazard ratio:", hr))

# --------------------------------------------
# SIMULATION
# --------------------------------------------

# 1. Select a* for single stage design requiring the smallest sample size
# --------------------------------------------
a_range <- select_a_range(a_min = 0.1, a_max = 5)
a <- find_a_star(a_range)
n <- rate * a + 1

message("OUTPUT DESIGN PARAMETERS:")
message("n:", n)  # in original code, n here is rounded down to the nearest int
message("a:", a)

# 2. Cacluate power, searching over values for a, tau, and c1
# --------------------------------------------

a0 <- a
ntau <- 1000

for (ia in 0:1000){
  a <- 0.95 * a0 + ia/rate
  if (a > 1.5*a0){
    break
  }
  new <- 0
  eam <- a + b
  ntau <- rate * (a + b)
  
  for (itau in 0:ntau){
    tau <- itau/rate
    
    for (ic1 in 0:1000){
      c1 <- -0.5 + 0.005 * ic1
      
      if(c1 < 1.5){
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
        
        c21=-3
        f1=type1(c21)
        c22=0
        f2=type1(c22)
        
      }
    }
  }
}
