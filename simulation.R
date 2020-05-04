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
med1 <- NA
med0 <- NA

beta <- 1 - pwr0

# SIMULATION PARAMETERS
idnum <- 1234
nsample <- 10000

message("Single-Stage and Two-Stage Designs")
message("Input Parameters:")
message(paste("alpha: ", alpha))
message(paste("1-beta: ", pwr0))
message(paste("accrual rate: ", rate))
message(paste("followup period: ", b))

za = znorm(alpha) # TODO why?
za = znorm(beta)

hr <- med1/med0

hz0 <- log(2)/med0 # TODO why this?
hz1 <- log(2)/med1

message(paste("hazard rate: ", hz0, ",", hz1))
message(paste("med: ", med0, ",", med1))
message(paste("hazard ratio:", hr))

# SIMULATION
# --------------------------------------------
a1 <- 0.1
a2 <- 5
f1 <- fun(a1)
f2 <- fun(a2)

if (f1 * f2 > 0){
  for (i in 1:10){
    a1 <- a1 / 2
    a2 <- a2 * 2
    f1 <- fun(a1)
    f2 <- fun(a2)
    
    if(f1 * f2 < 0){
      break
    }
  }
}
