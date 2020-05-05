# --------------------------------------------------------
# Minimax & Optimal 2-stage designs for 1-sampLe Log-rank test
# minimizing both expected N and study period
# Same as size.for, but replacing lambda1 with
# lambda=(lambda0+lambda1)/2 in calculating sigma_1^2
# --------------------------------------------------------
# Adapted from Kwak and Jung, 2014 
# Translated from Fortran to R

source("single_stage_functions.R")
source("two_stage_functions.R")

# --------------------------------------------------------
# INPUT PARAMETERS
# rate = expected accrual rate
# b = followup period
# alpha = probability of type 1 error
# power0 = Power
# med0, med1 = baseline and treatment group hazard?


# DEFAULTS
# These values are the defaults based on the comments 
alpha <- 0.1
power0 <- 0.9
rate <- 30
b <- 1
med0 <- 0.1
med1 <- 0.2

beta <- 1 - power0

# SIMULATION PARAMETERS
# --------------------------------------------
idnum <- 1234
nsample <- 10000

message("Single-Stage and Two-Stage Designs")
message("Input Parameters:")
message(paste("alpha: ", alpha))
message(paste("1-beta: ", power0))
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
a_limits <- select_a_limits(a_min = 0.1, a_max = 5)
a <- find_a_star(a_limits)
n <- rate * a + 1

message("OUTPUT DESIGN PARAMETERS:")
message("n:", n)  # in original code, n here is rounded down to the nearest int
message("a:", a)

# 2. Cacluate power, searching over values for a, tau, and c1 (stopping time)
# -------------------------------------------- 

# Some reservations here: 
# the values of a, tau, and c1 they propose to search over in the paper are not
# the space of parameters they search over in the code.  Here, I implement what
# the code says.

a0 <- a
a_values <- seq(0.95*a0, 1.5*a0, by = 1/rate)
# according to paper, should be 0.8*a0 to 1.5*a0
c1_values <- seq(-0.5, 1.5, by = 0.005)
# according to paper, should be -0.2 to 1.

result <- data.frame(n = NA, a = a_values, n1 = NA, tau_min = NA,
                     c1_min = NA, c_min = NA, exp_accrual_min = NA, en = NA,
                     p_early_termination_min = NA, power_min = NA,
                     d1 = NA, dd = NA)

any_solution <- FALSE

for (i in 1:length(a_values)){
  a <- a_values[i]
  
  exp_accrual_min <- a + b  # seems wrong to define this here?
  
  tau_values <- seq(1/rate, a + b, by = 1/rate)
  # according to paper, should be 0.2*a0 to 1.2*a0
  # no idea why they do 0 to a + b here.
  
  for (tau in tau_values){ 
    
    for (c1 in c1_values){
      
      params <- calculate_power(hz0, hz1, tau, a, b, c1, alpha, rate)
      
      power <- params$power
      c_star <- params$c_star
      
      if (power >= 1 - beta){
        any_solution <- TRUE
        p_early_termination <- 1 - pnorm(c1)
        exp_accrual <- a - p_early_termination * max(0, a - tau)
        
        # if we have discovered a new maximum, save the values
        if (exp_accrual < exp_accrual_min) {
          tau_min <- tau
          c1_min <- c1
          c_min <- c_star
          exp_accrual_min <- exp_accrual
          p_early_termination_min <- p_early_termination
          power_min <- power
      }
      } 
    }
  }
  if (any_solution == FALSE) { 
    warning("No solution found.")
    message(paste("a: ", a))
    row_i <- c(NA, a, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
  } else{
    n <- a * rate + 1
    n1 <- min(tau_min * rate + 1, n)
    
    en <- exp_accrual_min * rate
    d1 <- rate * tau_min * (1 - (1 - exp(-hz1 * tau_min)) / (hz1 * tau_min))
    dd <- n * (1 - exp(-hz1 * b) * (1 - exp(-hz1 * a)) / (hz1 * a))
    
    row_i <- c(n, a, n1, tau_min, c1_min, c_min, exp_accrual_min, en,
               p_early_termination_min, power_min, d1, dd) 
  }
  
  result[i,] <- row_i
}

