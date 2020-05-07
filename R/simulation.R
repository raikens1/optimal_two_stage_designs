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
# power0 = Required power
# med0, med1 = log2/hazard for under H0 and H1

t_start = proc.time()

# SIMULATION PARAMETERS
# --------------------------------------------------------
alpha <- 0.05
power0 <- 0.9
rate <- 30
b <- 1
med0 <- log(2)/0.7
med1 <- log(2)/0.5
# --------------------------------------------------------

# PRINT PARAMETERS
# --------------------------------------------------------

beta <- 1 - power0
message("Single-Stage and Two-Stage Designs")
message("Input Parameters:")
message(paste("alpha: ", alpha))
message(paste("1-beta: ", power0))
message(paste("accrual rate: ", rate))
message(paste("followup period: ", b))

za <- -qnorm(alpha)
zb <- -qnorm(beta)

hr <- med1/med0

hz0 <- log(2)/med0
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

message("SINGLE-STAGE DESIGN PARAMETERS:")
message("n:", floor(n)) 
message("a*:", a)

# 2. Cacluate power, searching over values for a, tau, and c1
# -------------------------------------------- 

# Some reservations here: 

# the values of a, tau, and c1 they propose to search over in the paper are not
# the space of parameters they search over in the code.  Here, I implement what
# the code says.

a0 <- a
a_values <- seq(0.95*a0, 1.5*a0, by = 1/rate) # what the fortran code does
#a_values <- seq(0.8*a0, 1.5*a0, by = 1/rate) # what the paper says

c1_values <- seq(-0.5, 1.5, by = 0.005) # what the fortran code does
#c1_values <- seq(-0.2, 1, by = 0.005) # what the paper says

result <- data.frame(n = NA, a = a_values, n1 = NA, tau_min = NA,
                     c1_min = NA, c_min = NA, EA_min = NA, en = NA,
                     PET_min = NA, power_min = NA,
                     d1 = NA, dd = NA)

any_solution <- FALSE

for (i in 1:length(a_values)){
  a <- a_values[i]
  
  EA_min <- a + b  # seems wrong to define this here?
  
  tau_values <- seq(1/rate, a + b, by = 1/rate) # what the fortran code does
  #tau_values <- seq(0.2 * a0, 1.2*a0, by = 1/rate) # what the paper says
  
  for (tau in tau_values){ 
    
    for (c1 in c1_values){
      
      params <- calculate_power(hz0, hz1, tau, a, b, c1, alpha, rate)
      
      power <- params[1]
      c_star <- params[2]
      
      # check if power is acceptable
      if (power >= 1 - beta){
        any_solution <- TRUE
        PET <- 1 - pnorm(c1)
        EA <- a - PET * max(0, a - tau)
        
        # if we have discovered a new maximum, save the values
        if (EA < EA_min) {
          tau_min <- tau
          c1_min <- c1
          c_min <- c_star
          EA_min <- EA
          PET_min <- PET
          power_min <- power
      }
      } 
    }
  }
  if (any_solution == FALSE) { 
    warning(paste("No solution found for a = ", a)) 
    # common for small a to give insufficient power
    row_i <- c(NA, a, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
  } else{
    # if a solution was found, calculate and save results
    n <- a * rate + 1
    n1 <- min(tau_min * rate + 1, n)
    
    en <- EA_min * rate
    d1 <- rate * tau_min * (1 - (1 - exp(-hz1 * tau_min)) / (hz1 * tau_min))
    dd <- floor(n) * (1 - exp(-hz1 * b) * (1 - exp(-hz1 * a)) / hz1 / a)
    
    row_i <- c(n, a, n1, tau_min, c1_min, c_min, EA_min, en,
               PET_min, power_min, d1, dd) 
  }
  
  result[i,] <- row_i
}

runtime <- proc.time() - t_start
