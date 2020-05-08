# --------------------------------------------------------
# Minimax & Optimal 2-stage designs for 1-sampLe Log-rank test
# minimizing both expected N and study period
# --------------------------------------------------------
# Adapted from Kwak and Jung, 2014 
# Translated from Fortran to R
# REFACTOR: Uses grid search implementation rather than `for` loops for speed
# much faster on more computationally-intensive set-ups

library(tidyverse)

source("single_stage_functions.R")
source("two_stage_functions.R")
source("refactor_functions.R")

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

# --------------------------------------------------------
# SIMULATION
# --------------------------------------------------------

# 1. Select a* for single stage design requiring the smallest sample size
# --------------------------------------------
a_limits <- select_a_limits(a_min = 0.1, a_max = 5)
a0 <- find_a_star(a_limits)
n <- rate * a0 + 1

message("SINGLE-STAGE DESIGN PARAMETERS:")
message("n:", floor(n))
message("a*:", a0)

# 2. Cacluate power, searching over values for a, tau, and c1
# -------------------------------------------- 

# Some reservations here:

# the values of a, tau, and c1 they propose to search over in the paper are not
# the space of parameters they search over in the supplemental code.  I include
# both versions here (search_power() function in refactor_functions.R sets up
# the tau values to search)

#a_values <- seq(0.95*a0, 1.5*a0, by = 1/rate) # what the fortran code does
a_values <- seq(0.8*a0, 1.5*a0, by = 1/rate) # what the paper says

#c1_values <- seq(-0.5, 1.5, by = 0.005) # what the fortran code does
c1_values <- seq(-0.2, 1, by = 0.005) # what the paper says


# -------------------------------------------- 

# set up grid of a and c1 values to search over
# (search_power() function sets up tau grid)
search_grid <- expand.grid(a_values, c1_values) %>%
  as_tibble() %>%
  rename(a = Var1, c1 = Var2)

# calculate power for all values of a, c1, and tau
result_raw <- search_grid %>% 
  group_by(a, c1) %>%
  do(search_power(.)) %>%
  ungroup()

# extract results with sufficient power
# for each a value, select the tau and c1 values which give minimal EA
result <- result_raw %>% 
  filter(power >= 1 - beta) %>%
  mutate(n = a * rate + 1,
         n1 = ifelse(tau * rate + 1 < n, tau * rate + 1, n),
         PET = 1 - pnorm(c1),
         EA = ifelse(a - tau <= 0, a, a - PET * (a - tau)),
         EN = EA * rate,
         D1 = rate * tau * (1 - (1 - exp(-hz1 * tau)) / (hz1 * tau)),
         D = floor(n) * (1 - exp(-hz1 * b) * (1 - exp(-hz1 * a)) / (hz1 * a))) %>%
  group_by(a) %>%
  filter(EA == min(EA)) %>%
  select(n, a, n1, tau, c1, c, EA, EN, PET, power, D1, D)

# runtimes of 10-20 minutes are common, depending on starting parameters
runtime <- proc.time() - t_start
