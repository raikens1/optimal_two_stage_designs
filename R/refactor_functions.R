#' Search Power
#'
#' Search over values of tau to find optimal power for a given a and c1.  NOTE:
#' depends on a0, hz0, hz1, b, alpha, and rate being correctly defined in the
#' global environment.
#'
#' @param a - accrual period
#' @param c1 - early stopping value
#'
#' @return
search_power <- function(row){
  a <- row$a 
  c1 <- row$c1
  
  tau_values <- seq(1/rate, a + b, by = 1/rate) # what the fortran code does
  #tau_values <- seq(0.2 * a0, 1.2*a0, by = 1/rate) # what the paper says
  
  tau_results <- sapply(tau_values, calculate_power,
                        hz0 = hz0, hz1 = hz1, a = a, b = b,
                        c1 = c1, alpha = alpha, rate = rate) %>%
    t() %>%
    as_tibble() %>% 
    rename(power = V1, c = V2) %>%
    mutate(tau = tau_values)
  
  return(tau_results)
}