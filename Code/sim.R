library(dplyr)
library(plyr)
library(tidyr)
library(ggplot2)


#' Generate BinAR(1) process

generate_binar1 <- function(n, m, p, r) {
  # Initialize
  I <- numeric(n)
  
  # Initial value from stationary distribution: Bin(m, p)
  I[1] <- rbinom(1, m, p)
  
  # Generate the process
  for (t in 2:n) {
    # BinAR(1) evolution: I_t = α∘I_{t-1} + ε_t
    # where α∘ is binomial thinning with parameter r
    # and ε_t ~ Bin(m, p(1-r)/(1-p*r))
    alpha_thinned <- rbinom(1, I[t-1], r)
    epsilon <- rbinom(1, m, p*(1-r)/(1-p*r)) # Ensures that the marginal distribution is Bin(m, p)
    I[t] <- alpha_thinned + epsilon
  }
  
  return(I)
}

# Test with small example
set.seed(123)
test_counts <- generate_binar1(n = 10, m = 3, p = 0.2, r = 0.35)
print(test_counts)

#' Generate amplitude-modulating process (missing observations)

generate_O <- function(n, pi) {
  rbinom(n, 1, pi)
}

params <- list(c(m = 3, p = 0.2, r = 0.35),
               c(m = 10, p = 0.45, r = 0.5),
               c(m = 3, p = 0.2, r = 0),
               c(m = 10, p = 0.45, r = 0))
sample_sizes <- c(50, 100, 200, 500, 1000)

pi_values <- c(1, 0.95, 0.9, 0.85, 0.8, 0.75, 0.5)


