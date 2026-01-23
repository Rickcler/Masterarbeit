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

marginal_probs <- function(m, data, n) {
  marg_probs <- numeric(m-1)
  for (i in 0:m-2) {
    marg_probs[i+1] <- sum(data <= i)*(1/n)
  } 
  return(marg_probs)
}

biv_probs <- function(m, data, n, h = 1) {
  biv_probs <- numeric(m-1)
  for (i in 0:m-2) {
    biv_probs[i+1] <- sum((data[1:(n-h)] <= i)*(data[(h+1):n] <= i))*(1/(n-h))
  } 
  return(biv_probs)
}


#' Simulation function
simulation <- function(n, m, p, r, pi, n_reps = 100) {
  results <- matrix(NA, nrow = n_reps, ncol = 4, dimnames = list(NULL, c("IOV", "Skew", "lag1_Cohen", "lag2_Cohen")))
  for (rep in 1:n_reps) {
  # Generate BinAR(1) process
  count_process <- generate_binar1(n, m, p, r)
  Missing_process <- generate_O(n, pi)
  observed_counts <- count_process[Missing_process ==1]
  CDF <- marginal_probs(m, observed_counts, length(observed_counts))
  results[rep, 1] <- (4/m)*(sum(CDF *(1-CDF)))
  results[rep, 2] <- (2/m)*sum(CDF-1)
  results[rep, 3] <- sum(biv_probs(m, observed_counts, length(observed_counts), h = 1) - CDF^2)/sum(CDF*(1-CDF))
  results[rep, 4] <- sum(biv_probs(m, observed_counts, length(observed_counts), h = 2) - CDF^2)/sum(CDF*(1-CDF))
  }
  return(rbind(colMeans(results, na.rm= T), apply(results, 2, sd)))
}



simulation(n = 100, m = 3, p = 0.2, r = 0.35, pi = 0.8, n_reps = 10)

# Alle Szenarien generieren
scenarios <- expand.grid(
  n = c(50, 100, 250, 500, 1000),
  m = c(3, 10),
  p = c(0.20, 0.45),
  r = c(0, 0.35, 0.50),
  pi = c(1, 0.95, 0.9, 0.85, 0.8, 0.75, 0.5)
)

# Filter für Paper-Szenarien
scenarios <- scenarios[
  (scenarios$m == 3 & scenarios$p == 0.20 & scenarios$r %in% c(0, 0.35)) |
  (scenarios$m == 10 & scenarios$p == 0.45 & scenarios$r %in% c(0, 0.50)),
]
results <- apply(scenarios, 1, function(row) {
  simulation(
    n = row["n"], m = row["m"], p = row["p"], 
    r = row["r"], pi = row["pi"], n_reps = 100
  )
})

target_scenarios <- scenarios[
  scenarios$m == 3 & 
  scenarios$p == 0.20 & 
  scenarios$r == 0.35 &
  scenarios$pi %in% c(1, 0.75) &
  scenarios$n %in% c(50, 100, 250, 500, 1000),
]
