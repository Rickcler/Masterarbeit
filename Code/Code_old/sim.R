library(dplyr)
library(tidyr)
library(ggplot2)

# ==============================================================================
# Data Generating Processes
# ==============================================================================

#' Generate BinAR(1) process
#' @param n    Length of the time series
#' @param m    Binomial parameter
#' @param p    Mean of the stationary distribution E[I_t]
#' @param r    Thinning parameter (must be in [0, 1])
generate_binar1 <- function(n, m, p, r) {
  alpha <- (p * (1 - r)) + r
  beta  <- p * (1 - r)
  I     <- numeric(n)
  I[1]  <- rbinom(1, m, p)
  for (t in 2:n) {
    I[t] <- rbinom(1, I[t-1], alpha) + rbinom(1, m - I[t-1], beta)
  }
  return(I)
}



#' Generate amplitude-modulating missingness process
#' @param n   Length of the time series
#' @param pi  Observation probability
generate_O <- function(n, pi) {
  rbinom(n, 1, pi)
}



#' Generate missingness process dependent on X_t (MAR)
#' 
#' The observation probability increases with X_t, so that higher
#' values of the process are more likely to be observed. This induces
#' dependence between (O_t) and (X_t), violating the MCAR assumption.
#'
#' @param x    The realized time series (X_t values)
#' @param m    Binomial parameter (max value of X_t)
#' @param pi_low   Observation probability when X_t = 0
#' @param pi_high  Observation probability when X_t = m
generate_O_MAR <- function(x, m, pi_low = 0.5, pi_high = 0.9) {
  pi_t <- pi_low + (pi_high - pi_low) * (x / m)
  rbinom(length(x), 1, pi_t)
}

generate_O_MAR_inv <- function(x, m, pi_low = 0.5, pi_high = 0.9) {
  # pi_t fällt mit x: hohe Werte seltener beobachtet
  pi_t <- pi_high - (pi_high - pi_low) * (x / m)
  rbinom(length(x), 1, pi_t)
}

#' Generate i.i.d. Binomial process
#' @param p  Success probability
#' @param m  Binomial parameter
#' @param n  Length of the time series
generate_iid <- function(p, m, n) {
  rbinom(n, m, p)
}


# ==============================================================================
# Estimation Functions
# ==============================================================================

#' Empirical marginal CDF
#' @param m     Number of categories minus 1
#' @param data  Observed data vector
#' @param n     Number of observations
marginal_probs_e <- function(m, data, n) {
  marg_probs <- numeric(m)
  for (i in 0:(m - 1)) {
    marg_probs[i + 1] <- sum(data <= i) * (1 / n)
  }
  return(marg_probs)
}

#' Empirical bivariate lag-h CDF (under missingness)
#' @param m   Number of categories minus 1
#' @param data  Full (unobserved) data vector
#' @param O   Missingness indicator vector
#' @param n   Length of the time series
#' @param h   Lag
biv_probs_e <- function(m, data, O, n, h = 1) {
  biv_probs <- numeric(m)
  valid <- (O[1:(n - h)] == 1) & (O[(h + 1):n] == 1)
  for (i in 0:(m - 1)) {
    biv_probs[i + 1] <- sum(
      (data[1:(n - h)] <= i) *
      (data[(h + 1):n] <= i) *
      valid
    ) / sum(valid)
  }
  return(biv_probs)
}


# ==============================================================================
# True Parameter Functions
# ==============================================================================

#' True IOV under BinAR(1) stationary distribution
#' @param m  Binomial parameter
#' @param p  Success probability
true_IOV <- function(m, p) {
  F <- pbinom(0:(m - 1), m, p)
  return((4 / m) * sum(F * (1 - F)))
}


# ==============================================================================
# Simulation Functions
# ==============================================================================

#' Simulation: returns summary statistics and estimated CDF across replications
#' @param n      Time series length
#' @param m      Binomial parameter
#' @param p      Success probability
#' @param r      Thinning parameter
#' @param pi     Observation probability
#' @param pi_h   Thinning parameter for missingness 
#' @param n_reps Number of replications
#' @return List with two matrices (summary, cdf), each with rows mean and sd
simulation <- function(n, m, p, r,  pi, pi_h, n_reps = 1000) {
  results <- matrix(
    NA, nrow = n_reps, ncol = 4,
    dimnames = list(NULL, c("IOV", "Skew", "lag1_Cohen", "lag2_Cohen"))
  )
  cdf_results <- matrix(
    NA, nrow = n_reps, ncol = m,
    dimnames = list(NULL, paste0("f_", 0:(m - 1)))
  )

  for (rep in 1:n_reps) {
    # BinAR(1) process with missingness
    count_process   <- generate_binar1(n, m, p, r)
    Missing_process <- generate_binar1(n, 1, pi, pi_h)
    observed_counts <- count_process[Missing_process == 1]
    CDF             <- marginal_probs_e(m, observed_counts, length(observed_counts))

    # i.i.d. process for Cohen's kappa
    iid_process         <- generate_iid(p, m, n)
    observed_counts_iid <- iid_process[Missing_process == 1]
    CDF_iid             <- marginal_probs_e(m, observed_counts_iid, length(observed_counts_iid))

    cdf_results[rep, ]  <- CDF
    results[rep, 1]     <- (4 / m) * sum(CDF * (1 - CDF))
    results[rep, 2]     <- (2 / m) * sum(CDF - 1)
    results[rep, 3]     <- sum(biv_probs_e(m, iid_process, Missing_process, n, h = 1) - CDF_iid^2) /
                           sum(CDF_iid * (1 - CDF_iid))
    results[rep, 4]     <- sum(biv_probs_e(m, iid_process, Missing_process, n, h = 2) - CDF_iid^2) /
                           sum(CDF_iid * (1 - CDF_iid))
  }

  list(
    summary = rbind(colMeans(results,     na.rm = TRUE), apply(results,     2, sd)),
    cdf     = rbind(colMeans(cdf_results, na.rm = TRUE), apply(cdf_results, 2, sd))
  )
}

#' Simulation under MAR (missingness depends on X_t)
#' @param n        Time series length
#' @param m        Binomial parameter
#' @param p        Success probability
#' @param r        Thinning parameter
#' @param pi_low   Observation probability when X_t = 0
#' @param pi_high  Observation probability when X_t = m
#' @param n_reps   Number of replications
simulation_MAR <- function(n, m, p, r, pi_low = 0.5, pi_high = 0.9,
                            n_reps = 1000) {
  results <- matrix(
    NA, nrow = n_reps, ncol = 4,
    dimnames = list(NULL, c("IOV", "Skew", "lag1_Cohen", "lag2_Cohen"))
  )

  for (rep in 1:n_reps) {
    count_process   <- generate_binar1(n, m, p, r)
    Missing_process <- generate_O_MAR(count_process, m, pi_low, pi_high)
    observed_counts <- count_process[Missing_process == 1]
    CDF             <- marginal_probs_e(m, observed_counts,
                                        length(observed_counts))

    iid_process         <- generate_iid(p, m, n)
    Missing_iid         <- generate_O_MAR(iid_process, m, pi_low, pi_high)
    observed_counts_iid <- iid_process[Missing_iid == 1]
    CDF_iid             <- marginal_probs_e(m, observed_counts_iid,
                                             length(observed_counts_iid))

    results[rep, 1] <- (4 / m) * sum(CDF * (1 - CDF))
    results[rep, 2] <- (2 / m) * sum(CDF - 1)
    results[rep, 3] <- sum(
      biv_probs_e(m, iid_process, Missing_iid, n, h = 1) - CDF_iid^2
    ) / sum(CDF_iid * (1 - CDF_iid))
    results[rep, 4] <- sum(
      biv_probs_e(m, iid_process, Missing_iid, n, h = 2) - CDF_iid^2
    ) / sum(CDF_iid * (1 - CDF_iid))
  }

  list(
    summary = rbind(colMeans(results, na.rm = TRUE),
                    apply(results, 2, sd))
  )
}

simulation_MAR_inv <- function(n, m, p, r, pi_low = 0.5, pi_high = 0.9,
                                n_reps = 1000) {
  results <- matrix(
    NA, nrow = n_reps, ncol = 4,
    dimnames = list(NULL, c("IOV", "Skew", "lag1_Cohen", "lag2_Cohen"))
  )

  for (rep in 1:n_reps) {
    count_process   <- generate_binar1(n, m, p, r)
    Missing_process <- generate_O_MAR_inv(count_process, m, pi_low, pi_high)
    observed_counts <- count_process[Missing_process == 1]
    CDF             <- marginal_probs_e(m, observed_counts,
                                        length(observed_counts))

    iid_process         <- generate_iid(p, m, n)
    Missing_iid         <- generate_O_MAR_inv(iid_process, m, pi_low, pi_high)
    observed_counts_iid <- iid_process[Missing_iid == 1]
    CDF_iid             <- marginal_probs_e(m, observed_counts_iid,
                                             length(observed_counts_iid))

    results[rep, 1] <- (4 / m) * sum(CDF * (1 - CDF))
    results[rep, 2] <- (2 / m) * sum(CDF - 1)
    results[rep, 3] <- sum(
      biv_probs_e(m, iid_process, Missing_iid, n, h = 1) - CDF_iid^2
    ) / sum(CDF_iid * (1 - CDF_iid))
    results[rep, 4] <- sum(
      biv_probs_e(m, iid_process, Missing_iid, n, h = 2) - CDF_iid^2
    ) / sum(CDF_iid * (1 - CDF_iid))
  }
 
  list(
    summary = rbind(colMeans(results, na.rm = TRUE),
                    apply(results, 2, sd))
  )
}



# ==============================================================================
# Scenarios
# ==============================================================================

scenarios <- expand.grid(
  n  = c(50, 100, 250, 500, 1000),
  m  = c(3, 10),
  p  = c(0.20, 0.45),
  r  = c(0, 0.35, 0.50),
  pi = c(1, 0.75),
  pi_h = c(0, 0.2, 0.75)
)
scenarios <- scenarios[
  (scenarios$m == 3  & scenarios$p == 0.20 & scenarios$r == 0.35 & scenarios$pi_h == 0) | 
  (scenarios$m == 3  & scenarios$p == 0.20 & scenarios$r == 0.35 & scenarios$pi_h == 0.2) |
  (scenarios$m == 3  & scenarios$p == 0.20 & scenarios$r == 0.35 & scenarios$pi_h == 0.75) |
  (scenarios$m == 10 & scenarios$p == 0.45 & scenarios$r == 0.50 & scenarios$pi_h == 0),
]
rownames(scenarios) <- NULL


# Szenarien unter MAR
scenarios_MAR <- expand.grid(
  n       = c(50, 100, 250, 500, 1000),
  m       = c(3, 10),
  p       = c(0.20, 0.45),
  r       = c(0.35, 0.50),
  pi_low  = 0.4,
  pi_high = 0.9
)
scenarios_MAR <- scenarios_MAR[
  (scenarios_MAR$m == 3  & scenarios_MAR$p == 0.20 & scenarios_MAR$r == 0.35) |
  (scenarios_MAR$m == 10 & scenarios_MAR$p == 0.45 & scenarios_MAR$r == 0.50),
]
rownames(scenarios_MAR) <- NULL

# ==============================================================================
# Run Simulations
# ==============================================================================

set.seed(123)

# Summary simulations across all scenarios
results <- apply(scenarios, 1, function(row) {
  simulation(
    n = row["n"], m = row["m"], p = row["p"],
    r = row["r"], pi = row["pi"], pi_h = row["pi_h"], n_reps = 1000
  )
})


set.seed(123)
results_MAR <- apply(scenarios_MAR, 1, function(row) {
  simulation_MAR(
    n       = row["n"],
    m       = row["m"],
    p       = row["p"],
    r       = row["r"],
    pi_low  = row["pi_low"],
    pi_high = row["pi_high"],
    n_reps  = 1000
  )
})

results_MAR_inv <- apply(scenarios_MAR, 1, function(row) {
  simulation_MAR_inv(
    n       = as.numeric(row["n"]),
    m       = as.numeric(row["m"]),
    p       = as.numeric(row["p"]),
    r       = as.numeric(row["r"]),
    pi_low  = as.numeric(row["pi_low"]),
    pi_high = as.numeric(row["pi_high"]),
    n_reps  = 1000
  )
})


# ==============================================================================
# Extract Results
# ==============================================================================

# Summary statistics per scenario (mean and sd rows)
sim_df <- do.call(rbind, lapply(seq_along(results), function(i) {
  s <- results[[i]]$summary
  data.frame(
    scenarios[i, ],
    type = "Simulation",
    mean_IOV   = s[1, "IOV"],
    sd_IOV     = s[2, "IOV"],
    lower_IOV  = s[1, "IOV"] - s[2, "IOV"],
    upper_IOV  = s[1, "IOV"] + s[2, "IOV"],
    mean_Skew  = s[1, "Skew"],
    sd_Skew    = s[2, "Skew"],
    lower_Skew = s[1, "Skew"] - s[2, "Skew"],
    upper_Skew = s[1, "Skew"] + s[2, "Skew"],
    mean_C     = s[1, "lag1_Cohen"],
    sd_C       = s[2, "lag1_Cohen"],
    lower_C    = s[1, "lag1_Cohen"] - s[2, "lag1_Cohen"],
    upper_C    = s[1, "lag1_Cohen"] + s[2, "lag1_Cohen"]
  )
}))

# Estimated CDFs per scenario

cdf_df <- do.call(rbind, lapply(which(scenarios$m == 3), function(i) {
  cdf <- results[[i]]$cdf
  data.frame(
    scenario  = i,
    statistic = c("mean", "sd"),
    scenarios[i, ],
    cdf
  )
}))

library(ggplot2)
library(dplyr)
library(tidyr)
# Wahre CDF




save.image("Masterarbeit.RData")
load("Masterarbeit.RData")