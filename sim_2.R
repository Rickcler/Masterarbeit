library(dplyr)
library(plyr)
library(tidyr)


# ==============================================================================
# Data Generating Processes
# ==============================================================================

#' Generate BinAR(1) process
#' @param n  Length of the time series
#' @param m  Binomial parameter
#' @param p  Mean of the stationary distribution E[I_t]
#' @param r  Thinning parameter (must be in [0, 1])
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
    marg_probs[i + 1] <- sum(data <= i) / n
  }
  return(marg_probs)
}

#' Empirical bivariate lag-h CDF (under missingness)
#' @param m     Number of categories minus 1
#' @param data  Full data vector
#' @param O     Missingness indicator vector
#' @param n     Length of the time series
#' @param h     Lag
biv_probs_e <- function(m, data, O, n, h = 1) {
  biv_probs <- numeric(m)
  valid <- (O[1:(n - h)] == 1) & (O[(h + 1):n] == 1)
  for (i in 0:(m - 1)) {
    biv_probs[i + 1] <- sum(
      (data[1:(n - h)] <= i) * (data[(h + 1):n] <= i) * valid
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

#' Simulation returning summary statistics and estimated CDF across replications
#' @param n      Time series length
#' @param m      Binomial parameter
#' @param p      Success probability
#' @param r      Thinning parameter
#' @param pi     Observation probability
#' @param n_reps Number of replications
#' @return List with two matrices (summary, cdf), each with rows: mean, sd
simulation <- function(n, m, p, r, pi, n_reps = 1000) {
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
    Missing_process <- generate_O(n, pi)
    observed_counts <- count_process[Missing_process == 1]
    CDF             <- marginal_probs_e(m, observed_counts, length(observed_counts))

    # i.i.d. process for Cohen's kappa
    iid_process         <- generate_iid(p, m, n)
    observed_counts_iid <- iid_process[Missing_process == 1]
    CDF_iid             <- marginal_probs_e(m, observed_counts_iid, length(observed_counts_iid))

    cdf_results[rep, ] <- CDF
    results[rep, 1]    <- (4 / m) * sum(CDF * (1 - CDF))
    results[rep, 2]    <- (2 / m) * sum(CDF - 1)
    results[rep, 3]    <- sum(biv_probs_e(m, iid_process, Missing_process, n, h = 1) - CDF_iid^2) /
                          sum(CDF_iid * (1 - CDF_iid))
    results[rep, 4]    <- sum(biv_probs_e(m, iid_process, Missing_process, n, h = 2) - CDF_iid^2) /
                          sum(CDF_iid * (1 - CDF_iid))
  }

  list(
    summary = rbind(colMeans(results,     na.rm = TRUE), apply(results,     2, sd)),
    cdf     = rbind(colMeans(cdf_results, na.rm = TRUE), apply(cdf_results, 2, sd))
  )
}

#' Simulation returning raw replication-level IOV and Skew estimates
#' @param n      Time series length
#' @param m      Binomial parameter
#' @param p      Success probability
#' @param r      Thinning parameter
#' @param pi     Observation probability
#' @param n_reps Number of replications
#' @return Data frame with one row per replication
simulation_raw <- function(n, m, p, r, pi, n_reps = 1000) {
  results <- data.frame(rep = 1:n_reps, IOV = NA, Skew = NA)

  for (rep in 1:n_reps) {
    count_process   <- generate_binar1(n, m, p, r)
    Missing_process <- generate_O(n, pi)
    observed_counts <- count_process[Missing_process == 1]
    CDF             <- marginal_probs_e(m, observed_counts, length(observed_counts))

    results$IOV[rep]  <- (4 / m) * sum(CDF * (1 - CDF))
    results$Skew[rep] <- (2 / m) * sum(CDF - 1)
  }

  results$n  <- n;  results$m  <- m
  results$p  <- p;  results$r  <- r
  results$pi <- pi

  return(results)
}


# ==============================================================================
# Scenarios
# ==============================================================================

scenarios <- expand.grid(
  n  = c(50, 100, 250, 500, 1000),
  m  = c(3, 10),
  p  = c(0.20, 0.45),
  r  = c(0, 0.35, 0.50),
  pi = c(1, 0.75)
)
scenarios <- scenarios[
  (scenarios$m == 3  & scenarios$p == 0.20 & scenarios$r == 0.35) |
  (scenarios$m == 10 & scenarios$p == 0.45 & scenarios$r == 0.50),
]
rownames(scenarios) <- NULL


# ==============================================================================
# Run Simulations
# ==============================================================================

set.seed(123)

# Summary simulations across all scenarios
results <- apply(scenarios, 1, function(row) {
  simulation(
    n = row["n"], m = row["m"], p = row["p"],
    r = row["r"], pi = row["pi"], n_reps = 1000
  )
})

# Raw replications for CLT / bias illustration
sim_10 <- simulation_raw(10,   10, 0.3, 0.2, 0.75, 2000)
sim_50   <- simulation_raw(50,   10, 0.3, 0.2, 0.75, 2000)
sim_200  <- simulation_raw(200,   10, 0.3, 0.2, 0.75, 2000)
sim_500  <- simulation_raw(200,   10, 0.3, 0.2, 0.75, 2000)
sim_1000 <- simulation_raw(1000, 10, 0.3, 0.2, 0.75, 2000)
sim_5000 <- simulation_raw(5000, 10, 0.3, 0.2, 0.75, 2000)

sim_data             <- rbind(sim_10, sim_50,  sim_1000)
sim_data_big         <- rbind(sim_50, sim_200, sim_500, sim_1000, sim_5000)
true_val             <- true_IOV(m = 10, p = 0.3)
sim_data$true_IOV    <- true_val
sim_data$diff        <- sim_data$IOV - true_val
sim_data$scaled_CLT  <- sqrt(sim_data$n) * sim_data$diff
sim_data$scaled_bias <- sim_data$n * sim_data$diff

sim_data_big$true_IOV    <- true_val
sim_data_big$diff        <- sim_data_big$IOV - true_val
sim_data_big$scaled_CLT  <- sqrt(sim_data_big$n) * sim_data_big$diff
sim_data_big$scaled_bias <- sim_data_big$n * sim_data_big$diff



# ==============================================================================
# Extract Results
# ==============================================================================

# Summary statistics per scenario (rows: mean, sd)
sim_df <- do.call(rbind, lapply(seq_along(results), function(i) {
  s <- results[[i]]$summary
  data.frame(
    scenarios[i, ],
    type       = "Simulated",
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

# Estimated CDFs for m = 3 scenarios only (consistent column count)
cdf_df <- do.call(rbind, lapply(which(scenarios$m == 3), function(i) {
  cdf <- results[[i]]$cdf
  data.frame(
    scenario  = i,
    statistic = c("mean", "sd"),
    scenarios[i, , drop = FALSE],
    cdf,
    row.names = NULL
  )
}))

save.image("Masterarbeit.RData")
