# ==============================================================================
# 01_functions.R
# Alle reinen Funktionen: DGP, Schätzer, asymptotische Theorie
# ==============================================================================

# ------------------------------------------------------------------------------
# Data Generating Processes
# ------------------------------------------------------------------------------

#' Generate BinAR(1) process
#' @param n    Length of the time series
#' @param m    Binomial parameter
#' @param p    Success probability (= stationary mean / m)
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

#' Generate i.i.d. Binomial process
#' @param n  Length of the time series
#' @param m  Binomial parameter
#' @param p  Success probability
generate_iid <- function(n, m, p) {
  rbinom(n, m, p)
}

#' Generate MCAR missingness indicator
#' @param n   Length of the time series
#' @param pi  Observation probability
generate_O <- function(n, pi) {
  rbinom(n, 1, pi)
}

#' Generate MAR missingness: observation probability increases with X_t
#' @param x        Realized time series
#' @param m        Binomial parameter
#' @param pi_low   Observation probability when X_t = 0
#' @param pi_high  Observation probability when X_t = m
generate_O_MAR <- function(x, m, pi_low = 0.5, pi_high = 0.9) {
  pi_t <- pi_low + (pi_high - pi_low) * (x / m)
  rbinom(length(x), 1, pi_t)
}

#' Generate MAR missingness (inverted): observation probability decreases with X_t
#' @param x        Realized time series
#' @param m        Binomial parameter
#' @param pi_low   Observation probability when X_t = m
#' @param pi_high  Observation probability when X_t = 0
generate_O_MAR_inv <- function(x, m, pi_low = 0.5, pi_high = 0.9) {
  pi_t <- pi_high - (pi_high - pi_low) * (x / m)
  rbinom(length(x), 1, pi_t)
}


# ------------------------------------------------------------------------------
# Estimation Functions
# ------------------------------------------------------------------------------

#' Empirical marginal CDF at 0, ..., m-1
#' @param m     Number of categories minus 1
#' @param data  Observed data vector
#' @param n     Number of observations (used as denominator)
marginal_probs_e <- function(m, data, n) {
  marg_probs <- numeric(m)
  for (i in 0:(m - 1)) {
    marg_probs[i + 1] <- sum(data <= i) / n
  }
  return(marg_probs)
}

#' Empirical bivariate lag-h CDF under pairwise-complete missingness
#' @param m     Number of categories minus 1
#' @param data  Full (partially unobserved) data vector
#' @param O     Missingness indicator vector (1 = observed)
#' @param n     Length of the time series
#' @param h     Lag
biv_probs_e <- function(m, data, O, n, h = 1) {
  biv_probs <- numeric(m)
  valid     <- (O[1:(n - h)] == 1) & (O[(h + 1):n] == 1)
  denom     <- sum(valid)
  for (i in 0:(m - 1)) {
    biv_probs[i + 1] <- sum(
      (data[1:(n - h)] <= i) * (data[(h + 1):n] <= i) * valid
    ) / denom
  }
  return(biv_probs)
}


# ------------------------------------------------------------------------------
# True Parameter Functions
# ------------------------------------------------------------------------------

#' True IOV under BinAR(1) stationary distribution
#' @param m  Binomial parameter
#' @param p  Success probability
true_IOV <- function(m, p) {
  F <- pbinom(0:(m - 1), m, p)
  (4 / m) * sum(F * (1 - F))
}

#' True Skewness under BinAR(1) stationary distribution
#' @param m  Binomial parameter
#' @param p  Success probability
true_Skew <- function(m, p) {
  F <- pbinom(0:(m - 1), m, p)
  (2 / m) * sum(F - 1)
}


# ------------------------------------------------------------------------------
# Joint Distribution Functions (BinAR(1))
# ------------------------------------------------------------------------------

#' Conditional distribution P(X_{t+h} = j | X_t = i) for BinAR(1)
#' @param m  Binomial parameter
#' @param p  Success probability
#' @param r  Thinning parameter
#' @param h  Lag
lag_h_conditional <- function(m, p, r, h = 1) {
  beta  <- p * (1 - r^h)
  alpha <- beta + r^h
  res   <- matrix(0, nrow = m + 1, ncol = m + 1)

  for (i in 0:m) {
    for (j in 0:m) {
      smaller <- min(i, j)
      bigger  <- max(0, i + j - m)
      res[i + 1, j + 1] <- sum(sapply(bigger:smaller, function(k) {
        choose(i, k) * choose(m - i, j - k) *
          alpha^k * (1 - alpha)^(i - k) *
          beta^(j - k) * (1 - beta)^(m - i - j + k)
      }))
    }
  }
  return(res)
}

#' Joint PMF of (X_t, X_{t+h}) for BinAR(1)
#' @param m  Binomial parameter
#' @param p  Success probability
#' @param r  Thinning parameter
#' @param h  Lag
lag_h_joint_pmf <- function(m, p, r, h = 1) {
  lag_h_conditional(m, p, r, h = h) * dbinom(0:m, m, p)
}

#' Convert joint PMF matrix to joint CDF matrix
#' @param pmf  (m+1) x (m+1) matrix of joint probabilities
pmf_to_cdf <- function(pmf) {
  cdf <- apply(pmf, 2, cumsum)
  cdf <- t(apply(cdf, 1, cumsum))
  return(cdf)
}

#' Joint CDF of (X_t, X_{t+h}) for BinAR(1)
#' @param m  Binomial parameter
#' @param p  Success probability
#' @param r  Thinning parameter
#' @param h  Lag
lag_h_joint_cdf <- function(m, p, r, h = 1) {
  pmf_to_cdf(lag_h_joint_pmf(m, p, r, h = h))
}


# ------------------------------------------------------------------------------
# Long-run Covariance Matrices
# ------------------------------------------------------------------------------


#' @param m     Binomial parameter of the count process
#' @param p     Success probability of the count process
#' @param r     Thinning parameter of the count process
#' @param pi    Marginal observation probability  P(R_t = 1)
#' @param pi_h  Lag-1 joint probability P(R_t = 1, R_{t+1} = 1); starting value
#'              for the missingness BinAR(1) recursion
#' @param H     Truncation lag (default 50)
Sigma_Star <- function(m, p, r, pi, pi_h, H = 50) {
  f     <- pbinom(0:(m - 1), m, p)
  Sigma <- matrix(0, nrow = m, ncol = m)

  for (i in 0:(m - 1)) {
    for (j in 0:(m - 1)) {
      smaller  <- min(i, j)
      iid_part <- (1 / pi) * (f[smaller + 1] - f[i + 1] * f[j + 1])

      lag_sum      <- 0
      pi_h_current <- pi_h   # reset for each (i,j) — fixes overwrite bug

      for (h in 1:H) {
        cdf_lag_h    <- lag_h_joint_cdf(m, p, r, h)
        pi_h_current <- lag_h_joint_pmf(1, pi, pi_h_current, h)[2, 2]
        lag_sum      <- lag_sum +
          pi_h_current * (
            cdf_lag_h[i + 1, j + 1] + cdf_lag_h[j + 1, i + 1] -
            2 * f[i + 1] * f[j + 1]
          )
      }

      Sigma[i + 1, j + 1] <- iid_part + (1 / pi^2) * lag_sum
    }
  }
  return(Sigma)
}


# ------------------------------------------------------------------------------
# Asymptotic Distribution Functions
# ------------------------------------------------------------------------------

#' Asymptotic expectation and standard deviation of the IOV estimator
#' @param Sigma        Long-run covariance matrix (m x m)
#' @param marginal_cdf Marginal CDF vector of length m
#' @param n            Sample size
#' @param m            Binomial parameter
#' @return Named numeric(2): c(expectation, sd)
IOV_asymptotic <- function(Sigma, marginal_cdf, n, m) {
  variance_sum <- 0
  for (i in 1:m) {
    for (j in 1:m) {
      variance_sum <- variance_sum +
        (1 - 2 * marginal_cdf[i]) * (1 - 2 * marginal_cdf[j]) * Sigma[i, j]
    }
  }

  Var_IOV        <- (1 / n) * (16 / m^2) * variance_sum
  IOV_real       <- (4 / m) * sum(marginal_cdf * (1 - marginal_cdf))
  Sigma_trace    <- sum(diag(Sigma))
  Expectation_IOV <- IOV_real - (1 / n) * (4 / m) * Sigma_trace

  c(expectation = Expectation_IOV, sd = sqrt(Var_IOV))
}

#' Asymptotic expectation and standard deviation of the Skewness estimator
#' @param Sigma        Long-run covariance matrix (m x m)
#' @param marginal_cdf Marginal CDF vector of length m
#' @param n            Sample size
#' @param m            Binomial parameter
#' @return Named numeric(2): c(expectation, sd)
Skew_asymptotic <- function(Sigma, marginal_cdf, n, m) {
  skew_real     <- (2 / m) * sum(marginal_cdf - 1)
  skew_variance <- (1 / n) * (4 / m^2) * sum(Sigma)
  c(expectation = skew_real, sd = sqrt(skew_variance))
}

#' Asymptotic expectation and standard deviation of Cohen's kappa (i.i.d. case)
#' @param n            Sample size
#' @param pi           Observation probability
#' @param m            Binomial parameter
#' @param marginal_cdf Marginal CDF vector of length m
#' @return Named numeric(2): c(expectation, sd)
Cohens_asymptotic_iid <- function(n, pi, m, marginal_cdf) {
  t_1 <- 0
  t_2 <- 0
  t_3 <- 0
  f   <- marginal_cdf

  for (i in 0:(m - 1)) {
    for (j in 0:(m - 1)) {
      smaller <- min(i, j)
      t_1     <- t_1 + (f[smaller + 1] - f[i + 1] * f[j + 1])^2
      t_2     <- t_2 + (f[i + 1] * f[j + 1] * (f[smaller + 1] - f[i + 1] * f[j + 1]))
    }
    t_3 <- t_3 + (f[i + 1] * (1 - f[i + 1]))
  }
  t_3 <- t_3^2

  Var_Cohens_K       <- ((1 / n) * (1 / pi^2) * (t_1 / t_3)) +
                        ((2 / n) * ((1 - pi) / pi^2) * (t_2 / t_3))
  Expectation_Cohens_K <- -(1 / (n * pi))

  c(expectation = Expectation_Cohens_K, sd = sqrt(Var_Cohens_K))
}


# ------------------------------------------------------------------------------
# Simulation Functions
# ------------------------------------------------------------------------------

#' Simulation: summary statistics across replications (MCAR / serially dependent)
#' @param n      Time series length
#' @param m      Binomial parameter
#' @param p      Success probability
#' @param r      Thinning parameter
#' @param pi     Marginal observation probability
#' @param pi_h   Lag-1 joint observation probability (0 = i.i.d. MCAR)
#' @param n_reps Number of replications
#' @return List with $summary (mean/sd x 4 statistics) and $cdf (mean/sd x m)
simulation <- function(n, m, p, r, pi, pi_h, n_reps = N_REPS) {
  results <- matrix(
    NA, nrow = n_reps, ncol = 4,
    dimnames = list(NULL, c("IOV", "Skew", "lag1_Cohen", "lag2_Cohen"))
  )
  cdf_results <- matrix(
    NA, nrow = n_reps, ncol = m,
    dimnames = list(NULL, paste0("f_", 0:(m - 1)))
  )

  for (rep in 1:n_reps) {
    count_process   <- generate_binar1(n, m, p, r)
    Missing_process <- generate_binar1(n, 1, pi, pi_h)
    observed_counts <- count_process[Missing_process == 1]
    CDF             <- marginal_probs_e(m, observed_counts, length(observed_counts))

    # i.i.d. process for Cohen's kappa
    iid_process         <- generate_iid(n, m, p)
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
    summary = rbind(
      mean = colMeans(results,     na.rm = TRUE),
      sd   = apply(results,     2, sd, na.rm = TRUE)
    ),
    cdf = rbind(
      mean = colMeans(cdf_results, na.rm = TRUE),
      sd   = apply(cdf_results, 2, sd, na.rm = TRUE)
    )
  )
}

#' Simulation under MAR (observation probability increases with X_t)
#' @param n        Time series length
#' @param m        Binomial parameter
#' @param p        Success probability
#' @param r        Thinning parameter
#' @param pi_low   Observation probability when X_t = 0
#' @param pi_high  Observation probability when X_t = m
#' @param n_reps   Number of replications
#' @return List with $summary (mean/sd x 4 statistics)
simulation_MAR <- function(n, m, p, r, pi_low = 0.5, pi_high = 0.9,
                            n_reps = N_REPS) {
  results <- matrix(
    NA, nrow = n_reps, ncol = 4,
    dimnames = list(NULL, c("IOV", "Skew", "lag1_Cohen", "lag2_Cohen"))
  )

  for (rep in 1:n_reps) {
    count_process   <- generate_binar1(n, m, p, r)
    Missing_process <- generate_O_MAR(count_process, m, pi_low, pi_high)
    observed_counts <- count_process[Missing_process == 1]
    CDF             <- marginal_probs_e(m, observed_counts, length(observed_counts))

    iid_process         <- generate_iid(n, m, p)
    Missing_iid         <- generate_O_MAR(iid_process, m, pi_low, pi_high)
    observed_counts_iid <- iid_process[Missing_iid == 1]
    CDF_iid             <- marginal_probs_e(m, observed_counts_iid, length(observed_counts_iid))

    results[rep, 1] <- (4 / m) * sum(CDF * (1 - CDF))
    results[rep, 2] <- (2 / m) * sum(CDF - 1)
    results[rep, 3] <- sum(biv_probs_e(m, iid_process, Missing_iid, n, h = 1) - CDF_iid^2) /
                       sum(CDF_iid * (1 - CDF_iid))
    results[rep, 4] <- sum(biv_probs_e(m, iid_process, Missing_iid, n, h = 2) - CDF_iid^2) /
                       sum(CDF_iid * (1 - CDF_iid))
  }

  list(
    summary = rbind(
      mean = colMeans(results, na.rm = TRUE),
      sd   = apply(results, 2, sd, na.rm = TRUE)
    )
  )
}

#' Simulation under inverted MAR (observation probability decreases with X_t)
#' @param n        Time series length
#' @param m        Binomial parameter
#' @param p        Success probability
#' @param r        Thinning parameter
#' @param pi_low   Observation probability when X_t = m
#' @param pi_high  Observation probability when X_t = 0
#' @param n_reps   Number of replications
#' @return List with $summary (mean/sd x 4 statistics)
simulation_MAR_inv <- function(n, m, p, r, pi_low = 0.5, pi_high = 0.9,
                                n_reps = N_REPS) {
  results <- matrix(
    NA, nrow = n_reps, ncol = 4,
    dimnames = list(NULL, c("IOV", "Skew", "lag1_Cohen", "lag2_Cohen"))
  )

  for (rep in 1:n_reps) {
    count_process   <- generate_binar1(n, m, p, r)
    Missing_process <- generate_O_MAR_inv(count_process, m, pi_low, pi_high)
    observed_counts <- count_process[Missing_process == 1]
    CDF             <- marginal_probs_e(m, observed_counts, length(observed_counts))

    iid_process         <- generate_iid(n, m, p)
    Missing_iid         <- generate_O_MAR_inv(iid_process, m, pi_low, pi_high)
    observed_counts_iid <- iid_process[Missing_iid == 1]
    CDF_iid             <- marginal_probs_e(m, observed_counts_iid, length(observed_counts_iid))

    results[rep, 1] <- (4 / m) * sum(CDF * (1 - CDF))
    results[rep, 2] <- (2 / m) * sum(CDF - 1)
    results[rep, 3] <- sum(biv_probs_e(m, iid_process, Missing_iid, n, h = 1) - CDF_iid^2) /
                       sum(CDF_iid * (1 - CDF_iid))
    results[rep, 4] <- sum(biv_probs_e(m, iid_process, Missing_iid, n, h = 2) - CDF_iid^2) /
                       sum(CDF_iid * (1 - CDF_iid))
  }

  list(
    summary = rbind(
      mean = colMeans(results, na.rm = TRUE),
      sd   = apply(results, 2, sd, na.rm = TRUE)
    )
  )
}

#' Simulation (raw): returns one row per replication for CLT / bias illustrations
#' @param n      Time series length
#' @param m      Binomial parameter
#' @param p      Success probability
#' @param r      Thinning parameter
#' @param pi     Observation probability
#' @param n_reps Number of replications
#' @return Data frame with columns rep, IOV, Skew, n, m, p, r, pi
simulation_raw <- function(n, m, p, r, pi, n_reps = N_REPS) {
  results <- data.frame(rep = 1:n_reps, IOV = NA_real_, Skew = NA_real_)

  for (rep in 1:n_reps) {
    count_process   <- generate_binar1(n, m, p, r)
    Missing_process <- generate_O(n, pi)
    observed_counts <- count_process[Missing_process == 1]
    CDF             <- marginal_probs_e(m, observed_counts, length(observed_counts))

    results$IOV[rep]  <- (4 / m) * sum(CDF * (1 - CDF))
    results$Skew[rep] <- (2 / m) * sum(CDF - 1)
  }

  results$n  <- n
  results$m  <- m
  results$p  <- p
  results$r  <- r
  results$pi <- pi
  return(results)
}
