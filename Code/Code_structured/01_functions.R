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

#' Compute rejection rate for kappa_ord test at alpha = 0.05
#' under H_A for a given scenario and n_grid
#' @param n       Time series length
#' @param m       Binomial parameter
#' @param p       Success probability
#' @param r       Thinning parameter (> 0 for H_A)
#' @param pi      Observation probability
#' @param h       Lag
#' @param alpha   Significance level
#' @param n_reps  Number of replications
rejection_rate <- function(n, m, p, r, pi, h = 1,
                            alpha = 0.05, n_reps = 1000) {
  # Kritischer Wert aus asymptotischer Varianz unter H_0
  sd_H0    <- sqrt(asymp_var_kappa_H0(m, p, pi) / n)
  crit_val <- qnorm(1 - alpha / 2) * sd_H0

  kappa_vals <- simulation_kappa_HA(n, m, p, r, pi,
                                     h = h, n_reps = n_reps)
  mean(abs(kappa_vals) > crit_val, na.rm = TRUE)
}


#' Rejection rate from raw simulation results
#' @param kappa_vals  Vector of simulated kappa values
#' @param crit        Critical value
rejection_from_vals <- function(kappa_vals, crit) {
  mean(abs(kappa_vals) > crit, na.rm = TRUE)
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

# True kappa unter H_A (aus geschlossener Form)
true_kappa_HA <- function(m, p, r, h = 1) {
  f    <- pbinom(0:(m - 1), size = m, prob = p)
  cdf  <- lag_h_joint_cdf(m, p, r, h = h)
  f_ii <- diag(cdf)[1:m]  # f_{ii}(h) fuer i = 0,...,m-1
  sum(f_ii - f^2) / sum(f * (1 - f))
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


# Computation of Cohens ordinal Kappa for given h
# κ_ord(h) = (∑_{j=0}^{m-1} (f_{jj}(h) - f_j^2)) / (∑_{i=0}^{m-1} f_i (1-f_i))
#' @param m  Binomial parameter
#' @param p  Success probability
#' @param r  Thinning parameter
#' @param h  Lag
kappa_ord <- function(m, p, r, h) {
  # gemeinsame CDF für Lag h
  joint_cdf <- lag_h_joint_cdf(m, p, r, h)
  # Randwahrscheinlichkeiten (stationär, daher gleich)
  marg_cdf <- pbinom(0:m, m, p)
  # Summe über j = 0 ... m-1
  idx <- 1:m   # entspricht Werten 0,1,...,m-1
  f_jj <- sapply(idx, function(j) joint_cdf[j, j])   # Diagonalelemente
  f_j  <- marg_cdf[idx]                              # Randwahrscheinlichkeiten
  # Zähler und Nenner
  numerator   <- sum(f_jj - f_j^2)
  denominator <- sum(f_j * (1 - f_j))
  return(numerator / denominator)
}

# ------------------------------------------------------------------------------
# Long-run Covariance Matrices
# ------------------------------------------------------------------------------


#' @param m     Binomial parameter of the count process
#' @param p     Success probability of the count process
#' @param r     Thinning parameter of the count process
#' @param pi    Marginal observation probability  P(R_t = 1)
#' @param r_pi  Lag-1 joint probability P(R_t = 1, R_{t+1} = 1); starting value
#'              for the missingness BinAR(1) recursion
#' @param H     Truncation lag (default 50)
Sigma_Star <- function(m, p, r, pi, r_pi, H = 50) {
  f     <- pbinom(0:(m - 1), m, p)
  Sigma <- matrix(0, nrow = m, ncol = m)

  for (i in 0:(m - 1)) {
    for (j in 0:(m - 1)) {
      smaller  <- min(i, j)
      iid_part <- (1 / pi) * (f[smaller + 1] - f[i + 1] * f[j + 1])

      lag_sum      <- 0
      pi_h_current <- r_pi   # reset for each (i,j) — fixes overwrite bug

      for (h in 1:H) {
        cdf_lag_h    <- lag_h_joint_cdf(m, p, r, h)
        r_pi_current <- lag_h_joint_pmf(1, pi, pi_h_current, h)[2, 2]
        lag_sum      <- lag_sum +
          r_pi_current * (
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


#' Asymptotic variance of kappa_ord under H_0 and i.i.d. O_t
#' @param m   Binomial parameter
#' @param p   Success probability
#' @param pi  Observation probability
asymp_var_kappa_H0 <- function(m, p, pi) {
  f <- pbinom(0:(m - 1), size = m, prob = p)
  B <- sum(f * (1 - f))

  # Doppelsumme Term 1: (f_min{i,j} - f_i f_j)^2
  term1 <- 0
  term2 <- 0
  for (i in 0:(m - 1)) {
    for (j in 0:(m - 1)) {
      fmin   <- f[min(i, j) + 1]
      fi     <- f[i + 1]
      fj     <- f[j + 1]
      delta  <- fmin - fi * fj
      term1  <- term1 + delta^2
      term2  <- term2 + fi * fj * delta
    }
  }

  var_kappa <- (1 / pi^2) * term1 / B^2 +
               2 * ((1 - pi) / pi^2) * term2 / B^2
  return(var_kappa)
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
#' @param r_pi   Lag-1 joint observation probability (0 = i.i.d. MCAR)
#' @param n_reps Number of replications
#' @return List with $summary (mean/sd x 4 statistics) and $cdf (mean/sd x m)
simulation <- function(n, m, p, r, pi, r_pi, n_reps = N_REPS) {
  results <- matrix(
    NA, nrow = n_reps, ncol = 6,
    dimnames = list(NULL, c("IOV", "Skew", "lag1_Cohen", "lag2_Cohen", "lag1_Cohen_bc", "lag2_Cohen_bc"))
  )
  cdf_results <- matrix(
    NA, nrow = n_reps, ncol = m,
    dimnames = list(NULL, paste0("f_", 0:(m - 1)))
  )

  for (rep in 1:n_reps) {
    count_process   <- generate_binar1(n, m, p, r)
    Missing_process <- generate_binar1(n, 1, pi, r_pi)
    observed_counts <- count_process[Missing_process == 1]
    CDF             <- marginal_probs_e(m, observed_counts, length(observed_counts))

    # i.i.d. process for Cohen's kappa
    iid_process         <- generate_iid(n, m, p)
    observed_counts_iid <- iid_process[Missing_process == 1]
    CDF_iid             <- marginal_probs_e(m, observed_counts_iid, length(observed_counts_iid))
    pi_est <- mean(Missing_process)
    cdf_results[rep, ] <- CDF
    results[rep, 1]    <- (4 / m) * sum(CDF * (1 - CDF))
    results[rep, 2]    <- (2 / m) * sum(CDF - 1)
    results[rep, 3]    <- sum(biv_probs_e(m, iid_process, Missing_process, n, h = 1) - CDF_iid^2) /
                          sum(CDF_iid * (1 - CDF_iid))
    results[rep, 4]    <- sum(biv_probs_e(m, iid_process, Missing_process, n, h = 2) - CDF_iid^2) /
                          sum(CDF_iid * (1 - CDF_iid))
    results[rep, 5]    <- results[rep, 3] + (1/(pi_est * n)) # Bias correction for kappa
    results[rep, 6]    <- results[rep, 4] + (1/(pi_est * n)) # Bias correction for kappa
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


#' Simulation: returns raw kappa estimates under H_A (r > 0)
#' @param n      Time series length
#' @param m      Binomial parameter
#' @param p      Success probability
#' @param r      Thinning parameter (> 0 for H_A)
#' @param pi     Observation probability
#' @param h      Lag for kappa
#' @param n_reps Number of replications
simulation_kappa_HA <- function(n, m, p, r, pi, h = 1, n_reps = 1000) {
  kappa_vals <- numeric(n_reps)

  for (rep in 1:n_reps) {
    count_process   <- generate_binar1(n, m, p, r)
    Missing_process <- generate_O(n, pi)

    observed        <- count_process[Missing_process == 1]
    CDF_obs         <- marginal_probs_e(m, observed, length(observed))
    biv             <- biv_probs_e(m, count_process, Missing_process, n, h = h)

    denom <- sum(CDF_obs * (1 - CDF_obs))
    if (denom > 0) {
      kappa_vals[rep] <- sum(biv - CDF_obs^2) / denom
    } else {
      kappa_vals[rep] <- NA
    }
  }
  return(kappa_vals)
}


# ------------------------------------------------------------------------------
# Plot Functions
# ------------------------------------------------------------------------------

#' Generischer Vergleichs-Plot (IOV, Skew oder Cohen's K)
#' @param subset        Vorbereitetes subset mit x_pos
#' @param y_var         String: Spaltenname der y-Variable (z.B. "mean_IOV")
#' @param ymin_var      String: untere Fehlerbalken-Spalte
#' @param ymax_var      String: obere Fehlerbalken-Spalte
#' @param true_val      Wahrer Parameterwert (horizontale Linie)
#' @param group_centers Mittelpunkte der x-Gruppen
#' @param x_labels      Beschriftungen der x-Gruppen
#' @param title         Plot-Titel
#' @param y_label       y-Achsenbeschriftung
#' @param subtitle      Plot-Untertitel
#' @param ylim_offset   c(unten, oben) relativ zum wahren Wert
#' @param group_var     String: Variable für Shape-Ästhetik
comparison_plot <- function(subset, y_var, ymin_var, ymax_var,
                            true_val, group_centers, x_labels,
                            title, y_label, subtitle,
                            ylim_offset = c(-0.125, 0.075),
                            group_var = "pi") {
  legend_label <- switch(
    group_var,
    "pi"   = expression(pi),
    "r_pi" = expression(r[pi]),
    "pi_h" = expression(pi[h]),
    group_var
  )
  n_groups   <- length(group_centers)
  rect_xmin  <- c(0.5, seq(2.5, by = 2, length.out = n_groups - 1))
  rect_xmax  <- c(seq(2.5, by = 2, length.out = n_groups - 1),
                  max(subset$x_pos) + 0.5)

  ggplot(subset, aes(x = x_pos, y = .data[[y_var]],
                     color = type,
                     shape = factor(.data[[group_var]]))) +
    annotate("rect",
             xmin  = rect_xmin[seq(1, n_groups, by = 2)],
             xmax  = rect_xmax[seq(1, n_groups, by = 2)],
             ymin  = -Inf, ymax = Inf,
             alpha = 0.05, fill = "gray90") +
    geom_hline(yintercept = true_val,
               color = "darkgreen", linetype = "dashed",
               linewidth = 1, alpha = 0.7) +
    geom_line(aes(group = interaction(n, .data[[group_var]])),
              color = "gray50", linetype = "dashed",
              alpha = 0.5,
              position = position_dodge(width = 0.2)) +
    geom_point(size = 3.5, position = position_dodge(width = 0.2)) +
    geom_errorbar(aes(ymin = .data[[ymin_var]], ymax = .data[[ymax_var]]),
                  width = 0.15, linewidth = 0.8,
                  position = position_dodge(width = 0.2)) +
    scale_x_continuous(breaks = group_centers, labels = x_labels,
                       expand = expansion(mult = 0.1)) +
    scale_color_manual(
      values = c("Asymptotic" = "#E41A1C", "Simulation" = "#377EB8"),
      name   = "Method"
    ) +
    scale_shape_manual(
      values = setNames(c(16, 17, 15, 18),
                        as.character(sort(unique(subset[[group_var]])))),
      name   = legend_label
    ) +
    coord_cartesian(ylim = c(true_val + ylim_offset[1],
                             true_val + ylim_offset[2])) +
    labs(title = title, subtitle = subtitle,
         x = "", y = "") +
    theme_minimal() +
    theme(
      plot.title         = element_text(hjust = 0.5, face = "bold", size = 20),
      plot.subtitle      = element_text(hjust = 0.5, color = "gray40"),
      axis.text.x        = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 12, color = "gray20"),
      axis.text.y        = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 14, color = "gray20"),
      legend.position    = "bottom",
      legend.box         = "horizontal",      # Legendenblöcke nebeneinander
      legend.direction   = "horizontal",        # Items innerhalb je vertikal
      legend.spacing.x   = unit(1, "cm"),     # Abstand zwischen den Blöcken
      legend.text        = element_text(size = 18),
      legend.title       = element_text(size = 20, face = "bold"),
      legend.key.size    = unit(0.6, "cm"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.border       = element_rect(color = "gray80", fill = NA,
                                        linewidth = 0.5)
    ) +
    guides(
      color = guide_legend(order = 1, override.aes = list(size = 4)),
      shape = guide_legend(order = 2, override.aes = list(size = 4))
    )
}
