library(dplyr)
library(tidyr)


# ==============================================================================
# Theoretical Joint Distribution Functions
# ==============================================================================

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
      res[i + 1, j + 1] <- sum(sapply(bigger:smaller, function(n) {
        choose(i, n) * choose(m - i, j - n) *
          alpha^n * (1 - alpha)^(i - n) *
          beta^(j - n) * (1 - beta)^(m - i - j + n)
      }))
    }
  }
  return(res)
}

#' Joint PMF matrix P(X_t = i, X_{t+h} = j) for BinAR(1)
#' @param m  Binomial parameter
#' @param p  Success probability
#' @param r  Thinning parameter
#' @param h  Lag
lag_h_joint_pmf <- function(m, p, r, h = 1) {
  lag_h_conditional(m, p, r, h = h) * dbinom(0:m, m, p)
}

#' Convert joint PMF matrix to joint CDF matrix
#' @param pmf  Joint PMF matrix
pmf_to_cdf <- function(pmf) {
  cdf <- apply(pmf, 2, cumsum)
  cdf <- t(apply(cdf, 1, cumsum))
  return(cdf)
}

#' Joint CDF matrix F(i, j; h) = P(X_t <= i, X_{t+h} <= j) for BinAR(1)
#' @param m  Binomial parameter
#' @param p  Success probability
#' @param r  Thinning parameter
#' @param h  Lag
lag_h_joint_cdf <- function(m, p, r, h = 1) {
  pmf_to_cdf(lag_h_joint_pmf(m, p, r, h = h))
}


# ==============================================================================
# Asymptotic Covariance Matrix
# ==============================================================================

#' Long-run covariance matrix Sigma* of the CDF estimator under missingness
#' @param m   Binomial parameter
#' @param p   Success probability
#' @param r   Thinning parameter
#' @param pi  Observation probability
#' @param H   Lag truncation (default 50)
Sigma_Star <- function(m, p, r, pi, H = 50) {
  f     <- pbinom(0:(m - 1), m, p)
  Sigma <- matrix(0, nrow = m, ncol = m)

  for (i in 0:(m - 1)) {
    for (j in 0:(m - 1)) {
      smaller  <- min(i, j)
      iid_part <- (1 / pi) * (f[smaller + 1] - f[i + 1] * f[j + 1])

      lag_sum <- 0
      for (h in 1:H) {
        cdf_h   <- lag_h_joint_cdf(m, p, r, h)
        lag_sum <- lag_sum + cdf_h[i + 1, j + 1] + cdf_h[j + 1, i + 1] - 2 * f[i + 1] * f[j + 1]
      }

      Sigma[i + 1, j + 1] <- iid_part + lag_sum
    }
  }
  return(Sigma)
}


# ==============================================================================
# Asymptotic Distribution Functions
# ==============================================================================

#' Asymptotic expectation and standard deviation of the IOV estimator
#' @param Sigma        Long-run covariance matrix
#' @param marginal_cdf Marginal CDF vector f_0, ..., f_{m-1}
#' @param n            Sample size
#' @param m            Binomial parameter
#' @return c(expectation, sd)
IOV_asymptotic <- function(Sigma, marginal_cdf, n, m) {
  variance_sum <- sum(
    outer(1 - 2 * marginal_cdf, 1 - 2 * marginal_cdf) * Sigma
  )
  Var_IOV          <- (16 / (n * m^2)) * variance_sum
  IOV_true         <- (4 / m) * sum(marginal_cdf * (1 - marginal_cdf))
  Expectation_IOV  <- IOV_true - (4 / (n * m)) * sum(diag(Sigma))
  return(c(Expectation_IOV, sqrt(Var_IOV)))
}

#' Asymptotic expectation and standard deviation of the Skew estimator
#' @param Sigma        Long-run covariance matrix
#' @param marginal_cdf Marginal CDF vector
#' @param n            Sample size
#' @param m            Binomial parameter
#' @return c(expectation, sd)
Skew_asymptotic <- function(Sigma, marginal_cdf, n, m) {
  skew_true     <- (2 / m) * sum(marginal_cdf - 1)
  skew_variance <- (4 / (n * m^2)) * sum(Sigma)
  return(c(skew_true, sqrt(skew_variance)))
}

#' Asymptotic expectation and standard deviation of Cohen's kappa (i.i.d. case)
#' @param n            Sample size
#' @param pi           Observation probability
#' @param m            Binomial parameter
#' @param marginal_cdf Marginal CDF vector
#' @return c(expectation, sd)
Cohens_asymptotic_iid <- function(n, pi, m, marginal_cdf) {
  f   <- marginal_cdf
  t_1 <- 0; t_2 <- 0; t_3 <- 0

  for (i in 0:(m - 1)) {
    for (j in 0:(m - 1)) {
      smaller <- min(i, j)
      t_1 <- t_1 + (f[smaller + 1] - f[i + 1] * f[j + 1])^2
      t_2 <- t_2 + (f[i + 1] * f[j + 1] * (f[smaller + 1] - f[i + 1] * f[j + 1]))
    }
    t_3 <- t_3 + f[i + 1] * (1 - f[i + 1])
  }
  t_3 <- t_3^2

  Var_Cohens_K      <- (1 / (n * pi^2)) * (t_1 / t_3) + (2 * (1 - pi) / (n * pi^2)) * (t_2 / t_3)
  Expectation_K     <- -(1 / (n * pi))
  return(c(Expectation_K, sqrt(Var_Cohens_K)))
}


# ==============================================================================
# Scenarios and Asymptotic Results
# ==============================================================================

unique_coeffs <- data.frame(rbind(
  c(m = 3,  p = 0.20, r = 0.35, pi = 1),
  c(m = 10, p = 0.45, r = 0.50, pi = 1),
  c(m = 3,  p = 0.20, r = 0.35, pi = 0.75),
  c(m = 10, p = 0.45, r = 0.50, pi = 0.75)
))
unique_n <- c(50, 100, 250, 500, 1000)

asymp_df <- do.call(rbind, lapply(1:nrow(unique_coeffs), function(idx) {
  params  <- unique_coeffs[idx, ]
  m_val   <- as.numeric(params["m"])
  p_val   <- as.numeric(params["p"])
  r_val   <- as.numeric(params["r"])
  pi_val  <- as.numeric(params["pi"])

  marginal <- pbinom(0:(m_val - 1), m_val, p_val)
  Sigma    <- Sigma_Star(m_val, p_val, r_val, pi_val)

  do.call(rbind, lapply(unique_n, function(n_val) {
    iov_res  <- IOV_asymptotic(Sigma, marginal, n_val, m_val)
    skew_res <- Skew_asymptotic(Sigma, marginal, n_val, m_val)
    C_res    <- Cohens_asymptotic_iid(n_val, pi_val, m_val, marginal)

    data.frame(
      n          = n_val,  pi   = pi_val,
      m          = m_val,  p    = p_val,  r = r_val,
      type       = "Asymptotic",
      mean_IOV   = iov_res[1],
      sd_IOV     = iov_res[2],
      lower_IOV  = iov_res[1]  - iov_res[2],
      upper_IOV  = iov_res[1]  + iov_res[2],
      mean_Skew  = skew_res[1],
      sd_Skew    = skew_res[2],
      lower_Skew = skew_res[1] - skew_res[2],
      upper_Skew = skew_res[1] + skew_res[2],
      mean_C     = C_res[1],
      sd_C       = C_res[2],
      lower_C    = C_res[1]    - C_res[2],
      upper_C    = C_res[1]    + C_res[2]
    )
  }))
}))

save.image("Masterarbeit.RData")
