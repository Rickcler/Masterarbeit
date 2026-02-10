library(dplyr)
library(plyr)
library(tidyr)
library(ggplot2)



m <- 3
p <- 0.2
r <- 0.35
pi <- 1

f <- pbinom(0:(m-1), m, p)

# Test mit Ihren Parametern
lag_h_conditional <- function(m, p, r, h=1) {
  beta <- p * (1 - r^h)
  alpha <- beta + r^h
  res <- matrix(0, nrow = m + 1, ncol = m + 1)
  
  for (i in 0:m) {
    for (j in 0:m) {
      smaller <- min(i, j)
      bigger <- max(0, i + j - m)
      
      res[i + 1, j + 1] <- sum(sapply(bigger:smaller, function(n) {
        choose(i, n) * choose(m - i, j - n) * 
          alpha^n * (1 - alpha)^(i - n) * 
          beta^(j - n) * (1 - beta)^(m - i - j + n)
      }))
    }
  }
  
  return(res)  
}


lag_h_joint_pmf <- function(m, p, r, h = 1) {
  lag_h_conditional(m, p, r, h = h) * dbinom(0:m, m, p) 
}


#-- Get CDF from PMF matrix
pmf_to_cdf <- function(pmf) {
  cdf <- apply(pmf, 2, cumsum)
  cdf <- t(apply(cdf, 1, cumsum))
  return(cdf)
}


lag_h_joint_cdf <- function(m, p, r, h = 1) {
  return(pmf_to_cdf(lag_h_joint_pmf(m,p, r, h=h)))
}

lag_h_joint_cdf(m,p,r, h = 1000)

Sigma_Star <- function(m, p, r, pi) {

  # marginale CDF auf 0,...,m-1
  f <- pbinom(0:(m-1), m, p)


  Sigma <- matrix(0, nrow = m, ncol = m)

  for (i in 0:(m-1)) {
    for (j in 0:(m-1)) {

      smaller <- min(i, j)

      # i.i.d.-Teil (Missingness)
      iid_part <- (1 / pi) * (f[smaller + 1] - f[i + 1] * f[j + 1])

      H <- 50  # oder adaptiv: bis r^H < 1e-8
      lag_sum <- 0

      for (h in 1:H) {
        cdf_lag_h <- lag_h_joint_cdf(m,p, r, h)

        lag_sum <- lag_sum + cdf_lag_h[i+1, j+1] + cdf_lag_h[j+1, i+1] - 2 * f[i+1] * f[j+1] 
      }

      Sigma[i + 1, j + 1] <- iid_part + lag_sum
    }
  }

  return(Sigma)
}

Sigma_Star(10, 0.45, 0.5, 0.75)

IOV_asymptotic <- function(Sigma, marginal_cdf, n, m) {

  variance_sum <- 0

  for (i in 1:m) {
    for (j in 1:m) {
      variance_sum <- variance_sum +
        (1 - 2 * marginal_cdf[i]) *
        (1 - 2 * marginal_cdf[j]) *
        Sigma[i, j]
    }
  }


  Var_IOV <- (1/n)*(16/m^2) * variance_sum

  IOV_real <- (4 / m) * sum(marginal_cdf * (1 - marginal_cdf))

  Sigma_trace <- sum(diag(Sigma))

  Expectation_IOV <- IOV_real - (1 / n) * (4 / m) * Sigma_trace

  return(c(Expectation_IOV, sqrt(Var_IOV)))
}

Skew_asymptotic <- function(Sigma, marginal_cdf, n,  m){
  skew_real <- (2/m) *  sum(marginal_cdf - 1)
  skew_expectation <- skew_real
  skew_variance <- (1/n)*(4/m^2)* sum(Sigma)
  return(c(skew_expectation, sqrt(skew_variance)))
}



scenarios <- expand.grid(
  n = c(50, 100, 250, 500, 1000),
  m = c(3, 10),
  p = c(0.20, 0.45),
  r = c(0, 0.35, 0.50),
  pi = c(1, 0.75)
)

scenarios <- scenarios[
  (scenarios$m == 3 & scenarios$p == 0.20 & scenarios$r  == 0.35) | 
  (scenarios$m == 10 & scenarios$p == 0.45 & scenarios$r  == 0.5),]

# Index reset
rownames(scenarios) <- NULL
scenarios

unique_coeffs <- data.frame(rbind(c(m= 3, p = 0.2, r = 0.35, pi = 1),
                       c(m= 10, p = 0.45, r = 0.5, pi = 1),
                       c(m= 3, p = 0.2, r = 0.35, pi = 0.75),
                       c(m= 10, p = 0.45, r = 0.5, pi = 0.75)))
unique_n <- c(50, 100, 250, 500, 1000)



'Asymp_results <- apply(unique_coeffs,1, function(params){
  marginal <- pbinom(0:(params["m"]-1), params["m"], params["p"])
  Sigma <- Sigma_Star(params["m"],params["p"], params["r"], params["pi"])
  VEC <- sapply(unique_n, function(row){
    return(c(IOV_asymptotic(Sigma, marginal, row, params["m"]), Skew_asymptotic(Sigma, marginal, row, params["m"])))
  })
  return(VEC)
})
'
# Create a clean data frame for results
asymp_df <- data.frame()

for (scenario_idx in 1:nrow(unique_coeffs)) {
  params <- unique_coeffs[scenario_idx, ]
  m_val <- as.numeric(params["m"])
  p_val <- as.numeric(params["p"])
  r_val <- as.numeric(params["r"])
  pi_val <- as.numeric(params["pi"])
  
  marginal <- pbinom(0:(m_val-1), m_val, p_val)
  Sigma <- Sigma_Star(m_val, p_val, r_val, pi_val)
  
  for (n_val in unique_n) {
    iov_result <- IOV_asymptotic(Sigma, marginal, n_val, m_val)
    skew_result <- Skew_asymptotic(Sigma, marginal, n_val, m_val)
    
    # Create a row for this combination
    row <- data.frame(
      n = n_val,
      pi = pi_val,
      m = m_val,
      p = p_val,
      r = r_val,
      type = "Asymptotic",
      mean_IOV = iov_result[1],
      sd_IOV = iov_result[2],
      lower_IOV = iov_result[1] - iov_result[2],
      upper_IOV = iov_result[1] + iov_result[2],
      mean_Skew = skew_result[1],
      sd_Skew = skew_result[2],
      lower_Skew = skew_result[1] - skew_result[2],
      upper_Skew = skew_result[1] + skew_result[2]
    )
    
    asymp_df <- rbind(asymp_df, row)
  }
}

save.image("Masterarbeit.RData")

