library(dplyr)
library(plyr)
library(tidyr)
library(ggplot2)
library(pracma)


m <- 3
p <- 0.2
r <- 0.35
pi <- 1

f <- pbinom(0:(m-1), m, p)

# Test mit Ihren Parametern
m <- 3
p <- 0.2
r <- 0.35

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
  
  return(res)  # Rückgabewert hinzufügen
}


lag_h_joint_pmf <- function(m, p, r, h = 1) {
  lag_h_conditional(m, p, r, h = h) * dbinom(0:m, m, p) 
}

lag_h_joint_pmf(m, p, r, h = 1000)

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
Skew_asymptotic(Sigma_Star(m, p, r, pi), pbinom(0:2, 3, 0.2), 50, 3)


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
unique_coeffs <- rbind(c(m= 3, p = 0.2, r = 0.35, pi = 1),
                       c(m= 10, p = 0.45, r = 0.5, pi = 1),
                       c(m= 3, p = 0.2, r = 0.35, pi = 0.75),
                       c(m= 10, p = 0.45, r = 0.5, pi = 0.75))
unique_n <- c(50, 100, 250, 500, 1000)
Asymp_results <- apply(unique_coeffs, 1, function(params){
  marginal <- pbinom(0:(params["m"]-1), params["m"], params["p"])
  Sigma <- Sigma_Star(params["m"],params["p"], params["r"], params["pi"])
  return(sapply(unique_n,function(row){
    return(c(IOV_asymptotic(Sigma, marginal, row["n"], params["m"]), Skew_asymptotic(Sigma, marginal, row["n"], params["m"])))
  }))
})
Asymp_results <- sapply(scenarios, 1, function(row){
  marginal <- pbinom(0:(row["m"]-1), row["m"], row["p"])
  Sigma <- Sigma_Star(row["m"],row["p"], row["r"], row["pi"])
  return(c(IOV_asymptotic(Sigma, marginal, row["n"], row["m"]), Skew_asymptotic(Sigma, marginal, row["n"], row["m"])))
})

Skew_asymp_results <- apply(scenarios, 1, function(row){
  marginal <- pbinom(0:(row["m"]-1), row["m"], row["p"])
  Sigma <- Sigma_Star(row["m"],row["p"], row["r"], row["pi"])
  return(Skew_asymptotic(Sigma, marginal, row["n"], row["m"]))
})




asymp_df <- data.frame(
  n = scenarios$n,
  pi = scenarios$pi,
  m = scenarios$m,
  p = scenarios$p,
  r = scenarios$r,
  type = "Asymptotic",
  mean_IOV = as.numeric(IOV_asymp_results[1, ]),
  sd_IOV = as.numeric(IOV_asymp_results[2, ]),
  lower_IOV = as.numeric(IOV_asymp_results[1, ]) - as.numeric(IOV_asymp_results[2, ]),
  upper_IOV = as.numeric(IOV_asymp_results[1, ]) + as.numeric(IOV_asymp_results[2, ]),
  mean_Skew = as.numeric(Skew_asymp_results[1, ]),
  sd_Skew = as.numeric(Skew_asymp_results[2, ]),
  lower_Skew = as.numeric(Skew_asymp_results[1, ]) - as.numeric(Skew_asymp_results[2, ]),
  upper_Skew = as.numeric(Skew_asymp_results[1, ]) + as.numeric(Skew_asymp_results[2, ])
)