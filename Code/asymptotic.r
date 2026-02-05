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

      H <- 12  # oder adaptiv: bis r^H < 1e-8
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
Sigma_Star(m, p, r, pi)
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

  Var_IOV <- (16 / (n * m^2)) * variance_sum

  IOV_real <- (4 / m) * sum(marginal_cdf * (1 - marginal_cdf))

  Sigma_trace <- sum(diag(Sigma))

  Expectation_IOV <- IOV_real - (1 / n) * (4 / m) * Sigma_trace

  return(c(Expectation_IOV, sqrt(Var_IOV)))
}

Skew_asymptotic <- function(Sigma, marginal_cdf, n,  m){
  skew_real <- (2/m) *  sum(marginal_cdf - 1)
  skew_expectation <- skew_real
  skew_variance <- (4 / n *m^2)* sum(Sigma)
  return(c(skew_expectation, skew_variance))
}



scenarios <- expand.grid(
  n = c(50, 100, 250, 500, 1000),
  m = c(3, 10),
  p = c(0.20, 0.45),
  r = c(0, 0.35, 0.50),
  pi = c(1, 0.95, 0.9, 0.85, 0.8, 0.75, 0.5)
)
scenarios <- scenarios[
  (scenarios$m == 3 & scenarios$p == 0.20 & scenarios$r %in% c(0, 0.35)) |
  (scenarios$m == 10 & scenarios$p == 0.45 & scenarios$r %in% c(0, 0.50)),
]
# Index reset
rownames(scenarios) <- NULL


target_scenarios <-  scenarios[
  scenarios$m == 3 & 
  scenarios$p == 0.20 & 
  scenarios$r == 0.35 &
  scenarios$pi %in% c(1, 0.75) &
  scenarios$n %in% c(50, 100, 250, 500, 1000),
]

IOV_asymp_results <- apply(target_scenarios, 1, function(row){
  marginal <- pbinom(0:(row["m"]-1), row["m"], row["p"])
  Sigma <- Sigma_Star(row["m"],row["p"], row["r"], row["pi"])
  return(IOV_asymptotic(Sigma, marginal, row["n"], row["m"]))
})

Skew_asymp_results <- apply(target_scenarios, 1, function(row){
  marginal <- pbinom(0:(row["m"]-1), row["m"], row["p"])
  Sigma <- Sigma_Star(row["m"],row["p"], row["r"], row["pi"])
  return(Skew_asymptotic(Sigma, marginal, row["n"], row["m"]))
})

Skew_asymp_results
# 1. Erstelle Dataframe für asymptotische Ergebnisse

asymp_df <- data.frame(
  n = sim_scenarios$n,
  pi = sim_scenarios$pi,
  m = sim_scenarios$m,
  p = sim_scenarios$p,
  r = sim_scenarios$r,
  type = "Asymptotic",
  mean = as.numeric(IOV_asymp_results[1, ]),
  sd = as.numeric(IOV_asymp_results[2, ]),
  lower = as.numeric(IOV_asymp_results[1, ]) - as.numeric(IOV_asymp_results[2, ]),
  upper = as.numeric(IOV_asymp_results[1, ]) + as.numeric(IOV_asymp_results[2, ])
)
asymp_df <- data.frame(
  n = target_scenarios$n,
  pi = target_scenarios$pi,
  m = target_scenarios$m,
  p = target_scenarios$p,
  r = target_scenarios$r,
  type = "Asymptotic",
  mean_IOV = as.numeric(IOV_asymp_results[1, ]),
  sd = as.numeric(IOV_asymp_results[2, ]),
  lower_IOV = as.numeric(IOV_asymp_results[1, ]) - as.numeric(IOV_asymp_results[2, ]),
  upper_IOV = as.numeric(IOV_asymp_results[1, ]) + as.numeric(IOV_asymp_results[2, ]),
  mean_Skew = as.numeric(Skew_asymp_results[1, ]),
  sd_Skew = as.numeric(Skew_asymp_results[2, ]),
  lower_Skew = as.numeric(Skew_asymp_results[1, ]) - as.numeric(Skew_asymp_results[2, ]),
  upper_Skew = as.numeric(Skew_asymp_results[1, ]) + as.numeric(Skew_asymp_results[2, ])
)