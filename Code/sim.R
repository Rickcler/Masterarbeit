library(dplyr)
library(plyr)
library(tidyr)



#' Generate BinAR(1) process
generate_binar1 <- function(n, m, p, r) {
  # n: Länge der Zeitreihe
  # m: Parameter der Binomialverteilung
  # p: Erwartungswert der stationären Verteilung = E[I_t]
  # r: Thinning-Parameter (muss in [0, 1] sein)
  alpha <- (p * (1-r)) + r
  beta <- p * (1-r)
  # Initialisiere
  I <- numeric(n)
  
  # Startwert aus stationärer Verteilung: Bin(m, p)
  I[1] <- rbinom(1, m, p)
  
  # Erzeuge den Prozess
  for (t in 2:n) {
    
    I[t] <- rbinom(1, I[t-1], alpha) + rbinom(1, m - I[t-1], beta)  
  }
  
  return(I)
}
#'' Generate Amplitude-modulating process
generate_O <- function(n, pi) {
  rbinom(n, 1, pi)
}
# Test with small example
set.seed(123)
test_counts <- generate_binar1(n = 10000, m = 3, p = 0.2, r = 0.35)

generate_iid <- function(p, m, n) {
  rbinom(n, m, p)
} 

biv_probs_e <- function(m, data, O, n, h = 1) {
  biv_probs <- numeric(m)
  valid <- (O[1:(n-h)] == 1) & (O[(h+1):n] == 1)
  for (i in 0:(m-1)) {
    biv_probs[i+1] <- sum(
      (data[1:(n-h)] <= i) *
      (data[(h+1):n] <= i) *
      valid
    ) / sum(valid)
  }
  return(biv_probs)
}


#-- Function to compute the empirical marginal cdf 
marginal_probs_e <- function(m, data, n) {
  marg_probs <- numeric(m)
  for (i in 0:(m-1)) {
    marg_probs[i+1] <- sum(data <= i)*(1/n)
  } 
  return(marg_probs)
}

#-- Function to compute the empirical lag(1) cdf
biv_probs_e <- function(m, data, n, h = 1) {
  biv_probs <- numeric(m)
  for (i in 0:(m-1)) {
    biv_probs[i+1] <- (1/(n-h)) * sum((data[1:(n-h)] <= i)*(data[(h+1):n] <= i))
  } 
  return(biv_probs)
}
biv_probs_e(3, generate_iid(0.2, 3, 1000), 1000)
pbinom(0:2, 3, 0.2)
#' Simulation function
simulation <- function(n, m, p, r, pi, n_reps = 100) {
  results <- matrix(NA, nrow = n_reps, ncol = 4, dimnames = list(NULL, c("IOV", "Skew", "lag1_Cohen", "lag2_Cohen")))
  for (rep in 1:n_reps) {
  # Generate BinAR(1) process
  count_process <- generate_binar1(n, m, p, r)
  Missing_process <- generate_O(n, pi) 
  observed_counts <- count_process[Missing_process ==1]
  CDF <- marginal_probs_e(m, observed_counts, length(observed_counts))
  
  iid_process <- generate_iid(p, m, n) # Für Cohens(k)
  observed_counts_iid <- iid_process[Missing_process == 1]
  CDF_iid <- marginal_probs_e(m, observed_counts_iid, length(observed_counts_iid))
  results[rep, 1] <- (4/m)*(sum(CDF *(1-CDF)))
  results[rep, 2] <- (2/m)*sum(CDF-1)
  results[rep, 3] <-   sum(biv_probs_e(m, iid_process, Missing_process, n, h = 1) - CDF_iid^2) /sum(CDF_iid*(1-CDF_iid))
  results[rep, 4] <- sum(biv_probs_e(m, iid_process, Missing_process, n, h = 2) - CDF_iid^2)/sum(CDF_iid*(1-CDF_iid))
  }
  return(rbind(colMeans(results, na.rm= T), apply(results, 2, sd)))
}

true_IOV <- function(m, p) {
  F <- pbinom(0:(m-1), m, p)
  return((4/m) * sum(F * (1 - F)))
}

simulation_raw <- function(n, m, p, r, pi, n_reps = 1000) {

  results <- data.frame(
    rep = 1:n_reps,
    IOV = NA,
    Skew = NA
  )

  for (rep in 1:n_reps) {

    count_process <- generate_binar1(n, m, p, r)
    Missing_process <- generate_O(n, pi)
    observed_counts <- count_process[Missing_process == 1]

    CDF <- marginal_probs_e(m, observed_counts, length(observed_counts))

    results$IOV[rep] <- (4/m) * sum(CDF * (1 - CDF))
    results$Skew[rep] <- (2/m) * sum(CDF - 1)
  }

  results$n <- n
  results$m <- m
  results$p <- p
  results$r <- r
  results$pi <- pi

  return(results)
}




# Alle Szenarien generieren
scenarios <- expand.grid(
  n = c(50, 100, 250, 500, 1000),
  m = c(3, 10),
  p = c(0.20, 0.45),
  r = c(0, 0.35, 0.50),
  pi = c(1, 0.75)
)
scenarios <- scenarios[
  (scenarios$m == 3 & scenarios$p == 0.20 & scenarios$r == 0.35) |
  (scenarios$m == 10 & scenarios$p == 0.45 & scenarios$r == 0.50),
]

# Index reset
rownames(scenarios) <- NULL

sim_50  <- simulation_raw(50, 10, 0.3, 0.2, 0.75, 2000)
sim_1000 <- simulation_raw(1000, 10, 0.3, 0.2, 0.75, 2000)
sim_5000 <- simulation_raw(5000, 10, 0.3, 0.2, 0.75, 2000)

sim_data <- rbind(sim_50, sim_1000, sim_5000)

true_val <- true_IOV(m = 10, p = 0.3)

sim_data$true_IOV <- true_val
sim_data$diff <- sim_data$IOV - true_val
sim_data$scaled_CLT <- sqrt(sim_data$n) * sim_data$diff
sim_data$scaled_bias <- sim_data$n * sim_data$diff


results <- apply(scenarios, 1, function(row) {
  simulation(
    n = row["n"], m = row["m"], p = row["p"], 
    r = row["r"], pi = row["pi"], n_reps = 1000
  )
})



save.image("Masterarbeit.RData")

IOV_sim <- results[1:2,]
Skew_sim <- results[3:4, ]
Co_1_sim <- results[5:6, ]
Co_2_sim <- results[7:8, ]


# Extrahiere die entsprechenden simulierten Ergebnisse
sim_IOV_means <- as.numeric(IOV_sim[1, ])  # Erste Zeile: Mittelwerte
sim_IOV_sds <- as.numeric(IOV_sim[2, ])    # Zweite Zeile: Standardabweichungen

sim_Skew_means <- as.numeric(Skew_sim[1, ])  # Erste Zeile: Mittelwerte
sim_Skew_sds <- as.numeric(Skew_sim[2, ])    # Zweite Zeile: Standardabweichungen

sim_Co_1_means <- as.numeric(Co_1_sim[1, ])  # Erste Zeile: Mittelwerte
sim_Co_1_sds <- as.numeric(Co_1_sim[2, ])    # Zweite Zeile: Standardabweichungen


sim_df <- data.frame(
  n = scenarios$n,
  pi = scenarios$pi,
  m = scenarios$m,
  p = scenarios$p,
  r = scenarios$r,
  type = "Simulated",
  mean_IOV = sim_IOV_means,
  sd_IOV = sim_IOV_sds,
  lower_IOV = sim_IOV_means - sim_IOV_sds,
  upper_IOV = sim_IOV_means + sim_IOV_sds,
  mean_Skew = sim_Skew_means,
  sd_Skew = sim_Skew_sds,
  lower_Skew = sim_Skew_means - sim_Skew_sds,
  upper_Skew = sim_Skew_means + sim_Skew_sds,
  mean_C = sim_Co_1_means,
  sd_C = sim_Co_1_sds,
  lower_C = sim_Co_1_means - sim_Co_1_sds,
  upper_C = sim_Co_1_means + sim_Co_1_sds
)