install.packages("tidyr")
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
#' @param n_reps Number of replications
#' @return List with two matrices (summary, cdf), each with rows mean and sd
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

#' Simulation: returns raw replication-level IOV and Skew estimates
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

  results$n  <- n
  results$m  <- m
  results$p  <- p
  results$r  <- r
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
sim_50   <- rbind(numeric(),sim_50)
sim_1000 <- simulation_raw(1000, 10, 0.3, 0.2, 0.75, 2000)
sim_5000 <- simulation_raw(5000, 10, 0.3, 0.2, 0.75, 2000)


sim_data          <- rbind(sim_50, sim_1000, sim_5000)
true_val          <- true_IOV(m = 10, p = 0.3)
sim_data$true_IOV <- true_val
sim_data$diff     <- sim_data$IOV - true_val
sim_data$scaled_CLT  <- sqrt(sim_data$n) * sim_data$diff
sim_data$scaled_bias <- sim_data$n       * sim_data$diff

sim_data_2 <- numeric()
for  (n in c(2, 3, 5, 10, 15, 25, 37 , 50,75, 100, 150, 200, 375, 500, 750, 1000)) {
  sim_data_2 <- rbind(sim_data_2, simulation_raw(n, 10, 0.3, 0.2, 0.75, 1000))
  }
true_val          <- true_IOV(m = 10, p = 0.3)
sim_data_2$true_IOV <- true_val
sim_data_2$diff     <- sim_data_2$IOV - true_val
sim_data_2$scaled_CLT  <- sqrt(sim_data_2$n) * sim_data_2$diff
sim_data_2$scaled_bias <- sim_data_2$n       * sim_data_2$diff


mean_df_2 <- sim_data_2 %>%
  dplyr::group_by(n) %>%
  dplyr::summarise(
    mean_diff        = mean(diff,        na.rm = TRUE),
    mean_scaled_bias = mean(scaled_bias, na.rm = TRUE),
    .groups = "drop"
  )

Sigma_raw <- Sigma_Star(m = 10, p = 0.3, r = 0.2, pi = 0.75)
# Theoretische Bias-Kurve
theory_df_2 <- data.frame(
  n    = seq(2, 1000, by = 1),
  diff = sapply(seq(2, 1000, by = 1), function(n) -(4 / (10 * n)) * sum(diag(Sigma_raw)))
)

# Plot 1: Mittlerer Bias
ggplot(mean_df_2, aes(x = n, y = mean_diff)) +
  geom_line(data = theory_df_2, aes(x = n, y = diff),
            color = "steelblue", linetype = "dashed", linewidth = 0.8) +
  geom_line() +
  geom_point(size = 2) +
  labs(
    title    = "Mean bias as a function of n",
    subtitle = "Dashed line: theoretical bias  -4/(mn) · tr(Σ)\nBinAR(1): m = 10, p = 0.3, r = 0.2, π = 0.75",
    x = "n",
    y = expression(Mean(hat(IOV) - IOV))
  ) +
  theme_minimal()

# Plot 2: Mit n normalisierter mittlerer Bias
ggplot(mean_df_2, aes(x = n, y = mean_scaled_bias)) +
  geom_hline(
    yintercept = -(4 / 10) * sum(diag(Sigma_raw)),
    color = "steelblue", linetype = "dashed", linewidth = 0.8
  ) +
  geom_line() +
  geom_point(size = 2) +
  labs(
    title    = "n-scaled mean bias as a function of n",
    subtitle = "Dashed line: theoretical limit  -4/m · tr(Σ)\nBinAR(1): m = 10, p = 0.3, r = 0.2, π = 0.75",
    x = "n",
    y = expression(n ~ Mean(hat(IOV) - IOV))
  ) +
  theme_minimal()

# ==============================================================================
# Extract Results
# ==============================================================================

# Summary statistics per scenario (mean and sd rows)
sim_df <- do.call(rbind, lapply(seq_along(results), function(i) {
  s <- results[[i]]$summary
  data.frame(
    scenarios[i, ],
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
true_cdf <- pbinom(0:2, size = 3, prob = 0.20)
true_df  <- data.frame(
  category = paste0("f_", 0:2),
  true_cdf = true_cdf
)

# Geschätzte CDFs umformen
est_df <- cdf_df %>%
  pivot_longer(cols = starts_with("f_"), names_to = "category", values_to = "value") %>%
  pivot_wider(names_from = statistic, values_from = value) %>%
  mutate(n = factor(n, levels = sort(unique(n))))

# Asymptotische Standardabweichungen berechnen
Sigma <- Sigma_Star(m = 3, p = 0.20, r = 0.35, pi = 0.75)

asymp_df <- do.call(rbind, lapply(unique(est_df$n), function(ni) {
  ni_num <- as.numeric(as.character(ni))
  data.frame(
    category = paste0("f_", 0:2),
    n        = factor(ni, levels = levels(est_df$n)),
    asymp_sd = sqrt(diag(Sigma) / ni_num)
  )
}))

# Wahre CDF als Referenz hinzufügen
asymp_df$true_cdf <- true_cdf

plot_df <- bind_rows(
  est_df %>% 
    mutate(type = "Simulated", ymin = mean - sd, ymax = mean + sd),
  asymp_df %>% 
    rename(mean = true_cdf) %>%
    mutate(type = "Asymptotic", ymin = mean - asymp_sd, ymax = mean + asymp_sd)
) %>%
  mutate(type = factor(type, levels = c("Simulated", "Asymptotic")))

ggplot(plot_df, aes(x = category, color = n, linetype = type)) +
  geom_bar(
    data = true_df,
    aes(x = category, y = true_cdf),
    inherit.aes = FALSE,
    stat = "identity", fill = "grey85", color = "grey60", width = 0.7
  ) +
  geom_linerange(
    aes(ymin = ymin, ymax = ymax),
    position = position_dodge(width = 0.7),
    linewidth = 0.8
  ) +
  geom_point(
    aes(y = mean),
    position = position_dodge(width = 0.7),
    size = 2
  ) +
  scale_color_grey(name = "n", start = 0.7, end = 0.1) +
  scale_linetype_manual(
    name   = "Std. Deviation",
    values = c("Simulated" = "solid", "Asymptotic" = "dotted")
  ) +
  guides(
    linetype = guide_legend(
      override.aes = list(
        linewidth = 1,
        size      = 0,
        color     = "black"
      )
    )
  ) +
  scale_y_continuous(limits = c(0, 1), name = "Cumulative Probability") +
  scale_x_discrete(name = "CDF Component") +
  labs(
    title    = "Estimated vs. True CDF",
    subtitle = expression(m == 3 ~ "," ~ p == 0.20 ~ "," ~ r == 0.35 ~ "," ~ pi == 0.75)
  ) +
  theme_minimal()
save.image("Masterarbeit.RData")
