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
# Raw replications for CLT / bias illustration
sim_50   <- simulation_raw(50, 10, 0.3, 0.2, 0.75, 2000)
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

# Theoretische Bias-Kurve
theory_df_2 <- data.frame(
  n    = seq(2, 1000, by = 1),
  diff = sapply(seq(2, 1000, by = 1), function(n) -(4 / (10 * n)) * sum(diag(Sigma_raw)))
)


mean_df_2 <- sim_data_2 %>%
  dplyr::group_by(n) %>%
  dplyr::summarise(
    mean_diff        = mean(diff,        na.rm = TRUE),
    mean_scaled_bias = mean(scaled_bias, na.rm = TRUE),
    .groups = "drop"
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

ggplot(sim_data, aes(x = scaled_CLT, fill = factor(n))) +
  geom_density(alpha = 0.4) +
  stat_function(fun = dnorm,
      args = list(mean = mean(sim_data$scaled_CLT),
                  sd = sd(sim_data$scaled_CLT)),
      linetype = "dashed") +
  labs(fill = "n",
       title = "√n (IOV_hat - IOV)") +
  theme_minimal()

mean_df <- aggregate(diff ~ n, data = sim_data, mean)

ggplot(mean_df, aes(x = n, y = diff)) +
  geom_point() +
  geom_line() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal() +
  labs(y = "Mean(IOV_hat - IOV)")


ggplot(sim_data, aes(x = scaled_bias, fill = factor(n))) +
  geom_density(alpha = 0.4) +
  stat_function(fun = dnorm,
      args = list(mean = mean(sim_data$scaled_bias),
                  sd = sd(sim_data$scaled_bias)),
      linetype = "dashed") +
  labs(fill = "n",
       title = "n (IOV_hat - IOV)") +
  theme_minimal()
