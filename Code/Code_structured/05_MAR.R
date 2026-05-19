# ==============================================================================
# 05_MAR.R
# MCAR vs. MAR vs. MAR (inv.) Vergleich
# Voraussetzung: Masterarbeit.RData geladen (erzeugt von 02 + 03)
# ==============================================================================

source("00_setup.R")
load("Masterarbeit.RData")

# ------------------------------------------------------------------------------
# Daten aufbereiten
# ------------------------------------------------------------------------------

# MCAR: direkt aus sim_df (m=3, p=0.20, r=0.35, pi=0.75, pi_h=0)
mcar_df <- sim_df %>%
  filter(m == 3, p == 0.20, r == 0.35, pi == 0.75, r_pi == 0) %>%
  mutate(mechanism = "MCAR") %>%
  select(n, mechanism, mean_IOV, sd_IOV, mean_Skew, sd_Skew, mean_C, sd_C)

# MAR: Szenarien für m=3, p=0.20, r=0.35 aus results_MAR
mar_idx <- which(
  scenarios_MAR$m == 3 &
  scenarios_MAR$p == 0.20 &
  scenarios_MAR$r == 0.35
)

mar_df <- do.call(rbind, lapply(seq_along(mar_idx), function(i) {
  idx <- mar_idx[i]
  s   <- results_MAR[[idx]]$summary
  data.frame(
    n         = scenarios_MAR$n[idx],
    mechanism = "MAR",
    mean_IOV  = s["mean", "IOV"],
    sd_IOV    = s["sd",   "IOV"],
    mean_Skew = s["mean", "Skew"],
    sd_Skew   = s["sd",   "Skew"],
    mean_C    = s["mean", "lag1_Cohen"],
    sd_C      = s["sd",   "lag1_Cohen"]
  )
}))

# MAR (inv.)
mar_inv_df <- do.call(rbind, lapply(seq_along(mar_idx), function(i) {
  idx <- mar_idx[i]
  s   <- results_MAR_inv[[idx]]$summary
  data.frame(
    n         = scenarios_MAR$n[idx],
    mechanism = "MAR (inv.)",
    mean_IOV  = s["mean", "IOV"],
    sd_IOV    = s["sd",   "IOV"],
    mean_Skew = s["mean", "Skew"],
    sd_Skew   = s["sd",   "Skew"],
    mean_C    = s["mean", "lag1_Cohen"],
    sd_C      = s["sd",   "lag1_Cohen"]
  )
}))

compare_df <- bind_rows(mcar_df, mar_df, mar_inv_df) %>%
  mutate(
    mechanism = factor(mechanism,
                       levels = c("MCAR", "MAR", "MAR (inv.)")),
    n = as.numeric(as.character(n))
  )

true_iov  <- true_IOV(m = 3, p = 0.20)
true_skew <- true_Skew(m = 3, p = 0.20)

compare_df$mean_bias        <- compare_df$mean_IOV - true_iov
compare_df$mean_scaled_bias <- compare_df$n * compare_df$mean_bias

# Manuelle Farbpalette (konsistent über alle Plots)
mech_colors    <- c("MCAR" = "black", "MAR" = "steelblue",
                    "MAR (inv.)" = "tomato")
mech_linetypes <- c("MCAR" = "solid", "MAR" = "dashed",
                    "MAR (inv.)" = "dotdash")

# ------------------------------------------------------------------------------
# Plot 1: Mittlerer IOV – MCAR vs. MAR
# ------------------------------------------------------------------------------

plt <- ggplot(compare_df,
       aes(x = n, y = mean_IOV, color = mechanism,
           linetype = mechanism, group = mechanism)) +
  geom_hline(yintercept = true_iov,
             color = "grey40", linetype = "dashed", linewidth = 0.7) +
  annotate("text", x = 55, y = true_iov + 0.005,
           label = "True IOV", color = "grey40", size = 3, hjust = 0) +
  geom_ribbon(aes(ymin = mean_IOV - sd_IOV,
                  ymax = mean_IOV + sd_IOV,
                  fill = mechanism),
              alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_color_manual(name   = "Missingness mechanism",
                     values = mech_colors) +
  scale_fill_manual(name    = "Missingness mechanism",
                    values  = mech_colors) +
  scale_linetype_manual(name   = "Missingness mechanism",
                        values = mech_linetypes) +
  scale_x_continuous(breaks = UNIQUE_N, trans = "log10") +
  labs(
    title    = "Mean estimated IOV under MCAR, MAR, and inverted MAR",
    subtitle = expression(
      m == 3 ~ "," ~ p == 0.20 ~ "," ~ r == 0.35 ~ "   Shaded: ±1 SD"
    ),
    x = "n (log scale)",
    y = expression(Mean(widehat(IOV)))
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave("Graphs/iov_mnar.png", width = 8, height = 5)
# ------------------------------------------------------------------------------
# Plot 2: Mittlerer Bias – MCAR vs. MAR
# ------------------------------------------------------------------------------

# Simulation Function

simulation_kappa_HA_MNAR <- function(
    n, m, p, r,
    pi_low = 0.5,
    pi_high = 0.9,
    h = 1,
    n_reps = 1000,
    mechanism = c("increasing", "decreasing")
) {

  mechanism <- match.arg(mechanism)

  kappa_vals <- numeric(n_reps)

  for(rep in 1:n_reps) {

    # abhängiger Prozess
    count_process <- generate_binar1(n, m, p, r)

    # MNAR Missingness
    O_t <- switch(
      mechanism,

      increasing =
        generate_O_MAR(
          count_process,
          m,
          pi_low,
          pi_high
        ),

      decreasing =
        generate_O_MAR_inv(
          count_process,
          m,
          pi_low,
          pi_high
        )
    )

    # geschätzte CDF
    CDF_hat <- marginal_probs_e(
      m,
      count_process[O_t == 1],
      sum(O_t)
    )

    # geschätzte biv probs
    biv_hat <- biv_probs_e(
      m,
      count_process,
      O_t,
      n,
      h = h
    )

    # Cohen-Kappa-artige Statistik
    num <- sum(biv_hat - CDF_hat^2)
    den <- sum(CDF_hat * (1 - CDF_hat))

    kappa_vals[rep] <- num / den
  }

  kappa_vals
}


# Rejection Rate
rejection_rate_MNAR <- function(
    n,
    m,
    p,
    r,
    pi_low = 0.5,
    pi_high = 0.9,
    h = 1,
    alpha = 0.05,
    n_reps = 1000,
    mechanism = c("increasing", "decreasing")
) {

  mechanism <- match.arg(mechanism)

  # asymptotischer kritischer Wert
  # (falsch unter MNAR -> genau das willst du zeigen)
  sd_H0 <- sqrt(asymp_var_kappa_H0(m, p, pi = 0.75) / n)

  crit_val <- qnorm(1 - alpha / 2) * sd_H0

  kappa_vals <- simulation_kappa_HA_MNAR(
    n         = n,
    m         = m,
    p         = p,
    r         = r,
    pi_low    = pi_low,
    pi_high   = pi_high,
    h         = h,
    n_reps    = n_reps,
    mechanism = mechanism
  )

  mean(abs(kappa_vals) > crit_val, na.rm = TRUE)
}

# Scenarios
scenarios_mnar <- expand.grid(
  r = c(0.00, 0.15, 0.35),
  mechanism = c("increasing", "decreasing"),
  stringsAsFactors = FALSE
)
# Simulation
n_grid_mnar <- c(50, 100, 250, 500, 1000, 2000)

rej_list_mnar <- lapply(1:nrow(scenarios_mnar), function(s) {

  sc <- scenarios_mnar[s, ]

  rates <- sapply(n_grid_mnar, function(n) {

    rejection_rate_MNAR(
      n          = n,
      m          = 3,
      p          = 0.20,
      r          = sc$r,
      pi_low     = 0.5,
      pi_high    = 0.9,
      h          = 1,
      alpha      = 0.05,
      n_reps     = 1000,
      mechanism  = sc$mechanism
    )
  })

  data.frame(
    n = n_grid_mnar,
    rate = rates,
    r = sc$r,
    mechanism = sc$mechanism
  )
})


rej_df_mnar <- do.call(rbind, rej_list_mnar)

# ------------------------------------------------------------------------------
# Plot 2: Mittlerer Bias – MCAR vs. MAR
# ------------------------------------------------------------------------------






Sigma_mcar <- Sigma_Star(m = 3, p = 0.20, r = 0.35, pi = 0.75, pi_h = 0)

theory_bias <- data.frame(
  n    = seq(50, 1000, by = 1),
  bias = sapply(seq(50, 1000, by = 1), function(ni)
    -(4 / (3 * ni)) * sum(diag(Sigma_mcar)))
)

ggplot(compare_df, aes(x = n, y = mean_bias,
                       color = mechanism, linetype = mechanism)) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.5) +
  geom_line(data = theory_bias, aes(x = n, y = bias),
            inherit.aes = FALSE,
            color = "black", linetype = "dotted", linewidth = 0.8) +
  annotate("text", x = 600, y = min(theory_bias$bias) - 0.002,
           label = "Theoretical bias (MCAR)",
           color = "black", size = 3) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_color_manual(name   = "Missingness mechanism",
                     values = mech_colors) +
  scale_linetype_manual(name   = "Missingness mechanism",
                        values = mech_linetypes) +
  scale_x_continuous(breaks = UNIQUE_N, trans = "log10") +
  labs(
    title    = "Mean bias of IOV estimator under MCAR vs. MAR",
    subtitle = expression(
      m == 3 ~ "," ~ p == 0.20 ~ "," ~ r == 0.35 ~
      "   Dotted: theoretical bias under MCAR"
    ),
    x = "n (log scale)",
    y = expression(Mean(widehat(IOV) - IOV))
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# ------------------------------------------------------------------------------
# Plot 3: Mittlere Skewness – MCAR vs. MAR
# ------------------------------------------------------------------------------

ggplot(compare_df,
       aes(x = n, y = mean_Skew, color = mechanism,
           linetype = mechanism, group = mechanism)) +
  geom_hline(yintercept = true_skew,
             color = "grey40", linetype = "dashed", linewidth = 0.7) +
  annotate("text", x = 55, y = true_skew + 0.005,
           label = "True Skew", color = "grey40", size = 3, hjust = 0) +
  geom_ribbon(aes(ymin = mean_Skew - sd_Skew,
                  ymax = mean_Skew + sd_Skew,
                  fill = mechanism),
              alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_color_manual(name   = "Missingness mechanism",
                     values = mech_colors) +
  scale_fill_manual(name    = "Missingness mechanism",
                    values  = mech_colors) +
  scale_linetype_manual(name   = "Missingness mechanism",
                        values = mech_linetypes) +
  scale_x_continuous(breaks = UNIQUE_N, trans = "log10") +
  labs(
    title    = "Mean estimated Skewness under MCAR, MAR, and inverted MAR",
    subtitle = expression(
      m == 3 ~ "," ~ p == 0.20 ~ "," ~ r == 0.35 ~ "   Shaded: ±1 SD"
    ),
    x = "n (log scale)",
    y = expression(Mean(widehat(Skew)))
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
