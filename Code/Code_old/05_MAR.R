# ==============================================================================
# 05_MAR.R
# MCAR vs. MAR vs. MAR (inv.) Vergleich
# Voraussetzung: Masterarbeit.RData geladen (erzeugt von 02 + 03)
# ==============================================================================

source("00_setup.R")
load("Masterarbeit.RData")
library(astsa)
library(ggplot2)


library(ggplot2)
library(dplyr)

# Proband 1 extrahieren, NA-Zeilen behalten (sichtbar machen)

df_plot <- sleep1[[1]] %>%
  mutate(
    # Gruppe bricht bei NA auf – jede zusammenhängende Sequenz ohne NA
    # bekommt eine eigene Gruppen-ID
    gruppe = cumsum(is.na(state) | c(FALSE, is.na(head(state, -1))))
  ) %>%
  filter(!is.na(state))

sleep_plot <- ggplot(df_plot, aes(x = min, y = state, group = gruppe)) +
  geom_line(linewidth = 0.4, color = "steelblue") +
  geom_point(size = 5, color = "steelblue") +
  scale_y_continuous(
    breaks = 1:6,
    labels = c("Deep", "Light 2", "Light 1", "REM", "Wake", "Movement")
  ) +
  labs(
    x       = "Time (minutes)",
    y       = "Sleep state",
    caption = ""
  ) +
  theme_minimal()
print(sleep_plot)
ggsave("Graphs/sleep_example.png", sleep_plot, width = 8, height = 4)
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
print(plt)
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
  r = c(0.00, 0.15, 0.35, 0.6),
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


library(ggplot2)
library(dplyr)
library(scales)

rej_df_mnar <- rej_df_mnar %>%
  mutate(
    r_fac = factor(
      paste0("r = ", r),
      levels = c("r = 0", "r = 0.15", "r = 0.35", "r = 0.6")
    ),
    mech_fac = factor(
      mechanism,
      levels = c("increasing", "decreasing"),
      labels = c(
        "Increasing",
        "Decreasing"
      )
    )
  )

mech_fac = factor(
  mechanism,
  levels = c("increasing", "decreasing"),
  labels = c(
    "Increasing",
    "Decreasing"
  )
)

ggplot(
  rej_df_mnar,
  aes(
    x = n,
    y = rate,
    color = r_fac,
    linetype = mech_fac,
    group = interaction(r_fac, mech_fac)
  )
) +

  # Nominalniveau
  geom_hline(
    yintercept = 0.05,
    color = "grey50",
    linetype = "dotted",
    linewidth = 0.7
  ) +

  annotate(
    "text",
    x = max(rej_df_mnar$n) * 0.8,
    y = 0.08,
    label = "α = 5%",
    color = "grey50",
    size = 4
  ) +

  # Linien + Punkte
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.8) +

  # Farben nach r
  scale_color_manual(
    name = "Serial dependence",
    values = c(
      "r = 0"    = "black",
      "r = 0.15" = "#5B8DB8",
      "r = 0.35" = "#2C3E6B",
      "r = 0.6"  = "#C0392B"
    )
  ) +

  # Linientyp nach MNAR-Mechanismus
  scale_linetype_manual(
  name = "MNAR mechanism",
  values = c(
    "Increasing" = "solid",
    "Decreasing" = "dashed"
  ),
  labels = c(
    expression(P(O[t] == 1) %prop% X[t]),
    expression(P(O[t] == 1) %prop% -X[t])
  )
) +

  scale_x_continuous(
    trans  = "log10",
    breaks = c(50, 100, 250, 500, 1000, 2000),
    name   = "n (log scale)"
  ) +

  scale_y_continuous(
    limits = c(0, 1.05),
    breaks = seq(0, 1, by = 0.25),
    labels = percent,
    name   = "Rejection rate"
  ) +

  labs(
    title = expression(
      paste(
        "Rejection rate under MNAR missingness"
      )
    ),

    subtitle =
      "m = 3 , p = 0.2 , h = 1    π = 0.75"
  ) +
  guides(
  color    = guide_legend(order = 1, keywidth = unit(1, "cm")),
  linetype = guide_legend(order = 2, keywidth = unit(1, "cm"))
  ) +

  theme_minimal() +

  theme(
    plot.title = element_text(
      size = 13,
      face = "bold"
    ),
    axis.text.x        = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 12, color = "gray20"),
    axis.text.y        = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 12, color = "gray20"),
    plot.subtitle = element_text(
      size = 9,
      color = "grey40"
    ),

    legend.position = "right",

    legend.title = element_text(
      size = 8,
      face = "bold"
    ),

    legend.text = element_text(size = 8),

    panel.grid.minor = element_blank(),

    axis.ticks.length = unit(2.5, "mm")
  )

  
ggsave("Graphs/MNAR_kappa_rejection_rate.png", width = 5.5, height = 8)

# ------------------------------------------------------------------------------
# Konfidenzintervalle 
# ------------------------------------------------------------------------------
set.seed(42)

kappa_inc_list <- lapply(n_grid, function(n) {
  vals <- simulation_kappa_HA_MNAR(
    n, m_val, p_val, r_val, pi_val,
    h = h_val,
    n_reps = 1000,
    mechanism = "increasing"
  )
  data.frame(n = n, kappa_hat = vals, mech = "Increasing")
})

kappa_dec_list <- lapply(n_grid, function(n) {
  vals <- simulation_kappa_HA_MNAR(
    n, m_val, p_val, r_val, pi_val,
    h = h_val,
    n_reps = 1000,
    mechanism = "decreasing"
  )
  data.frame(n = n, kappa_hat = vals, mech = "Decreasing")
})

kappa_MNAR_df <- do.call(rbind, c(kappa_inc_list, kappa_dec_list)) %>%
  mutate(n = factor(n, levels = n_grid))

kappa_MNAR_summary <- kappa_MNAR_df %>%
  group_by(n, mech) %>%
  summarise(
    mean_kappa = mean(kappa_hat),
    sd_kappa   = sd(kappa_hat),
    .groups = "drop"
  ) %>%
  mutate(n_num = as.numeric(as.character(n)))

ci_df <- data.frame(n_num = n_grid) %>%
  rowwise() %>%
  mutate(
    sd_H0    = sqrt(asymp_var_kappa_H0(m_val, p_val, pi_val) / n_num),
    ci_lower = -1.96 * sd_H0,
    ci_upper =  1.96 * sd_H0
  ) %>%
  ungroup()
ggplot() +

  # H0 baseline
  geom_hline(yintercept = 0,
             color = "grey50",
             linetype = "dotted") +
  annotate("text",
         x = max(n_grid) * 0.25, 
         y = ci_df$ci_upper[nrow(ci_df)] - 0.01,
         label = expression(paste("95% CI under ", H[0])),
         color = "grey50", size = 5) +
  # H0 CI ribbon
  geom_ribbon(
    data = ci_df,
    aes(x = n_num, ymin = ci_lower, ymax = ci_upper),
    fill = "Grey20",
    alpha = 0.15
  ) +

  # MNAR ribbons (two scenarios)
  geom_ribbon(
    data = kappa_MNAR_summary,
    aes(
      x = n_num,
      ymin = mean_kappa - 1.96 * sd_kappa,
      ymax = mean_kappa + 1.96 * sd_kappa,
      fill = mech
    ),
    alpha = 0.15
  ) +

  # MNAR mean lines
  geom_line(
    data = kappa_MNAR_summary,
    aes(x = n_num, y = mean_kappa, color = mech),
    linewidth = 0.9
  ) +

  geom_point(
    data = kappa_MNAR_summary,
    aes(x = n_num, y = mean_kappa, color = mech),
    size = 2
  ) +

  # H0 bounds lines
  geom_line(
    data = ci_df,
    aes(x = n_num, y = ci_upper),
    color = "steelblue",
    linetype = "dashed"
  ) +

  geom_line(
    data = ci_df,
    aes(x = n_num, y = ci_lower),
    color = "steelblue",
    linetype = "dashed"
  ) +

  geom_hline(
    yintercept = true_kappa,
    color = "black",
    linetype = "dashed"
  ) +
  annotate("text",
           x = max(n_grid) * 0.15,
           y = true_kappa - 0.02,
           label = expression(paste("True ", kappa[ord](h))),
           color = "black", size = 5)  +
  scale_x_continuous(
    trans = "log10",
    breaks = n_grid,
    name = "n"
  ) +

  scale_y_continuous(
    name = expression(hat(kappa)[ord](1))
  ) +

scale_fill_manual(
  name   = "MNAR mechanism",
  values = c(
    "Increasing" = "steelblue",
    "Decreasing" = "#C0392B"
  ),
  labels = c(
    expression(P(O[t] == 1) %prop% X[t]),
    expression(P(O[t] == 1) %prop% -X[t])
  )
) +
guides(
  color = "none",    # entfernt die color-Legende
  fill  = guide_legend(title = "MNAR mechanism")  # behält nur fill
) +

  theme_minimal() +
  theme(
    legend.position = "right"
  )

ggsave("Graphs/MNAR_Kappa_H_A.png", width = 5.5, height = 8)
# ------------------------------------------------------------------------------
# Plot 2: Mittlerer Bias – MCAR vs. MAR
# ------------------------------------------------------------------------------







Sigma_mcar <- Sigma_Star(m = 3, p = 0.20, r = 0.35, pi = 0.75, r_pi = 0)

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
