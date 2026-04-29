# ==============================================================================
# Consistency of kappa_ord test under H_A
# ==============================================================================

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

# True kappa unter H_A (aus geschlossener Form)
true_kappa_HA <- function(m, p, r, h = 1) {
  f    <- pbinom(0:(m - 1), size = m, prob = p)
  cdf  <- lag_h_joint_cdf(m, p, r, h = h)
  f_ii <- diag(cdf)[1:m]  # f_{ii}(h) fuer i = 0,...,m-1
  sum(f_ii - f^2) / sum(f * (1 - f))
}

# Szenarien: verschiedene n, festes r > 0
n_grid  <- c(50, 100, 250, 500, 1000, 2000)
m_val   <- 3
p_val   <- 0.20
r_val   <- 0.35
pi_val  <- 0.75
h_val   <- 1

true_kappa <- true_kappa_HA(m_val, p_val, r_val, h = h_val)

set.seed(42)
kappa_HA_list <- lapply(n_grid, function(n) {
  vals <- simulation_kappa_HA(n, m_val, p_val, r_val, pi_val,
                               h = h_val, n_reps = 1000)
  data.frame(
    n         = n,
    kappa_hat = vals
  )
})

kappa_HA_df <- do.call(rbind, kappa_HA_list) %>%
  mutate(n = factor(n, levels = n_grid))

# Zusammenfassung: Mittelwert und SD pro n
kappa_summary <- kappa_HA_df %>%
  group_by(n) %>%
  summarise(
    mean_kappa = mean(kappa_hat, na.rm = TRUE),
    sd_kappa   = sd(kappa_hat,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(n_num = as.numeric(as.character(n)))

# ==============================================================================
# Plot 1: Mittlerer kappa_hat als Funktion von n
# (Zeigt Konsistenz: kappa_hat -> true kappa)
# ==============================================================================

ggplot(kappa_summary, aes(x = n_num, y = mean_kappa)) +
  geom_hline(yintercept = true_kappa,
             color = "steelblue", linetype = "dashed", linewidth = 0.8) +
  geom_hline(yintercept = 0,
             color = "grey60", linetype = "dotted", linewidth = 0.6) +
  geom_ribbon(aes(ymin = mean_kappa - sd_kappa,
                  ymax = mean_kappa + sd_kappa),
              fill = "grey80", alpha = 0.5) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  annotate("text", x = max(n_grid) * 0.6,
           y = true_kappa + 0.015,
           label = expression(paste("True ", kappa[ord](h))),
           color = "steelblue", size = 4) +
  annotate("text", x = max(n_grid) * 0.6,
           y = 0 + 0.015,
           label = expression(H[0]: kappa[ord](h) == 0),
           color = "grey50", size = 3.5) +
  scale_x_continuous(
    breaks = n_grid,
    name   = "n"
  ) +
  scale_y_continuous(name = expression(Mean(hat(kappa)[ord](h)))) +
  labs(
    title    = expression(
      paste("Consistency of ", hat(kappa)[ord](h), " under ", H[A])
    ),
    subtitle = bquote(
      m == .(m_val) ~ "," ~
      p == .(p_val) ~ "," ~
      r == .(r_val) ~ "," ~
      pi == .(pi_val) ~ "," ~
      h == .(h_val) ~
      "   Shaded: ±1 SD"
    )
  ) +
  theme_minimal() +
  theme(
    plot.title    = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "grey40")
  )

# ==============================================================================
# Plot 2: Teststatistik sqrt(n) * |kappa_hat| divergiert unter H_A
# (Direktes Konsistenzargument)
# ==============================================================================

kappa_stat_df <- kappa_HA_df %>%
  mutate(
    n_num    = as.numeric(as.character(n)),
    test_stat = sqrt(n_num) * abs(kappa_hat)
  ) %>%
  group_by(n, n_num) %>%
  summarise(
    mean_stat = mean(test_stat, na.rm = TRUE),
    sd_stat   = sd(test_stat,   na.rm = TRUE),
    .groups = "drop"
  )

# Theoretische Wachstumskurve: sqrt(n) * |true_kappa|
theory_stat <- data.frame(
  n_num     = seq(min(n_grid), max(n_grid), by = 1),
  true_stat = sqrt(seq(min(n_grid), max(n_grid), by = 1)) * abs(true_kappa)
)

ggplot(kappa_stat_df, aes(x = n_num, y = mean_stat)) +
  geom_line(data = theory_stat,
            aes(x = n_num, y = true_stat),
            color = "steelblue", linetype = "dashed", linewidth = 0.8) +
  geom_ribbon(aes(ymin = mean_stat - sd_stat,
                  ymax = mean_stat + sd_stat),
              fill = "grey80", alpha = 0.5) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  annotate("text",
           x = max(n_grid) * 0.55,
           y = sqrt(max(n_grid)) * abs(true_kappa) * 0.85,
           label = expression(sqrt(n) ~ "|" ~ kappa[ord](h) ~ "|"),
           color = "steelblue", size = 4) +
  scale_x_continuous(breaks = n_grid, name = "n") +
  scale_y_continuous(
    name = expression(Mean(sqrt(n) ~ "|" ~ hat(kappa)[ord](h) ~ "|"))
  ) +
  labs(
    title    = expression(
      paste("Divergence of test statistic ", sqrt(n), "|",
            hat(kappa)[ord](h), "| under ", H[A])
    ),
    subtitle = bquote(
      m == .(m_val) ~ "," ~
      p == .(p_val) ~ "," ~
      r == .(r_val) ~ "," ~
      pi == .(pi_val) ~ "," ~
      h == .(h_val) ~
      "   Dashed: " ~ sqrt(n) ~ "|" ~ kappa[ord](h) ~ "|"
    )
  ) +
  theme_minimal() +
  theme(
    plot.title    = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "grey40")
  )

ggsave("Graphs/kappa_consistency_HA_mean.pdf",
       width = 8, height = 5)
ggsave("Graphs/kappa_consistency_HA_stat.pdf",
       width = 8, height = 5)

# ==============================================================================
# Asymptotische Varianz von kappa_hat unter H_0 und i.i.d. O_t
# ==============================================================================

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

# ==============================================================================
# Konsistenz-Plot mit 95%-Konfidenzband unter H_0
# ==============================================================================

# Parameter
m_val  <- 3
p_val  <- 0.20
r_val  <- 0.35
pi_val <- 0.75
h_val  <- 1

n_grid <- c(50, 100, 250, 500, 1000, 2000)

# Wahres kappa unter H_A
true_kappa <- true_kappa_HA(m_val, p_val, r_val, h = h_val)

# Simulationen unter H_A
set.seed(42)
kappa_HA_list <- lapply(n_grid, function(n) {
  vals <- simulation_kappa_HA(n, m_val, p_val, r_val, pi_val,
                               h = h_val, n_reps = 1000)
  data.frame(n = n, kappa_hat = vals)
})

kappa_HA_df <- do.call(rbind, kappa_HA_list) %>%
  mutate(n = factor(n, levels = n_grid))

kappa_summary <- kappa_HA_df %>%
  group_by(n) %>%
  summarise(
    mean_kappa = mean(kappa_hat, na.rm = TRUE),
    sd_kappa   = sd(kappa_hat,   na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  mutate(n_num = as.numeric(as.character(n)))

# 95%-Konfidenzband unter H_0 für jedes n
ci_df <- data.frame(
  n_num = n_grid
) %>%
  rowwise() %>%
  mutate(
    var_kappa = asymp_var_kappa_H0(m_val, p_val, pi_val) / n_num,
    sd_H0     = sqrt(var_kappa),
    ci_lower  = -1.96 * sd_H0,
    ci_upper  =  1.96 * sd_H0
  ) %>%
  ungroup()

# ==============================================================================
# Plot
# ==============================================================================

ggplot() +
  # H_0-Nulllinie
  geom_hline(yintercept = 0,
             color = "grey50", linetype = "dotted", linewidth = 0.6) +
  # 95%-Konfidenzband unter H_0
  geom_ribbon(data = ci_df,
              aes(x = n_num, ymin = ci_lower, ymax = ci_upper),
              fill = "steelblue", alpha = 0.15) +
  geom_line(data = ci_df,
            aes(x = n_num, y = ci_upper),
            color = "steelblue", linetype = "dashed", linewidth = 0.7) +
  geom_line(data = ci_df,
            aes(x = n_num, y = ci_lower),
            color = "steelblue", linetype = "dashed", linewidth = 0.7) +
  annotate("text",
         x = max(n_grid) * 0.55, 
         y = ci_df$ci_upper[nrow(ci_df)] + 0.01,
         label = expression(paste("95% CI under ", H[0])),
         color = "steelblue", size = 3.5) +
  # Wahres kappa unter H_A
  geom_hline(yintercept = true_kappa,
             color = "#C0392B", linetype = "dashed", linewidth = 0.8) +
  annotate("text",
           x = max(n_grid) * 0.1,
           y = true_kappa + 0.015,
           label = expression(paste("True ", kappa[ord](h))),
           color = "#C0392B", size = 3.5) +
  # Simulierte Mittelwerte unter H_A
  geom_ribbon(data = kappa_summary,
              aes(x = n_num,
                  ymin = mean_kappa - sd_kappa,
                  ymax = mean_kappa + sd_kappa),
              fill = "grey40", alpha = 0.2) +
  geom_line(data = kappa_summary,
            aes(x = n_num, y = mean_kappa),
            linewidth = 0.9, color = "black") +
  geom_point(data = kappa_summary,
             aes(x = n_num, y = mean_kappa),
             size = 2.5, color = "black") +
  scale_x_continuous(breaks = n_grid, name = "n") +
  scale_y_continuous(
    name = expression(hat(kappa)[ord](h))
  ) +
  labs(
    title    = expression(
      paste("Consistency of ", hat(kappa)[ord](h),
            " under ", H[A], " — 95% CI under ", H[0])
    ),
    subtitle = bquote(
      m == .(m_val) ~ "," ~
      p == .(p_val) ~ "," ~
      r == .(r_val) ~ "," ~
      pi == .(pi_val) ~ "," ~
      h == .(h_val) ~
      "   Black: mean ± 1SD under" ~ H[A] ~
      "   Blue band: 95% CI under" ~ H[0]
    ),
    caption = "Rejection occurs when the black line exits the blue band"
  ) +
  theme_minimal() +
  theme(
    plot.title    = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    plot.caption  = element_text(size = 9, color = "grey50", face = "italic")
  )

ggsave("Graphs/kappa_consistency_CI.pdf", width = 9, height = 5.5)