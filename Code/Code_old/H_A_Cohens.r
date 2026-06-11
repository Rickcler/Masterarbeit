# ==============================================================================
# Consistency of kappa_ord test under H_A
# ==============================================================================


save.image("Masterarbeit.RData")


# Szenarien: verschiedene n, festes r > 0
n_grid  <- c(50, 100, 250, 500, 1000, 2000)
m_val   <- 3
p_val   <- 0.20
r_val   <- 0.35
pi_val  <- 0.75
h_val   <- 1

true_kappa <- true_kappa_HA(m_val, p_val, r_val, h = h_val)

set.seed(SEED)
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
# Asymptotische Varianz von kappa_hat unter H_0 und i.i.d. O_t
# ==============================================================================


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
         x = max(n_grid) * 0.25, 
         y = ci_df$ci_upper[nrow(ci_df)] - 0.01,
         label = expression(paste("95% CI under ", H[0])),
         color = "steelblue", size = 5) +
  # Wahres kappa unter H_A
  geom_hline(yintercept = true_kappa,
             color = "#C0392B", linetype = "dashed", linewidth = 0.8) +
  annotate("text",
           x = max(n_grid) * 0.15,
           y = true_kappa + 0.05,
           label = expression(paste("True ", kappa[ord](h))),
           color = "#C0392B", size = 5) +
  # Simulierte Mittelwerte unter H_A
  geom_ribbon(data = kappa_summary,
              aes(x = n_num,
                  ymin = mean_kappa - (1.96 * sd_kappa),
                  ymax = mean_kappa + (1.96 * sd_kappa)),
              fill = "grey40", alpha = 0.2) +
  geom_line(data = kappa_summary,
            aes(x = n_num, y = mean_kappa),
            linewidth = 0.9, color = "black") +
  geom_point(data = kappa_summary,
             aes(x = n_num, y = mean_kappa),
             size = 2.5, color = "black") +
  scale_x_continuous(
    breaks = n_grid_rej,
    trans  = "log10",
    name   = ""
  ) +
  scale_y_continuous(
    name = expression(hat(kappa)[ord](1))
  ) +
  labs(
    title    = expression(
      paste("Consistency of ", hat(kappa)[ord](1),
            " under ", H[A], " — 95% CI under ", H[0])
    ),
    subtitle = bquote(
      m == .(m_val) ~ "," ~
      p == .(p_val) ~ "," ~
      r == .(r_val) ~ "," ~
      pi == .(pi_val) ~ "," ~
      h == .(h_val) ~
      "   Black: mean ± 1.96SD under" ~ H[A] ~
      "   Blue band: 95% CI under" ~ H[0]
    )
  ) +
  theme_minimal() +
  theme(
    axis.text.x        = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 12, color = "gray20"),
    axis.text.y        = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 14, color = "gray20"),
    plot.title    = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    plot.caption  = element_text(size = 9, color = "grey50", face = "italic"),
    axis.ticks.length = unit(2.5, "mm")
  )

ggsave("Graphs/Kappa_H_A.png", width = 5.5, height = 8)

# ==============================================================================
# Rejection Rate Plot — mehrere Szenarien
# Variiere r (Stärke der Abhängigkeit) und pi (Missingness)
# um unterschiedliche Konvergenzgeschwindigkeiten zu zeigen
# ==============================================================================

#' Compute rejection rate for kappa_ord test at alpha = 0.05
#' under H_A for a given scenario and n_grid
#' @param n       Time series length
#' @param m       Binomial parameter
#' @param p       Success probability
#' @param r       Thinning parameter (> 0 for H_A)
#' @param pi      Observation probability
#' @param h       Lag
#' @param alpha   Significance level
#' @param n_reps  Number of replications
rejection_rate <- function(n, m, p, r, pi, h = 1,
                            alpha = 0.05, n_reps = 1000) {
  # Kritischer Wert aus asymptotischer Varianz unter H_0
  sd_H0    <- sqrt(asymp_var_kappa_H0(m, p, pi) / n)
  crit_val <- qnorm(1 - alpha / 2) * sd_H0

  kappa_vals <- simulation_kappa_HA(n, m, p, r, pi,
                                     h = h, n_reps = n_reps)
  mean(abs(kappa_vals) > crit_val, na.rm = TRUE)
}

# ==============================================================================
# Szenarien: variiere r und pi
# r:  schwache (0.15), moderate (0.35), starke (0.60) Abhängigkeit
# pi: vollständige (1.0) und partielle (0.75) Beobachtung
# ==============================================================================

scenarios_rej <- expand.grid(
  r  = c(0.15, 0.35, 0.60),
  pi = c(1.00, 0.75)
) %>%
  mutate(
    label = paste0(
      "r = ", r, ", π = ", pi
    ),
    # Linientyp nach pi
    lty = ifelse(pi == 1.00, "solid", "dashed"),
    # Farbe nach r
    col = case_when(
      r == 0.15 ~ "#5B8DB8",
      r == 0.35 ~ "#2C3E6B",
      r == 0.60 ~ "#C0392B"
    )
  )

n_grid_rej <- c(50, 100, 250, 500, 1000, 2000)
m_val      <- 3
p_val      <- 0.20
h_val      <- 1

set.seed(42)
rej_list <- lapply(1:nrow(scenarios_rej), function(s) {
  sc <- scenarios_rej[s, ]
  message("Scenario: r = ", sc$r, ", pi = ", sc$pi)

  rates <- sapply(n_grid_rej, function(n) {
    rejection_rate(
      n      = n,
      m      = m_val,
      p      = p_val,
      r      = sc$r,
      pi     = sc$pi,
      h      = h_val,
      alpha  = 0.05,
      n_reps = 1000
    )
  })

  data.frame(
    n     = n_grid_rej,
    rate  = rates,
    r     = sc$r,
    pi    = sc$pi,
    label = sc$label,
    col   = sc$col,
    lty   = sc$lty
  )
})

rej_df <- do.call(rbind, rej_list) %>%
  mutate(
    label = factor(label, levels = scenarios_rej$label),
    r_fac = factor(paste0("r = ", r),
                   levels = c("r = 0.15", "r = 0.35", "r = 0.6")),
    pi_fac = factor(paste0("π = ", pi),
                    levels = c("π = 1", "π = 0.75"))
  )

# ==============================================================================
# Plot: Rejection Rate als Funktion von n — mehrere Szenarien
# ==============================================================================

ggplot(rej_df,
       aes(x     = n,
           y     = rate,
           color = r_fac,
           linetype = pi_fac,
           group = label)) +
  # 5%-Nominalniveau
  geom_hline(yintercept = 0.05,
             color = "grey50", linetype = "dotted", linewidth = 0.6) +
  annotate("text",
           x = max(n_grid_rej) * 0.85, y = 0.07,
           label = "α = 5%",
           color = "grey50", size = 3.5) +
  # 100%-Referenzlinie
  geom_hline(yintercept = 1.00,
             color = "grey80", linetype = "solid", linewidth = 0.4) +
  # Kurven
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  scale_color_manual(
    name   = "Serial dependence",
    values = c("r = 0.15" = "#5B8DB8",
               "r = 0.35" = "#2C3E6B",
               "r = 0.6" = "#C0392B"),
    #labels = c(
    #  expression(r == 0.15 ~ "(weak)"),
    #  expression(r == 0.35 ~ "(moderate)"),
    #  expression(r == 0.60 ~ "(strong)")
    #)
  ) +
  scale_linetype_manual(
    name   = "Observation probability",
    values = c("π = 1"    = "solid",
               "π = 0.75" = "dashed")
    )+
  scale_x_continuous(
    breaks = n_grid_rej,
    trans  = "log10",
    name   = ""
  ) +
  scale_y_continuous(
    limits = c(0, 1.05),
    breaks = seq(0, 1, by = 0.25),
    labels = scales::percent,
    name   = "Rejection rate"
  ) +
  labs(
    title    = expression(
      paste("Rejection rate of ", hat(kappa)[ord](h),
            " test under ", H[A], "  (α = 5%)")
    ),
    subtitle = bquote(
      m == .(m_val) ~ "," ~
      p == .(p_val) ~ "," ~
      h == .(h_val) ~
      "   Solid: π = 1   Dashed: π = 0.75   Dotted: nominal level"
    )
  ) +
  theme_minimal() +
  guides(
  color    = guide_legend(order = 1, keywidth = unit(1, "cm")),
  linetype = guide_legend(order = 2, keywidth = unit(1, "cm"))
  ) +
  theme(
    plot.title    = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9,  color = "grey40"),
    axis.text.x        = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 12, color = "gray20"),
    axis.text.y        = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 12, color = "gray20"),
    legend.position = "right",
    legend.title    = element_text(size = 8, face = "bold"),
    legend.text     = element_text(size = 7),
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(2.5, "mm")
  )

ggsave("Graphs/kappa_rejection_rate.png", width = 5.5, height = 8)