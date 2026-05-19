# ==============================================================================
# Rejection Rate unter H_0 — Level Control
# Nutzt simulation() direkt mit r = 0
# ==============================================================================

n_grid_H0  <- c(25, 37, 50, 75, 100,  200, 375, 500, 750, 1000, 2500, 5000)
alpha_test <- 0.05

# Szenarien A und B mit r = 0 (H_0: i.i.d.)
scenarios_H0 <- data.frame(
  m     = c(3,    10),
  p     = c(0.20, 0.45),
  pi    = c(0.75, 0.75),
  pi_h  = c(0,    0),

  label = c("Scenario A  (m = 3, p = 0.20)",
            "Scenario B  (m = 10, p = 0.45)")
)

#' Rejection rate from raw simulation results
#' @param kappa_vals  Vector of simulated kappa values
#' @param crit        Critical value
rejection_from_vals <- function(kappa_vals, crit) {
  mean(abs(kappa_vals) > crit, na.rm = TRUE)
}

set.seed(42)
rej_H0_list <- lapply(1:nrow(scenarios_H0), function(s) {
  sc <- scenarios_H0[s, ]

  do.call(rbind, lapply(n_grid_H0, function(n) {

    # Asymptotische Varianz unter H_0 für kritischen Wert
    asymp <- Cohens_asymptotic_iid(
      n,
      sc$pi,
      sc$m,
      pbinom(0:(sc$m - 1), sc$m, sc$p)
    )

    crit <- qnorm(1 - alpha_test / 2) * asymp["sd"]


    # Wir brauchen rohe Werte — simulation() gibt nur mean/sd zurück
    # Daher raw-Replikation inline
    kappa_raw <- replicate(5000, {
      iid_proc <- generate_iid(n, sc$m, sc$p)
      O_proc   <- generate_binar1(n, 1, sc$pi, sc$pi_h)
      obs_iid  <- iid_proc[O_proc == 1]
      CDF_iid  <- marginal_probs_e(sc$m, obs_iid, length(obs_iid))
      pi_est   <- mean(O_proc)
      biv      <- biv_probs_e(sc$m, iid_proc, O_proc, n, h = 1)
      denom    <- sum(CDF_iid * (1 - CDF_iid))
      if (denom > 0 && pi_est > 0) {
        kappa    <- sum(biv - CDF_iid^2) / denom
        kappa_bc <- kappa + 1 / (pi_est * n)
        c(kappa, kappa_bc)
      } else {
        c(NA, NA)
      }
    })

    data.frame(
      n          = n,
      scenario   = sc$label,
      m          = sc$m,
      p          = sc$p,
      pi         = sc$pi,
      rej_uncorr = rejection_from_vals(kappa_raw[1, ], crit),
      rej_bc     = rejection_from_vals(kappa_raw[2, ], crit)
    )
  }))
})

rej_H0_df <- do.call(rbind, rej_H0_list) %>%
  pivot_longer(
    cols      = c(rej_uncorr, rej_bc),
    names_to  = "estimator",
    values_to = "rejection_rate"
  ) %>%
  mutate(
    estimator = factor(
      estimator,
      levels = c("rej_uncorr", "rej_bc"),
      labels = c("Uncorrected", "Bias-corrected")
    ),
    scenario = factor(scenario, levels = scenarios_H0$label)
  )

# ==============================================================================
# Plot
# ==============================================================================

ggplot(rej_H0_df,
       aes(x        = n,
           y        = rejection_rate,
           color    = estimator,
           group    = estimator)) +
  # Nominalniveau
  geom_hline(yintercept = alpha_test,
             color     = "grey40",
             linetype  = "dashed",
             linewidth = 0.7) +
  annotate("text",
           x     = 40,
           y     = alpha_test - 0.008,
           label = "Nominal level (5%)",
           color = "grey40",
           size  = 3.5) +
  # Kurven
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  facet_wrap(~ scenario, ncol = 2) +
  scale_color_manual(
    name   = "Test statistic",
    values = c("Uncorrected"    = "#2C3E6B",
               "Bias-corrected" = "#C0392B")
  )  +
  scale_x_log10(name = "log(n)") +
  scale_y_continuous(
    limits = c(0, 0.1),
    breaks = seq(0, 0.20, by = 0.05),
    labels = scales::percent_format(accuracy = 1),
    name   = "Empirical rejection frequency"
  ) +
  labs(
    title    = expression(
      paste("Empirical rejection frequency under ",
            H[0], "  (", alpha, " = 5%)")
    )
  ) +
  theme_minimal() +
  theme(
    plot.title       = element_text(size = 13, face = "bold"),
    plot.subtitle    = element_text(size = 11,  color = "grey40"),
    plot.caption     = element_text(size = 9,  color = "grey50",
                                    face = "italic"),
    legend.position  = "bottom",
    legend.title     = element_text(size = 10, face = "bold"),
    legend.text      = element_text(size = 9),
    strip.text       = element_text(size = 11, face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  guides(
    color    = guide_legend(order = 1),
    linetype = guide_legend(order = 2)
  )

ggsave("Graphs/rejection_H0_level.png", width = 10, height = 5)