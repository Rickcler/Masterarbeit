# ==============================================================================
# Vergleich MCAR vs. MAR: Aufbereitung - korrigiert
# ==============================================================================

# MCAR direkt aus sim_df extrahieren (bereits sauber aufbereitet)
mcar_df <- sim_df %>%
  filter(m == 3, p == 0.20, r == 0.35, pi == 0.75) %>%
  mutate(mechanism = "MCAR") %>%
  select(n, mechanism, mean_IOV, sd_IOV, mean_Skew, sd_Skew, mean_C, sd_C)

# MAR Ergebnisse neu aufbereiten
mar_idx <- which(
  scenarios_MAR$m == 3 &
  scenarios_MAR$p == 0.20 &
  scenarios_MAR$r == 0.35
)

mar_df <- do.call(rbind, lapply(seq_along(mar_idx), function(i) {
  idx <- mar_idx[i]
  s   <- results_MAR[[i]]$summary
  data.frame(
    n          = scenarios_MAR$n[idx],
    mechanism  = "MAR",
    mean_IOV   = s[1, "IOV"],
    sd_IOV     = s[2, "IOV"],
    mean_Skew  = s[1, "Skew"],
    sd_Skew    = s[2, "Skew"],
    mean_C     = s[1, "lag1_Cohen"],
    sd_C       = s[2, "lag1_Cohen"]
  )
}))

mar_inv_idx <- which(
  scenarios_MAR$m == 3 &
  scenarios_MAR$p == 0.20 &
  scenarios_MAR$r == 0.35
)

mar_inv_df <- do.call(rbind, lapply(seq_along(mar_inv_idx), function(i) {
  idx <- mar_inv_idx[i]
  s   <- results_MAR_inv[[i]]$summary
  data.frame(
    n         = scenarios_MAR$n[idx],
    mechanism = "MAR (inv.)",
    mean_IOV  = s[1, "IOV"],
    sd_IOV    = s[2, "IOV"],
    mean_Skew = s[1, "Skew"],
    sd_Skew   = s[2, "Skew"],
    mean_C    = s[1, "lag1_Cohen"],
    sd_C      = s[2, "lag1_Cohen"]
  )
}))




# Alle drei zusammenführen
compare_df <- bind_rows(mcar_df, mar_df, mar_inv_df) %>%
  mutate(
    mechanism = factor(mechanism,
                       levels = c("MCAR", "MAR", "MAR (inv.)")),
    n = as.numeric(as.character(n))
  )

true_iov <- true_IOV(m = 3, p = 0.20)

# ==============================================================================
# Plot 1: Mittlerer IOV - MCAR vs. MAR
# ==============================================================================

ggplot(compare_df,
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
  scale_color_manual(
    name   = "Missingness mechanism",
    values = c("MCAR"       = "black",
               "MAR"        = "steelblue",
               "MAR (inv.)" = "tomato")
  ) +
  scale_fill_manual(
    name   = "Missingness mechanism",
    values = c("MCAR"       = "black",
               "MAR"        = "steelblue",
               "MAR (inv.)" = "tomato")
  ) +
  scale_linetype_manual(
    name   = "Missingness mechanism",
    values = c("MCAR"       = "solid",
               "MAR"        = "dashed",
               "MAR (inv.)" = "dotdash")
  ) +
  scale_x_continuous(
    breaks = c(50, 100, 250, 500, 1000),
    trans  = "log10"
  ) +
  labs(
    title    = "Mean estimated IOV under MCAR, MAR, and inverted MAR",
    subtitle = expression(
      m == 3 ~ "," ~ p == 0.20 ~ "," ~ r == 0.35 ~
      "   Shaded: ±1 SD"
    ),
    x = "n (log scale)",
    y = expression(Mean(widehat(IOV)))
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")# ==============================================================================
# Plot 2: Mittlerer Bias (skaliert mit n) - MCAR vs. MAR
# ==============================================================================

compare_df$mean_bias        <- compare_df$mean_IOV - true_iov
compare_df$mean_scaled_bias <- compare_df$n * compare_df$mean_bias

# Theoretische Bias-Kurve unter MCAR (pi = 0.75)
Sigma_mcar <- Sigma_Star(m = 3, p = 0.20, r = 0.35, pi = 0.75)
theory_bias <- data.frame(
  n    = seq(50, 1000, by = 1),
  bias = sapply(seq(50, 1000, by = 1),
                function(ni) -(4 / (3 * ni)) * sum(diag(Sigma_mcar)))
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
  scale_color_manual(
    name   = "Missingness",
    values = c("MCAR" = "black", "MAR (MNAR)" = "steelblue")
  ) +
  scale_linetype_manual(
    name   = "Missingness",
    values = c("MCAR" = "solid", "MAR (MNAR)" = "dashed")
  ) +
  scale_x_continuous(
    breaks = c(50, 100, 250, 500, 1000),
    trans  = "log10"
  ) +
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

# ==============================================================================
# Plot 3: Geschätzte CDF-Komponenten - MCAR vs. MAR (festes n = 500)
# ==============================================================================

# Hier brauchst du cdf_results aus simulation_MAR —
# ergänze analog zu simulation() ein cdf_results-Objekt in simulation_MAR(),
# dann kannst du direkt vergleichen. Alternativ mit den vorhandenen
# summary-Statistiken:

cdf_compare <- data.frame(
  category  = rep(paste0("f_", 0:2), 2),
  mechanism = rep(c("MCAR", "MAR (MNAR)"), each = 3),
  true_val  = rep(true_cdf, 2)
)

# MCAR: n=500, pi=0.75 — Index in mcar_df suchen
mcar_500 <- mcar_df[mcar_df$n == 500 & mcar_df$pi == 0.75, ]
mar_500  <- mar_df[mar_df$n == 500, ]

# Asymptotische SD unter MCAR
asymp_sd_mcar <- sqrt(diag(Sigma_mcar) / 500)

cdf_compare$mean <- c(
  # MCAR: aus cdf_df extrahieren
  as.numeric(cdf_df[
    cdf_df$scenario == mcar_idx[which(scenarios$n[mcar_idx] == 500 &
                                       scenarios$pi[mcar_idx] == 0.75)] &
    cdf_df$statistic == "mean",
    c("f_0", "f_1", "f_2")
  ]),
  # MAR: Platzhalter — ergänze wenn cdf_results in simulation_MAR verfügbar
  rep(NA, 3)
)

# Hinweis: Plot 3 setzt voraus dass simulation_MAR() ebenfalls
# cdf_results zurückgibt. Ergänze dazu in simulation_MAR():
#
#   cdf_results <- matrix(NA, nrow = n_reps, ncol = m, ...)
#   cdf_results[rep, ] <- CDF
#   # und im return():
#   list(summary = ..., cdf = rbind(colMeans(cdf_results), apply(...)))