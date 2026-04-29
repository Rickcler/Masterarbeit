# ==============================================================================
# 04_visuals.R
# Alle Plots: IOV / Skew / Cohen's K Vergleich, Bias, CLT, CDF, Verteilungen
# Voraussetzung: Masterarbeit.RData geladen (erzeugt von 02 + 03)
# ==============================================================================

source("00_setup.R")
load("Masterarbeit.RData")

# ------------------------------------------------------------------------------
# Hilfsfunktionen
# ------------------------------------------------------------------------------

#' Bereitet subset + x-Positionen für Vergleichs-Plots vor
#' @param data       combined_df
#' @param filter_fn  Eine dplyr-filter-Bedingung als Funktion
#' @param group_var  Variable für x-Gruppierung (Symbol, z.B. quote(pi))
prepare_comparison_subset <- function(data, filter_expr, group_var_name) {
  subset <- data %>% filter({{ filter_expr }})

  unique_groups <- unique(subset[, c("n", group_var_name, "type")])
  unique_groups <- unique_groups[order(
    unique_groups$n,
    unique_groups[[group_var_name]],
    factor(unique_groups$type, levels = c("Simulation", "Asymptotic"))
  ), ]

  key <- paste(unique_groups$n, unique_groups[[group_var_name]],
               unique_groups$type, sep = "_")
  x_positions        <- setNames(seq_len(nrow(unique_groups)), key)
  subset$x_pos       <- x_positions[paste(subset$n, subset[[group_var_name]],
                                           subset$type, sep = "_")]
  subset$x_pos       <- as.numeric(subset$x_pos)

  group_labels <- unique(subset[, c("n", group_var_name)])
  group_labels <- group_labels[order(group_labels$n,
                                     group_labels[[group_var_name]]), ]

  group_centers <- sapply(seq_len(nrow(group_labels)), function(i) {
    mean(subset$x_pos[
      subset$n == group_labels$n[i] &
      subset[[group_var_name]] == group_labels[[group_var_name]][i]
    ])
  })

  x_labels <- paste0("n = ", group_labels$n, "\n",
                     group_var_name, " = ", group_labels[[group_var_name]])

  list(subset = subset, group_centers = group_centers, x_labels = x_labels)
}

#' Generischer Vergleichs-Plot (IOV, Skew oder Cohen's K)
#' @param subset        Vorbereitetes subset mit x_pos
#' @param y_var         String: Spaltenname der y-Variable (z.B. "mean_IOV")
#' @param ymin_var      String: untere Fehlerbalken-Spalte
#' @param ymax_var      String: obere Fehlerbalken-Spalte
#' @param true_val      Wahrer Parameterwert (horizontale Linie)
#' @param group_centers Mittelpunkte der x-Gruppen
#' @param x_labels      Beschriftungen der x-Gruppen
#' @param title         Plot-Titel
#' @param y_label       y-Achsenbeschriftung
#' @param subtitle      Plot-Untertitel
#' @param ylim_offset   c(unten, oben) relativ zum wahren Wert
#' @param group_var     String: Variable für Shape-Ästhetik
comparison_plot <- function(subset, y_var, ymin_var, ymax_var,
                            true_val, group_centers, x_labels,
                            title, y_label, subtitle,
                            ylim_offset = c(-0.125, 0.075),
                            group_var = "pi") {

  n_groups   <- length(group_centers)
  rect_xmin  <- c(0.5, seq(2.5, by = 2, length.out = n_groups - 1))
  rect_xmax  <- c(seq(2.5, by = 2, length.out = n_groups - 1),
                  max(subset$x_pos) + 0.5)

  ggplot(subset, aes(x = x_pos, y = .data[[y_var]],
                     color = type,
                     shape = factor(.data[[group_var]]))) +
    annotate("rect",
             xmin  = rect_xmin[seq(1, n_groups, by = 2)],
             xmax  = rect_xmax[seq(1, n_groups, by = 2)],
             ymin  = -Inf, ymax = Inf,
             alpha = 0.05, fill = "gray90") +
    geom_hline(yintercept = true_val,
               color = "darkgreen", linetype = "dashed",
               linewidth = 1, alpha = 0.7) +
    geom_line(aes(group = interaction(n, .data[[group_var]])),
              color = "gray50", linetype = "dashed",
              alpha = 0.5,
              position = position_dodge(width = 0.2)) +
    geom_point(size = 3.5, position = position_dodge(width = 0.2)) +
    geom_errorbar(aes(ymin = .data[[ymin_var]], ymax = .data[[ymax_var]]),
                  width = 0.15, linewidth = 0.8,
                  position = position_dodge(width = 0.2)) +
    scale_x_continuous(breaks = group_centers, labels = x_labels,
                       expand = expansion(mult = 0.1)) +
    scale_color_manual(
      values = c("Asymptotic" = "#E41A1C", "Simulation" = "#377EB8"),
      name   = "Method"
    ) +
    scale_shape_manual(
      values = setNames(c(16, 17, 15, 18),
                        as.character(sort(unique(subset[[group_var]])))),
      name   = group_var
    ) +
    coord_cartesian(ylim = c(true_val + ylim_offset[1],
                             true_val + ylim_offset[2])) +
    labs(title = title, subtitle = subtitle,
         x = "Scenario", y = y_label) +
    theme_minimal() +
    theme(
      plot.title         = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle      = element_text(hjust = 0.5, color = "gray40"),
      axis.text.x        = element_text(angle = 0, hjust = 0.5, vjust = 1),
      legend.position    = "bottom",
      legend.box         = "vertical",
      legend.spacing.y   = unit(0.2, "cm"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.border       = element_rect(color = "gray80", fill = NA,
                                        linewidth = 0.5)
    )
}

# ------------------------------------------------------------------------------
# combined_df aufbauen
# ------------------------------------------------------------------------------

combined_df <- rbind(asymp_df, sim_df)

combined_df <- combined_df %>%
  mutate(
    true_IOV  = mapply(true_IOV,  m = m, p = p),
    true_Skew = mapply(true_Skew, m = m, p = p),
    true_C1   = 0
  )

# ==============================================================================
# Plot-Block A: MCAR / pi_h == 0, beide pi-Werte
# ==============================================================================

prep_A <- prepare_comparison_subset(
  combined_df,
  filter_expr   = (p == 0.3 & pi_h == 0),
  group_var_name = "pi"
)
sub_A         <- prep_A$subset
t_IOV_A       <- sub_A$true_IOV[1]
t_Skew_A      <- sub_A$true_Skew[1]
m_A <- sub_A$m[1]; p_A <- sub_A$p[1]; r_A <- sub_A$r[1]
sub_A$true_C1 <- 0

# IOV
print(comparison_plot(
  sub_A, "mean_IOV", "lower_IOV", "upper_IOV",
  true_val      = t_IOV_A,
  group_centers = prep_A$group_centers,
  x_labels      = prep_A$x_labels,
  y_label       = "IOV Value",
  title = sprintf("BinAR(1): m = %d, p = %.2f, r = %.2f  |  MCAR",
                          m_A, p_A, r_A),
  subtitle = "",
  ylim_offset   = c(-0.07, 0.03),
  group_var     = "pi"
))

# Skew
comparison_plot(
  sub_A, "mean_Skew", "lower_Skew", "upper_Skew",
  true_val      = t_Skew_A,
  group_centers = prep_A$group_centers,
  x_labels      = prep_A$x_labels,
  title         = "Comparison of Asymptotic vs Simulated Skewness Results",
  y_label       = "Skewness Value",
  subtitle      = sprintf("BinAR(1): m = %d, p = %.2f, r = %.2f  |  MCAR",
                          m_A, p_A, r_A),
  ylim_offset   = c(-0.15, 0.15),
  group_var     = "pi"
)

# Cohen's K
comparison_plot(
  sub_A, "mean_C", "lower_C", "upper_C",
  true_val      = 0,
  group_centers = prep_A$group_centers,
  x_labels      = prep_A$x_labels,
  title         = "Comparison of Asymptotic vs Simulated Cohen's κ lag(1)",
  y_label       = "Cohen's κ lag(1)",
  subtitle      = sprintf("BinAR(1): m = %d, p = %.2f, r = %.2f  |  MCAR",
                          m_A, p_A, r_A),
  ylim_offset   = c(-0.3, 0.3),
  group_var     = "pi"
)

# ==============================================================================
# Plot-Block B: Serielle Abhängigkeit in Missingness, pi == 0.75
# ==============================================================================

prep_B <- prepare_comparison_subset(
  combined_df,
  filter_expr    = (pi_h %in% c(0.2, 0.75) & pi == 0.75),
  group_var_name = "pi_h"
)
sub_B    <- prep_B$subset
t_IOV_B  <- sub_B$true_IOV[1]
t_Skew_B <- sub_B$true_Skew[1]
m_B <- sub_B$m[1]; p_B <- sub_B$p[1]; r_B <- sub_B$r[1]

# IOV
comparison_plot(
  sub_B, "mean_IOV", "lower_IOV", "upper_IOV",
  true_val      = t_IOV_B,
  group_centers = prep_B$group_centers,
  x_labels      = prep_B$x_labels,
  title         = "Comparison of Asymptotic vs Simulated IOV Results\n(Serially Dependent Missingness)",
  y_label       = "IOV Value",
  subtitle      = sprintf("BinAR(1): m = %d, p = %.2f, r = %.2f  |  π = 0.75",
                          m_B, p_B, r_B),
  ylim_offset   = c(-0.125, 0.075),
  group_var     = "pi_h"
)

# Skew
comparison_plot(
  sub_B, "mean_Skew", "lower_Skew", "upper_Skew",
  true_val      = t_Skew_B,
  group_centers = prep_B$group_centers,
  x_labels      = prep_B$x_labels,
  title         = "Comparison of Asymptotic vs Simulated Skewness Results\n(Serially Dependent Missingness)",
  y_label       = "Skewness Value",
  subtitle      = sprintf("BinAR(1): m = %d, p = %.2f, r = %.2f  |  π = 0.75",
                          m_B, p_B, r_B),
  ylim_offset   = c(-0.15, 0.15),
  group_var     = "pi_h"
)

# Cohen's K
comparison_plot(
  sub_B, "mean_C", "lower_C", "upper_C",
  true_val      = 0,
  group_centers = prep_B$group_centers,
  x_labels      = prep_B$x_labels,
  title         = "Comparison of Asymptotic vs Simulated Cohen's κ lag(1)\n(Serially Dependent Missingness)",
  y_label       = "Cohen's κ lag(1)",
  subtitle      = sprintf("BinAR(1): m = %d, p = %.2f, r = %.2f  |  π = 0.75",
                          m_B, p_B, r_B),
  ylim_offset   = c(-0.3, 0.3),
  group_var     = "pi_h"
)

# ==============================================================================
# Plot-Block C: Bias und CLT
# ==============================================================================

# Theoretische Bias-Kurve
theory_df_2 <- data.frame(
  n    = seq(2, 1000, by = 1),
  diff = sapply(seq(2, 1000, by = 1), function(n_val)
    -(4 / (10 * n_val)) * sum(diag(Sigma_raw)))
)

mean_df_2 <- sim_data_2 %>%
  dplyr::group_by(n) %>%
  dplyr::summarise(
  mean_scaled_bias = mean(scaled_bias, na.rm = TRUE),
  mean_diff        = mean(diff,        na.rm = TRUE),
  se = sd(scaled_bias, na.rm = TRUE) / sqrt(n()),
  .groups = "drop"
)

# Plot C1: Mittlerer Bias
ggplot(mean_df_2, aes(x = n, y = mean_diff)) +
  geom_line(data = theory_df_2, aes(x = n, y = diff),
            color = "steelblue", linetype = "dashed", linewidth = 0.8) +
  geom_line() +
  geom_point(size = 2) +
  labs(
    title    = "Mean bias as a function of n",
    subtitle = "Dashed: theoretical bias  −4/(mn) · tr(Σ)\nBinAR(1): m = 10, p = 0.3, r = 0.2, π = 0.75",
    x = "n",
    y = expression(Mean(hat(IOV) - IOV))
  ) +
  theme_minimal()



# Plot C3: CLT-Dichte
ggplot(sim_data, aes(x = scaled_CLT, fill = factor(n))) +
  geom_density(alpha = 0.4) +
  stat_function(
    fun  = dnorm,
    args = list(mean = mean(sim_data$scaled_CLT),
                sd   = sd(sim_data$scaled_CLT)),
    linetype = "dashed"
  ) +
  labs(fill = "n", title = expression(sqrt(n) ~ (hat(IOV) - IOV))) +
  theme_minimal()

# Plot C4: Bias-Dichte
ggplot(sim_data, aes(x = scaled_bias, fill = factor(n))) +
  geom_density(alpha = 0.4) +
  stat_function(
    fun  = dnorm,
    args = list(mean = mean(sim_data$scaled_bias),
                sd   = sd(sim_data$scaled_bias)),
    linetype = "dashed"
  ) +
  labs(fill = "n", title = expression(n ~ (hat(IOV) - IOV))) +
  theme_minimal()

# Plot C5: CLT-Dichte unscaled
ggplot(sim_data, aes(x = diff, fill = factor(n))) +
  geom_density(alpha = 0.4) +
  stat_function(
    fun  = dnorm,
    args = list(mean = mean(sim_data$diff),
                sd   = sd(sim_data$diff)),
    linetype = "dashed"
  ) +
  labs(fill = "n", title = expression(sqrt(n) ~ (hat(IOV) - IOV))) +
  theme_minimal()

# ==============================================================================
# Plot-Block D: Geschätzte CDF-Komponenten
# ==============================================================================

true_cdf <- pbinom(0:2, size = 3, prob = 0.20)
true_df  <- data.frame(
  category = paste0("f_", 0:2),
  true_cdf = true_cdf
)

Sigma_cdf <- Sigma_Star(m = 3, p = 0.20, r = 0.35, pi = 0.75, pi_h = 0)

est_df <- cdf_df %>%
  pivot_longer(cols = starts_with("f_"),
               names_to = "category", values_to = "value") %>%
  pivot_wider(names_from = statistic, values_from = value) %>%
  mutate(n = factor(n, levels = sort(unique(as.numeric(as.character(n))))))

asymp_cdf_df <- do.call(rbind, lapply(unique(est_df$n), function(ni) {
  ni_num <- as.numeric(as.character(ni))
  data.frame(
    category = paste0("f_", 0:2),
    n        = factor(ni, levels = levels(est_df$n)),
    asymp_sd = sqrt(diag(Sigma_cdf) / ni_num),
    true_cdf = true_cdf
  )
}))

plot_cdf_df <- bind_rows(
  est_df %>%
    mutate(type = "Simulation", ymin = mean - sd, ymax = mean + sd),
  asymp_cdf_df %>%
    rename(mean = true_cdf) %>%
    mutate(type = "Asymptotic",
           ymin = mean - asymp_sd,
           ymax = mean + asymp_sd)
) %>%
  mutate(type = factor(type, levels = c("Simulation", "Asymptotic")))

ggplot(plot_cdf_df, aes(x = category, color = n, linetype = type)) +
  geom_bar(
    data = true_df,
    aes(x = category, y = true_cdf),
    inherit.aes = FALSE,
    stat = "identity", fill = "grey85", color = "grey60", width = 0.7
  ) +
  geom_linerange(aes(ymin = ymin, ymax = ymax),
                 position = position_dodge(width = 0.7),
                 linewidth = 0.8) +
  geom_point(aes(y = mean),
             position = position_dodge(width = 0.7),
             size = 2) +
  scale_color_grey(name = "n", start = 0.7, end = 0.1) +
  scale_linetype_manual(
    name   = "Std. Deviation",
    values = c("Simulation" = "solid", "Asymptotic" = "dotted")
  ) +
  guides(linetype = guide_legend(
    override.aes = list(linewidth = 1, size = 0, color = "black")
  )) +
  scale_y_continuous(limits = c(0, 1), name = "Cumulative Probability") +
  scale_x_discrete(name = "CDF Component") +
  labs(
    title    = "Estimated vs. True CDF",
    subtitle = expression(m == 3 ~ "," ~ p == 0.20 ~ "," ~ r == 0.35 ~
                          "," ~ pi == 0.75)
  ) +
  theme_minimal()

# ==============================================================================
# Plot-Block E: Marginale Verteilungen der Simulationsszenarien
# ================ge==============================================================

# Masterarbeit
scen_A <- data.frame(
  scenario = "Scenario A\nm = 3, p = 0.20, r = 0.35",
  category = 0:3,
  pmf      = dbinom(0:3,  size = 3,  prob = 0.20),
  cdf      = pbinom(0:3,  size = 3,  prob = 0.20)
)

scen_B <- data.frame(
  scenario = "Scenario B\nm = 10, p = 0.45, r = 0.50",
  category = 0:10,
  pmf      = dbinom(0:10, size = 10, prob = 0.45),
  cdf      = pbinom(0:10, size = 10, prob = 0.45)
)

plot_dist_df <- bind_rows(scen_A, scen_B) %>%
  mutate(
    scenario = factor(scenario, levels = c(
      "Scenario A\nm = 3, p = 0.20, r = 0.35",
      "Scenario B\nm = 10, p = 0.45, r = 0.50"
    )),
    category = factor(category)
  )

p_pmf <- ggplot(plot_dist_df, aes(x = category, y = pmf)) +
  geom_col(fill = "grey70", color = "grey40", width = 0.6) +
  geom_text(aes(label = round(pmf, 3)),
            vjust = -0.4, size = 2.8, color = "grey30") +
  facet_wrap(~ scenario, scales = "free_x") +
  scale_y_continuous(limits = c(0, 0.60), name = "Probability") +
  scale_x_discrete(name = "Category") +
  labs(title = "Marginal PMF") +
  theme_minimal() +
  theme(strip.text = element_text(size = 10, face = "bold"),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(size = 11))

p_cdf_dist <- ggplot(plot_dist_df, aes(x = category, y = cdf)) +
  geom_col(fill = "grey70", color = "grey40", width = 0.6) +
  geom_text(aes(label = round(cdf, 3)),
            vjust = -0.4, size = 2.8, color = "grey30") +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             color = "steelblue", linewidth = 0.6) +
  annotate("text", x = 0.6, y = 0.52,
           label = "0.5", color = "steelblue", size = 3) +
  facet_wrap(~ scenario, scales = "free_x") +
  scale_y_continuous(limits = c(0, 1.05), name = "Cumulative Probability") +
  scale_x_discrete(name = "Category") +
  labs(title = "Marginal CDF") +
  theme_minimal() +
  theme(strip.text = element_text(size = 10, face = "bold"),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(size = 11))

p_pmf / p_cdf_dist +
  plot_annotation(
    title    = "Marginal distributions of the two simulation scenarios",
    subtitle = "Top: PMF    Bottom: CDF    Dashed line: 0.5",
    theme    = theme(
      plot.title    = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "grey40")
    )
  )
# Vortrag

# Farbpalette für Vortrag
col_primary   <- "#2C3E6B"   # dunkles Blau
col_secondary <- "#5B8DB8"   # mittleres Blau
col_accent    <- "#E8EEF4"   # sehr helles Blau (Hintergrund Balken)
col_text      <- "#2C2C2C"
col_grid      <- "#E5E5E5"


scen_A <- data.frame(
  scenario = "Scenario A",
  category = 0:3,
  pmf      = dbinom(0:3,  size = 3,  prob = 0.20),
  cdf      = pbinom(0:3,  size = 3,  prob = 0.20)
)

scen_B <- data.frame(
  scenario = "Scenario B",
  category = 0:10,
  pmf      = dbinom(0:10, size = 10, prob = 0.45),
  cdf      = pbinom(0:10, size = 10, prob = 0.45)
)

plot_dist_df <- bind_rows(scen_A, scen_B) %>%
  mutate(
    scenario = factor(scenario, levels = c("Scenario A", "Scenario B")),
    category = factor(category)
  )
# Einheitliches Vortrag-Theme
theme_talk <- function() {
  theme_minimal(base_size = 14) +
    theme(
      # Hintergrund
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = col_grid, linewidth = 0.4),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      # Titel
      plot.title    = element_text(size = 16, face = "bold",
                                   color = col_primary, margin = margin(b = 4)),
      plot.subtitle = element_text(size = 12, color = "grey50",
                                   margin = margin(b = 10)),
      # Achsen
      axis.title    = element_text(size = 12, color = col_text),
      axis.text     = element_text(size = 11, color = col_text),
      axis.ticks    = element_blank(),
      # Facetten
      strip.text    = element_text(size = 13, face = "bold",
                                   color = col_primary),
      strip.background = element_rect(fill = col_accent, color = NA),
      # Legende
      legend.position  = "bottom",
      legend.title     = element_text(size = 11, face = "bold"),
      legend.text      = element_text(size = 11),
      # Ränder
      plot.margin = margin(12, 16, 12, 16)
    )
}

p_pmf <- ggplot(plot_dist_df, aes(x = category, y = pmf)) +
  geom_col(fill = col_secondary, color = "white",
           width = 0.65, alpha = 0.85) +
  geom_text(aes(label = sprintf("%.3f", pmf)),
            vjust = -0.5, size = 3.5,
            color = col_primary, fontface = "bold") +
  facet_wrap(~ scenario, scales = "free_x") +
  scale_y_continuous(
    limits = c(0, 0.68),
    expand = c(0, 0),
    name   = "Probability"
  ) +
  scale_x_discrete(name = "Category") +
  labs(title = "Marginal PMF") +
  theme_talk()

p_cdf_dist <- ggplot(plot_dist_df, aes(x = category, y = cdf)) +
  geom_col(fill = col_secondary, color = "white",
           width = 0.65, alpha = 0.85) +
  geom_text(aes(label = sprintf("%.3f", cdf)),
            vjust = -0.5, size = 3.5,
            color = col_primary, fontface = "bold") +
  annotate("text",
           x = 0.55, y = 0.54,
           label = "0.5",
           color = "#C0392B",
           size  = 4,
           fontface = "bold") +
  facet_wrap(~ scenario, scales = "free_x") +
  scale_y_continuous(
    limits = c(0, 1.10),
    expand = c(0, 0),
    name   = "Cumulative Probability"
  ) +
  scale_x_discrete(name = "Category") +
  labs(title = "Marginal CDF") +
  theme_talk()

p_pmf / p_cdf_dist +
  plot_annotation(
    theme = theme(
      plot.title = element_text(
        size     = 17,
        face     = "bold",
        color    = col_primary,
        margin   = margin(b = 8)
      ),
      plot.background = element_rect(fill = "white", color = NA)
    )
  ) &
  theme(plot.background = element_rect(fill = "white", color = NA))

ggsave("Graphs/marginal_distributions_talk.pdf",
       width = 11, height = 7, units = "in", dpi = 300)

