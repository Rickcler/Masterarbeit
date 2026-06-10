# ==============================================================================
# 04_visuals.R
# Alle Plots: IOV / Skew / Cohen's K Vergleich, Bias, CLT, CDF, Verteilungen
# Voraussetzung: Masterarbeit.RData geladen (erzeugt von 02 + 03)
# ==============================================================================

source("00_setup.R")
load("Masterarbeit.RData")
#------------------------------------------------------------------------------
# Figure 4.1: Cohen's Kappa of BinAR(1) for different values of r
#------------------------------------------------------------------------------

m <- 10         # Anzahl der Versuche (Zustandsraum {0,...,m})
p <- 0.45        # Erfolgswahrscheinlichkeit der marginalen Binomialverteilung
r_vals <- c(-0.5, -0.2, 0.2, 0.4, 0.6, 0.8, 0.95)   # verschiedene Autokorrelationsparameter
max_h <- 100  

kappa_matrix <- matrix(NA, nrow = max_h, ncol = length(r_vals))
colnames(kappa_matrix) <- paste0("r = ", r_vals)

for (i in seq_along(r_vals)) {
  r <- r_vals[i]
  for (h in 1:max_h) {
    kappa_matrix[h, i] <- kappa_ord(m, p, r, h)
  }
}

h_werte <- 1:max_h 
df <- data.frame(h = h_werte, kappa_matrix, check.names = FALSE)   # Spalten: h, r1, r2, ...


df_long <- pivot_longer(df, 
                        cols = -h,          # alle Spalten außer 'h' werden umgeformt
                        names_to = "r", 
                        values_to = "kappa")

# 2. Plot
Cohens_Plot <- ggplot(df_long, aes(x = h, y = kappa, color = r)) +
  geom_line() +          # Linien
  geom_point() +         # Punkte (optional)
  labs(x = "h", 
       y = "\u03BA(h)",
       color = "r-Values",
       title = "\u03BA(h) in relation to h",
       subtitle = "BinAR(1) with p = 0.45, m = 10 over multiple r-Values") +
  theme_minimal() +
  coord_fixed(ratio = 40)
print(Cohens_Plot)
ggsave("Graphs/Cohens_Kappa_by_h.png", Cohens_Plot + theme(legend.position = "none"), width = 8, height = 6)

#------------------------------------------------------------------------------
# Figure 4.2: Maginal Distribution 
#------------------------------------------------------------------------------

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

MarginalsPlot <- p_pmf / p_cdf_dist +
  plot_annotation(
    title    = "Marginal distributions of the two simulation scenarios",
    subtitle = "Top: PMF    Bottom: CDF    Dashed line: 0.5",
    theme    = theme(
      plot.title    = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "grey40")
    )
  )

print(MarginalsPlot)
ggsave("Graphs/MarginalsPlot.png", MarginalsPlot, width = 8, height = 8)
#------------------------------------------------------------------------------
# Figures 4.3 and 4.4: Bias and CLT 
#------------------------------------------------------------------------------
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
  se = sd(scaled_bias, na.rm = TRUE),
  .groups = "drop"
) %>%
  dplyr::mutate(mean_scaled_diff = (mean_diff * n))
Marginal <- pbinom(0:9, size = 10, prob = 0.3)


Var_Iov <- 0
for (i in 1:10) {
  for (j in 1:10) {
    Var_Iov <- Var_Iov + (1- 2 * Marginal[i])*(1- 2 * Marginal[j])*Sigma_raw[i, j]
  }
}

SD_Iov <- sqrt((16 / 10^2) * Var_Iov)


# Figure 4.3 left: CLT-Dichte unscaled
unscaled_CLT_plot <- ggplot(sim_data, aes(x = diff, fill = factor(n))) +
  geom_density(alpha = 0.4) +
  labs(x = expression(Bias(IOV))) +
  theme_minimal() +
  theme(legend.position = "none")
print(unscaled_CLT_plot)
ggsave("Graphs/CLT_unscaled.png", unscaled_CLT_plot + theme(legend.position = "none"), width = 8, height = 5)


# Figure 4.3 right: CLT-Dichte scaled
scaled_CLT_plot <- ggplot(sim_data, aes(x = scaled_CLT, fill = factor(n))) +
  geom_density(alpha = 0.4) +
  stat_function(
    fun  = dnorm,
    args = list(mean = 0,
                sd   = SD_Iov),
    linetype = "dashed"
  ) +
  labs(x = expression(Bias(IOV))) +
  theme_minimal() +
  theme(legend.position = "none")
print(scaled_CLT_plot)
ggsave("Graphs/CLT_scaled.png", scaled_CLT_plot + theme(legend.position = "none"), width = 8, height = 5)


# Figure 4.4 left: Mittlerer Bias unscaled
unscaled_bias_plot <- ggplot(mean_df_2, aes(x = n, y = mean_diff)) +
  geom_line(data = theory_df_2, aes(x = n, y = diff),
            color = "steelblue", linetype = "dashed", linewidth = 0.8) +
  geom_line() +
  geom_point(size = 2) +
  labs(
    title    = "Mean bias as a function of n",
    subtitle = "",
    x = "n",
    y = expression(Mean(hat(IOV) - IOV))
  ) +
  theme_minimal()

print(unscaled_bias_plot)
ggsave("Graphs/Bias_unscaled.png", unscaled_bias_plot + theme(legend.position = "none"), width = 8, height = 5)


# Figure 4.4 right: Mittlerer Bias scaled

scaled_bias_plot <- ggplot(mean_df_2, aes(x = n, y = mean_scaled_diff)) +
  
  # Punkte
  geom_point(size = 1.5, color = "black") +
  
  # schwarze Linie (empirisch)
  geom_line(color = "black", linewidth = 0.5) +
  
  # rote Glättung (wie im Plot)
  geom_smooth(method = "loess", se = FALSE,
              color = "red", linewidth = 1) +
  
  # theoretischer Grenzwert (horizontal)
  geom_hline(
    yintercept = -(4 / 10) * sum(diag(Sigma_raw)),
    linetype = "dashed",
    color = "steelblue"
  ) +
  
  labs(
    x = "n",
    y = expression(n %.% Mean(hat(IOV) - IOV))
  ) +
  
  theme_minimal()

print(scaled_bias_plot)
ggsave("Graphs/Bias_scaled.png", scaled_bias_plot + theme(legend.position = "none"), width = 8, height = 5)

# ------------------------------------------------------------------------------
# Figures 4.5 - 4.8: Comparisonplots IOV, Skew, Cohen's K 
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

  x_labels <- mapply(
  function(n_val, g_val) {

    if (group_var_name == "pi") {
      bquote(atop(n == .(n_val),
                   pi == .(g_val)))

    } else if (group_var_name == "r_pi") {
      bquote(atop(n == .(n_val),
                   r[pi] == .(g_val)))

    } else if (group_var_name == "pi_h") {
      bquote(atop(n == .(n_val),
                   pi[h] == .(g_val)))

    } else {
      bquote(atop(n == .(n_val),
                   .(group_var_name) == .(g_val)))
    }

  },
  group_labels$n,
  group_labels[[group_var_name]],
  SIMPLIFY = FALSE
)

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
  legend_label <- switch(
    group_var,
    "pi"   = expression(pi),
    "r_pi" = expression(r[pi]),
    "pi_h" = expression(pi[h]),
    group_var
  )
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
      name   = legend_label
    ) +
    coord_cartesian(ylim = c(true_val + ylim_offset[1],
                             true_val + ylim_offset[2])) +
    labs(title = title, subtitle = subtitle,
         x = "", y = "") +
    theme_minimal() +
    theme(
      plot.title         = element_text(hjust = 0.5, face = "bold", size = 20),
      plot.subtitle      = element_text(hjust = 0.5, color = "gray40"),
      axis.text.x        = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 12, color = "gray20"),
      axis.text.y        = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 14, color = "gray20"),
      legend.position    = "bottom",
      legend.box         = "horizontal",      # Legendenblöcke nebeneinander
      legend.direction   = "horizontal",        # Items innerhalb je vertikal
      legend.spacing.x   = unit(1, "cm"),     # Abstand zwischen den Blöcken
      legend.text        = element_text(size = 18),
      legend.title       = element_text(size = 20, face = "bold"),
      legend.key.size    = unit(0.6, "cm"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.border       = element_rect(color = "gray80", fill = NA,
                                        linewidth = 0.5)
    ) +
    guides(
      color = guide_legend(order = 1, override.aes = list(size = 4)),
      shape = guide_legend(order = 2, override.aes = list(size = 4))
    )
}


# Combined_df aufbauen


combined_df <- rbind(asymp_df, sim_df)

combined_df <- combined_df %>%
  mutate(
    true_IOV  = mapply(true_IOV,  m = m, p = p),
    true_Skew = mapply(true_Skew, m = m, p = p),
    true_C1   = 0
  )


# Plot-Block A (Figure 4.5, 4.7, 4.8): Vergleichsplots für pi in {1, 0.75} und r_pi == 0 (MCAR)

### Scenario A

prep_A1 <- prepare_comparison_subset(
  combined_df,
  filter_expr   = (p == 0.2 & r_pi == 0), # hier p 0.2 oder 0.45 wählen
  group_var_name = "pi"
)

sub_A1         <- prep_A1$subset
t_IOV_A1       <- sub_A1$true_IOV[1]
t_Skew_A1      <- sub_A1$true_Skew[1]
m_A1 <- sub_A1$m[1]; p_A1 <- sub_A1$p[1]; r_A1 <- sub_A1$r[1]
sub_A1$true_C1 <- 0
sub_A1$IOV_centered <- sub_A1$mean_IOV - t_IOV_A1
sub_A1$IOV_lower_centered <- sub_A1$lower_IOV - t_IOV_A1
sub_A1$IOV_upper_centered <- sub_A1$upper_IOV - t_IOV_A1

# IOV
IOV_plot <- comparison_plot(
  sub_A1, "IOV_centered", "IOV_lower_centered", "IOV_upper_centered",
  true_val      = 0,
  group_centers = prep_A1$group_centers,
  x_labels      = prep_A1$x_labels,
  y_label       = "IOV Value",
  title = "Scenario A",
  subtitle = sprintf("BinAR(1): m = %d, p = %.2f, r = %.2f  |  MCAR",
                          m_A1, p_A1, r_A1),
  ylim_offset   = c(-0.10, 0.06), # anpassen je nach p
  group_var     = "pi"
)
print(IOV_plot)
ggsave(sprintf("Graphs/iov_m%d.png", m_A1), IOV_plot + theme(legend.position = "none"), width = 8, height = 5)

# Skew
Skew_plot <- comparison_plot(
  sub_A1, "mean_Skew", "lower_Skew", "upper_Skew",
  true_val      = t_Skew_A1,
  group_centers = prep_A1$group_centers,
  x_labels      = prep_A1$x_labels,
  title         = "Scenario A",
  y_label       = "Skewness Value",
  subtitle      = sprintf("BinAR(1): m = %d, p = %.2f, r = %.2f  |  MCAR",
                          m_A1, p_A1, r_A1),
  ylim_offset   = c(-0.10, 0.1), # anpassen je nach p
  group_var     = "pi"
)
print(Skew_plot)
ggsave(sprintf("Graphs/skew_m%d.png", m_A1), Skew_plot + theme(legend.position = "none"), width = 8, height = 5)

Cohens_plot <- comparison_plot(
  sub_A1, "mean_C", "lower_C", "upper_C",
  true_val      = 0,
  group_centers = prep_A1$group_centers,
  x_labels      = prep_A1$x_labels,
  title         = "Scenario A",
  y_label       = "Cohen's κ lag(1)",
  subtitle      = sprintf("BinAR(1): m = %d, p = %.2f, r = 0  |  MCAR",
                          m_A1, p_A1, r_A1),
  ylim_offset   = c(-0.275, 0.26), # anpassen je nach p
  group_var     = "pi"
)
print(Cohens_plot)
ggsave(sprintf("Graphs/kappa_m%d.png", m_A1), Cohens_plot + theme(legend.position = "none"), width = 8, height = 5)


### Scenario B

prep_B1 <- prepare_comparison_subset(
  combined_df,
  filter_expr   = (p == 0.45 & r_pi == 0), # hier p 0.2 oder 0.45 wählen
  group_var_name = "pi"
)
sub_B1         <- prep_B1$subset
t_IOV_B1       <- sub_B1$true_IOV[1]
t_Skew_B1      <- sub_B1$true_Skew[1]
m_B1 <- sub_B1$m[1]; p_B1 <- sub_B1$p[1]; r_B1 <- sub_B1$r[1]
sub_B1$true_C1 <- 0
sub_B1$IOV_centered <- sub_B1$mean_IOV - t_IOV_B1
sub_B1$IOV_lower_centered <- sub_B1$lower_IOV - t_IOV_B1
sub_B1$IOV_upper_centered <- sub_B1$upper_IOV - t_IOV_B1

# IOV
IOV_plot <- comparison_plot(
  sub_B1, "IOV_centered", "IOV_lower_centered", "IOV_upper_centered",
  true_val      = 0,
  group_centers = prep_B1$group_centers,
  x_labels      = prep_B1$x_labels,
  y_label       = "IOV Value",
  title = "Scenario B",
  subtitle = sprintf("BinAR(1): m = %d, p = %.2f, r = %.2f  |  MCAR",
                     m_B1, p_B1, r_B1),
  ylim_offset   = c(-0.10, 0.06), # anpassen je nach p
  group_var     = "pi"
)
print(IOV_plot)
ggsave(sprintf("Graphs/iov_m%d.png", m_B1), IOV_plot + theme(legend.position = "none"), width = 8, height = 5)

# Skew
Skew_plot <- comparison_plot(
  sub_B1, "mean_Skew", "lower_Skew", "upper_Skew",
  true_val      = t_Skew_B1,
  group_centers = prep_B1$group_centers,
  x_labels      = prep_B1$x_labels,
  title         = "Scenario B",
  subtitle      = sprintf("BinAR(1): m = %d, p = %.2f, r = %.2f  |  MCAR",
                          m_B1, p_B1, r_B1),
  y_label       = "Skewness Value",
  ylim_offset   = c(-0.10, 0.1), # anpassen je nach p
  group_var     = "pi"
)
print(Skew_plot)
ggsave(sprintf("Graphs/skew_m%d.png", m_B1), Skew_plot + theme(legend.position = "none"), width = 8, height = 5)

Cohens_plot <- comparison_plot(
  sub_B1, "mean_C", "lower_C", "upper_C",
  true_val      = 0,
  group_centers = prep_B1$group_centers,
  x_labels      = prep_B1$x_labels,
  title         = "Scenario B",
  y_label       = "Cohen's κ lag(1)",
  subtitle      = sprintf("BinAR(1): m = %d, p = %.2f, r = 0  |  MCAR",
                          m_B1, p_B1, r_B1),
  ylim_offset   = c(-0.275, 0.26), # anpassen je nach p
  group_var     = "pi"
)
print(Cohens_plot)
ggsave(sprintf("Graphs/kappa_m%d.png", m_B1), Cohens_plot + theme(legend.position = "none"), width = 8, height = 5)

# Legende
p <- Skew_plot
legend <- get_legend(
  p + theme(legend.position = "bottom")
)
legend_plot_IOV_sd <- ggdraw(legend)

ggsave(
  "Graphs/legend_plot_IOV_Skew.png",
  legend_plot_IOV_sd,
  width = 8,
  height = 1
)


# Plot-Block B (Figure 4.6): IOV with serial dependence in Missingness, pi == 0.75


### Scenario A
prep_A2 <- prepare_comparison_subset(
  combined_df,
  filter_expr    = (r_pi %in% c(0.2, 0.75) & pi == 0.75 & m == 3),
  group_var_name = "r_pi"
)
sub_A2 <- prep_A2$subset
t_IOV_A2  <- sub_A2$true_IOV[1]
t_Skew_A2 <- sub_A2$true_Skew[1]
m_A2 <- sub_A2$m[1]; p_A2 <- sub_A2$p[1]; r_A2 <- sub_A2$r[1]

sub_A2$IOV_centered <- sub_A2$mean_IOV - t_IOV_A2
sub_A2$IOV_lower_centered <- sub_A2$lower_IOV - t_IOV_A2  
sub_A2$IOV_upper_centered <- sub_A2$upper_IOV - t_IOV_A2
# IOV
IOV_sd_plot <- comparison_plot(
  sub_A2, "IOV_centered", "IOV_lower_centered", "IOV_upper_centered",
  true_val      = 0,
  group_centers = prep_A2$group_centers,
  x_labels      = prep_A2$x_labels,
  title         = "Scenario A", # C
  y_label       = "IOV Value",
  subtitle      = sprintf("BinAR(1): m = %d, p = %.2f, r = %.2f  |  π = 0.75",
                          m_A2, p_A2, r_A2),
  ylim_offset   = c(-0.10, 0.06), # m = 3 -> , m = 10 -> c(-0.07, 0.035)
  group_var     = "r_pi"
)
print(IOV_sd_plot + theme(legend.position = "none"))
ggsave(sprintf("Graphs/IOV_sd_m%d.png", m_A2), IOV_sd_plot + theme(legend.position = "none"), width = 8, height = 5)


### Scenario B
prep_B2 <- prepare_comparison_subset(
  combined_df,
  filter_expr    = (r_pi %in% c(0.2, 0.75) & pi == 0.75 & m == 10),
  group_var_name = "r_pi"
)
sub_B2    <- prep_B2$subset
t_IOV_B2  <- sub_B2$true_IOV[1]
t_Skew_B2 <- sub_B2$true_Skew[1]
m_B2 <- sub_B2$m[1]; p_B2 <- sub_B2$p[1]; r_B2 <- sub_B2$r[1]

sub_B2$IOV_centered <- sub_B2$mean_IOV - t_IOV_B2
sub_B2$IOV_lower_centered <- sub_B2$lower_IOV - t_IOV_B2
sub_B2$IOV_upper_centered <- sub_B2$upper_IOV - t_IOV_B2
# IOV
IOV_sd_plot <- comparison_plot(
  sub_B2, "IOV_centered", "IOV_lower_centered", "IOV_upper_centered",
  true_val      = 0,
  group_centers = prep_B2$group_centers,
  x_labels      = prep_B2$x_labels,
  title         = "Scenario B", # C
  y_label       = "IOV Value",
  subtitle      = sprintf("BinAR(1): m = %d, p = %.2f, r = %.2f  |  π = 0.75",
                          m_B2, p_B2, r_B2),
  ylim_offset   = c(-0.10, 0.06), # m = 3 -> , m = 10 -> c(-0.07, 0.035)
  group_var     = "r_pi"
)
print(IOV_sd_plot)
ggsave(sprintf("Graphs/IOV_sd_m%d.png", m_B2), IOV_sd_plot + theme(legend.position = "none"), width = 8, height = 5)



# Legende
p <- IOV_sd_plot
legend <- get_legend(
  p + theme(legend.position = "bottom")
)
legend_plot_IOV_sd <- ggdraw(legend)

ggsave(
  "Graphs/legend_plot_IOV_sd.png",
  legend_plot_IOV_sd,
  width = 8,
  height = 1
)

#====================================================
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
  scale_x_discete(name = "Category") +
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

