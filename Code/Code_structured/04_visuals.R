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
ggsave("Graphs/Cohens_Kappa_by_h.png", Cohens_Plot, width = 8, height = 6)

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

band_df_2 <- sim_data_2 %>%
  group_by(n) %>%
  summarise(
    B                = sum(!is.na(IOV)),
    mean_scaled_diff = unique(n) * mean(IOV - IOV_true, na.rm = TRUE),
    sd_IOV           = sd(IOV, na.rm = TRUE),
    se_scaled_diff   = unique(n) * sd_IOV / sqrt(B)
  ) %>%
  ungroup()

IOV_true <- (4 / 10) * sum(f_true * (1 - f_true))
scaled_bias_plot <- ggplot(band_df_2, aes(x = n, y = mean_scaled_diff)) +
  
  # 95%-Konfidenzintervall des Monte-Carlo-Mittelwerts (vertikale Linien)
  geom_errorbar(
    aes(ymin = mean_scaled_diff - 1.96 * se_scaled_diff,
        ymax = mean_scaled_diff + 1.96 * se_scaled_diff),
    width = 0, color = "grey50"
  ) +
  
  # Punkte
  geom_point(size = 1.5, color = "black") +
  
  # schwarze Linie (empirisch)
  geom_line(color = "black", linewidth = 0.5) +
  
  # rote Glättung
  geom_smooth(method = "loess", se = FALSE,
              color = "red", linewidth = 1) +
  
  # theoretischer Grenzwert
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



# Combined_df aufbauen


combined_df <- rbind(asymp_df, sim_df)

combined_df <- combined_df %>%
  mutate(
    true_IOV  = mapply(true_IOV,  m = m, p = p),
    true_Skew = mapply(true_Skew, m = m, p = p),
    true_C1   = 0
  )


# Plot-Block A (Figure 4.5, 4.7, 4.8): Comparisonplots for pi in {1, 0.75} and r_pi == 0 (MCAR)

### Scenario A

##### Datapreparation 
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



# IOV (Figure 4.5)
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

# Skew (Figure 4.7)
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

# Cohen's κ (Figure 4.8)
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

##### Datapreparation

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

# IOV (Figure 4.5)
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

# Skew (Figure 4.7)
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

# Cohen's κ (Figure 4.8)
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

# Legende (Figures 4.5, 4.7, 4.8)
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

##### Datapreparation
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


# IOV (Figure 4.6)
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

##### Datapreparation
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

# IOV (Figure 4.6)
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
# ==============================================================================
# Figure 4.9: Rejection rates under H_0 for uncorrected vs. bias-corrected test statistic
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

# ------------------------------------------------------------------------------
# Figure 4.10: Consistency of Cohen's κ under H_A 
# ------------------------------------------------------------------------------

### Left Plot (with Confidence Ribbons)
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


### Right Plot (Rejection Rates)

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


# ------------------------------------------------------------------------------
# Figure 4.11: IOV under MNAR 
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

MNAR_IOV_plt <- ggplot(compare_df,
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
print(MNAR_IOV_plt)

ggsave(
  "Graphs/iov_mnar.png",
  MNAR_IOV_plt,
  width = 8,
  height = 5
)
# ------------------------------------------------------------------------------
# Figure 4.12: Serial dependence under MNAR 
# ------------------------------------------------------------------------------

# Left Plot - Confidence Intervalls 
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
    name = ""
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
    plot.title = element_text(
      size = 13,
      face = "bold"
    ),
    axis.text.x        = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 12, color = "gray20"),
    axis.text.y        = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 12, color = "gray20"),
    plot.subtitle = element_text(
      size = 9,
      color = "grey40"
    ))

ggsave("Graphs/MNAR_Kappa_H_A.png", width = 5.5, height = 8)

# Right Plot - Rejection Rates 
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
    name   = ""
  ) +

  scale_y_continuous(
    limits = c(0, 1.05),
    breaks = seq(0, 1, by = 0.25),
    labels = percent,
    name   = "Rejection rate"
  ) +

  labs(
    title = "",

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
