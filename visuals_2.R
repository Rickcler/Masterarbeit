library(ggplot2)
library(dplyr)
library(tidyr)


# ==============================================================================
# Helper: Shared Plot Theme and Scenario Preparation
# ==============================================================================

thesis_theme <- function() {
  theme_minimal() +
    theme(
      plot.title       = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle    = element_text(hjust = 0.5, color = "gray40"),
      axis.text.x      = element_text(angle = 0, hjust = 0.5, vjust = 1),
      legend.position  = "bottom",
      legend.box       = "vertical",
      legend.spacing.y = unit(0.2, "cm"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.border       = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
    )
}

#' Prepare combined asymptotic + simulated data for a given scenario
#' Adds true parameter values, x-positions, and axis labels
#' @param p_val  Value of p to filter on (0.20 or 0.45)
prepare_subset <- function(combined_df, p_val) {
  subset <- combined_df %>%
    filter(p == p_val) %>%
    mutate(
      true_IOV  = mapply(function(m, p) { f <- pbinom(0:(m-1), m, p); (4/m) * sum(f * (1-f)) }, m, p),
      true_Skew = mapply(function(m, p) { f <- pbinom(0:(m-1), m, p); (2/m) * sum(f - 1) },     m, p),
      true_C1   = 0
    )

  unique_groups <- unique(subset[, c("n", "pi", "type")])
  unique_groups <- unique_groups[order(
    unique_groups$n, -unique_groups$pi,
    factor(unique_groups$type, levels = c("Simulated", "Asymptotic"))
  ), ]

  x_positions   <- setNames(
    1:nrow(unique_groups),
    paste(unique_groups$n, unique_groups$pi, unique_groups$type, sep = "_")
  )
  subset$x_pos  <- x_positions[paste(subset$n, subset$pi, subset$type, sep = "_")]
  subset$x_pos  <- as.numeric(as.character(subset$x_pos))

  return(subset)
}

#' Compute group center x-positions and labels for the x-axis
#' @param subset  Data frame prepared by prepare_subset()
make_x_axis <- function(subset) {
  group_labels  <- unique(subset[, c("n", "pi")])
  group_labels  <- group_labels[order(group_labels$n, -group_labels$pi), ]
  group_centers <- sapply(1:nrow(group_labels), function(i) {
    mean(subset$x_pos[subset$n == group_labels$n[i] & subset$pi == group_labels$pi[i]])
  })
  x_labels <- paste0("n = ", group_labels$n, "\nπ = ", group_labels$pi)
  list(breaks = group_centers, labels = x_labels)
}

#' Background shading rectangles for n-groups
background_rects <- function() {
  annotate("rect",
    xmin  = c(0.5, 2.5, 4.5, 6.5, 8.5),
    xmax  = c(2.5, 4.5, 6.5, 8.5, 10.5),
    ymin  = -Inf, ymax = Inf,
    alpha = 0.05, fill = "gray90"
  )
}



# ==============================================================================
# Plot 1: Convergence at different rates
# ==============================================================================

scaling_df <- bind_rows(
  sim_data %>% dplyr::mutate(scaled = diff,        type = "IOV_hat - IOV"),
  sim_data %>% dplyr::mutate(scaled = scaled_CLT,  type = "sqrt(n) · (IOV_hat - IOV)"),
  sim_data %>% dplyr::mutate(scaled = scaled_bias, type = "n · (IOV_hat - IOV)")
) %>%
  dplyr::mutate(type = factor(type, levels = c(
    "IOV_hat - IOV",
    "sqrt(n) · (IOV_hat - IOV)",
    "n · (IOV_hat - IOV)"
  )))

mean_lines <- scaling_df %>%
  dplyr::group_by(type, n) %>%
  dplyr::summarise(mean_scaled = mean(scaled, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(type = factor(type, levels = c(
    "IOV_hat - IOV",
    "sqrt(n) · (IOV_hat - IOV)",
    "n · (IOV_hat - IOV)"
  )))

facet_labels <- c(
  "IOV_hat - IOV"             = "hat(IOV) - IOV ~ (Order ~ 1)",
  "sqrt(n) · (IOV_hat - IOV)" = "sqrt(n) ~ (hat(IOV) - IOV) ~ (Order ~ n^{-1/2})",
  "n · (IOV_hat - IOV)"       = "n ~ (hat(IOV) - IOV) ~ (Order ~ n^{-1})"
)

ggplot(scaling_df, aes(x = scaled, fill = factor(n), color = factor(n))) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_vline(
    data = mean_lines,
    aes(xintercept = mean_scaled, color = factor(n)),
    linetype = "solid", linewidth = 0.8
  ) +
  facet_wrap(~ type, scales = "free", nrow = 1, labeller = as_labeller(facet_labels, label_parsed)) +
  coord_cartesian(ylim = c(0, NA), expand = FALSE) +
  labs(fill = "n", color = "n", x = NULL, title = "Convergence at different rates") +
  theme_minimal() +
  theme(aspect.ratio = 0.9)
# ==============================================================================
# Plot 2: Bias als Funktion von 1/n mit theoretischer Kurve
# ==============================================================================

# Theoretischer Bias: -4/(m*n) * sum(diag(Sigma))
Sigma_raw <- Sigma_Star(m = 10, p = 0.3, r = 0.2, pi = 0.75)
theoretical_bias <- function(n) -(4 / (10 * n)) * sum(diag(Sigma_raw))
theory_df <- data.frame(
  n    = seq(50, 5000, by = 10),
  diff = sapply(seq(50, 5000, by = 10), theoretical_bias)
)

mean_df <- aggregate(diff ~ n, data = sim_data_big, mean)

ggplot(mean_df, aes(x = n, y = diff)) +
  geom_line(data = theory_df, aes(x = n, y = diff),
            color = "steelblue", linetype = "dashed", linewidth = 0.8) +
  geom_point(size = 3) +
  geom_line() +
  scale_x_continuous(
    sec.axis = sec_axis(~ 1 / ., name = "1/n",
                        breaks = 1 / c(50, 1000, 5000),
                        labels = c("1/50", "1/1000", "1/5000"))
  ) +
labs(
  fill     = "n", color = "n",
  x        = NULL,
  title    = "Convergence at different rates",
  subtitle = expression(atop(
  BinAR(1) ~ ":" ~ m == 10 ~ "," ~ p == 0.3 ~ "," ~ r == 0.2 ~ "," ~ pi == 0.75,
  "Dashed line: theoretical bias  -4/(mn) · tr(Σ)"
  ))
) +
  theme_minimal()


# ==============================================================================
# Plot 4: Estimated vs. True CDF (m = 3, p = 0.20, r = 0.35, pi = 0.75)
# ==============================================================================

# True CDF
true_cdf <- pbinom(0:2, size = 3, prob = 0.20)
true_df  <- data.frame(category = paste0("f_", 0:2), true_cdf = true_cdf)

# Estimated CDFs from simulation
i_cdf  <- which(scenarios$m == 3 & scenarios$p == 0.20 &
                scenarios$r == 0.35 & scenarios$pi == 0.75)
est_df <- cdf_df %>%
  filter(pi == 0.75) %>%
  pivot_longer(cols = starts_with("f_"), names_to = "category", values_to = "value") %>%
  pivot_wider(names_from = statistic, values_from = value) %>%
  mutate(n = factor(n, levels = sort(unique(n))))

# Asymptotic standard deviations
Sigma_cdf <- Sigma_Star(m = 3, p = 0.20, r = 0.35, pi = 0.75)
asymp_cdf_df <- do.call(rbind, lapply(unique(est_df$n), function(ni) {
  ni_num <- as.numeric(as.character(ni))
  data.frame(
    category = paste0("f_", 0:2),
    n        = factor(ni, levels = levels(est_df$n)),
    asymp_sd = sqrt(diag(Sigma_cdf) / ni_num),
    true_cdf = true_cdf
  )
}))

# Combine for plotting
plot_cdf_df <- bind_rows(
  est_df %>%
    mutate(type = "Simulated", ymin = mean - sd, ymax = mean + sd),
  asymp_cdf_df %>%
    dplyr::rename(mean = true_cdf) %>%
    mutate(type = "Asymptotic", ymin = mean - asymp_sd, ymax = mean + asymp_sd)
) %>%
  dplyr::mutate(type = factor(type, levels = c("Simulated", "Asymptotic")))

ggplot(plot_cdf_df, aes(x = category, color = n, linetype = type)) +
  geom_bar(
    data = true_df, aes(x = category, y = true_cdf),
    inherit.aes = FALSE,
    stat = "identity", fill = "grey85", color = "grey60", width = 0.7
  ) +
  geom_linerange(
    aes(ymin = ymin, ymax = ymax),
    position = position_dodge(width = 0.7), linewidth = 0.8
  ) +
  geom_point(
    aes(y = mean),
    position = position_dodge(width = 0.7), size = 2
  ) +
  scale_color_grey(name = "n", start = 0.7, end = 0.1) +
  scale_linetype_manual(
    name   = "Std. Deviation",
    values = c("Simulated" = "solid", "Asymptotic" = "dotted")
  ) +
  guides(linetype = guide_legend(override.aes = list(linewidth = 1, size = 0, color = "black"))) +
  scale_y_continuous(limits = c(0, 1), name = "Cumulative Probability") +
  scale_x_discrete(name = "CDF Component") +
  labs(
    title    = "Estimated vs. True CDF",
    subtitle = expression(m == 3 ~ "," ~ p == 0.20 ~ "," ~ r == 0.35 ~ "," ~ pi == 0.75)
  ) +
  theme_minimal()


# ==============================================================================
# Plots 5–7: IOV / Skew / Cohen's K — Asymptotic vs. Simulated
# ==============================================================================

combined_df <- rbind(asymp_df, sim_df)

# Change p_val to 0.45 for the second scenario
p_val  <- 0.20
subset <- prepare_subset(combined_df, p_val)
x_axis <- make_x_axis(subset)

t_IOV  <- subset$true_IOV[1]
t_Skew <- subset$true_Skew[1]
t_C1   <- subset$true_C1[1]
m_val  <- subset$m[1]
r_val  <- subset$r[1]

shared_layers <- list(
  scale_x_continuous(breaks = x_axis$breaks, labels = x_axis$labels,
                     expand = expansion(mult = 0.1)),
  scale_color_manual(values = c("Asymptotic" = "#E41A1C", "Simulated" = "#377EB8"),
                     name = "Method"),
  scale_shape_manual(values = c("1" = 16, "0.75" = 17), name = expression(pi)),
  background_rects(),
  thesis_theme()
)

# IOV
ggplot(subset, aes(x = x_pos, y = mean_IOV, color = type, shape = factor(pi))) +
  geom_hline(yintercept = t_IOV, color = "darkgreen", linetype = "dashed",
             linewidth = 1, alpha = 0.7) +
  geom_point(size = 3.5, position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin = lower_IOV, ymax = upper_IOV), width = 0.15,
                linewidth = 0.8, position = position_dodge(width = 0.2)) +
  geom_line(aes(group = interaction(n, pi)), color = "gray50", linetype = "dashed",
            alpha = 0.5, position = position_dodge(width = 0.2)) +
  coord_cartesian(ylim = c(t_IOV - 0.125, t_IOV + 0.075)) +
  labs(
    title    = "Asymptotic vs. Simulated IOV",
    subtitle = sprintf("BinAR(1): m = %d, p = %.2f, r = %.2f", m_val, p_val, r_val),
    x = "Scenario (n and π)", y = "IOV"
  ) +
  shared_layers

# Skew
ggplot(subset, aes(x = x_pos, y = mean_Skew, color = type, shape = factor(pi))) +
  geom_hline(yintercept = t_Skew, color = "darkgreen", linetype = "dashed",
             linewidth = 1, alpha = 0.7) +
  geom_point(size = 3.5, position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin = lower_Skew, ymax = upper_Skew), width = 0.15,
                linewidth = 0.8, position = position_dodge(width = 0.2)) +
  geom_line(aes(group = interaction(n, pi)), color = "gray50", linetype = "dashed",
            alpha = 0.5, position = position_dodge(width = 0.2)) +
  coord_cartesian(ylim = c(t_Skew - 0.15, t_Skew + 0.15)) +
  labs(
    title    = "Asymptotic vs. Simulated Skew",
    subtitle = sprintf("BinAR(1): m = %d, p = %.2f, r = %.2f", m_val, p_val, r_val),
    x = "Scenario (n and π)", y = "Skew"
  ) +
  shared_layers

# Cohen's K lag(1)
ggplot(subset, aes(x = x_pos, y = mean_C, color = type, shape = factor(pi))) +
  geom_hline(yintercept = t_C1, color = "darkgreen", linetype = "dashed",
             linewidth = 1, alpha = 0.7) +
  geom_point(size = 3.5, position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin = lower_C, ymax = upper_C), width = 0.15,
                linewidth = 0.8, position = position_dodge(width = 0.2)) +
  geom_line(aes(group = interaction(n, pi)), color = "gray50", linetype = "dashed",
            alpha = 0.5, position = position_dodge(width = 0.2)) +
  coord_cartesian(ylim = c(-0.3, 0.3)) +
  labs(
    title    = "Asymptotic vs. Simulated Cohen's K lag(1)",
    subtitle = sprintf("BinAR(1): m = %d, p = %.2f, r = %.2f", m_val, p_val, r_val),
    x = "Scenario (n and π)", y = "Cohen's K lag(1)"
  ) +
  shared_layers
