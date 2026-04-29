library(ggplot2)
library(dplyr)
library(plyr)
library(tidyr)



combined_df <- rbind(asymp_df, sim_df)
head(combined_df)
# 4. Erstelle eine kombinierte x-Position
# Wir wollen: n → pi → type (Asymptotic/Simulated)
combined_df$group_id <- interaction(
  combined_df$p,
  combined_df$n, 
  combined_df$pi, 
  combined_df$type,
  sep = "_"
)

# Assuming your data has different p values
# First, let's add the true values for each p
combined_df <- combined_df %>%
  mutate(
    true_IOV = apply(combined_df, 1, function(row) {
      m_val <- as.numeric(row["m"])
      p_val <- as.numeric(row["p"])
      marginal <- pbinom(0:(m_val-1), m_val, p_val)
      (4/m_val) * sum(marginal * (1 - marginal))
    }),
    true_Skew = apply(combined_df, 1, function(row) {
      m_val <- as.numeric(row["m"])
      p_val <- as.numeric(row["p"])
      marginal <- pbinom(0:(m_val-1), m_val, p_val)
      (2/m_val) * sum(marginal - 1)
    }),
    true_C1 = 0
  )



# ==============================================================================
# Asymptotic vs Simulated - no serial dependence in O_t
# ==============================================================================

c = 0.2 # Ändere {0.2, 0.45} 

subset_pi <- combined_df %>% filter(p == c, pi_h == 0, pi == 0.75)


unique_groups <- unique(subset_pi[, c("n", "pi", "type")])
unique_groups <- unique_groups[order(unique_groups$n, -unique_groups$pi, 
                                    factor(unique_groups$type, levels = c("Simulation","Asymptotic" ))), ]

# Weise x-Positionen zu
x_positions <- setNames(1:nrow(unique_groups), 
                      paste(unique_groups$n, unique_groups$pi, unique_groups$type, sep = "_"))

subset_pi$x_pos <- x_positions[as.character(paste(subset_pi$n, subset_pi$pi, subset_pi$type, sep = "_"))]

# 5. Erstelle x-Achsen-Beschriftungen
# Gruppiere nach n und pi
group_labels <- unique(subset_pi[, c("n", "pi")])
group_labels <- group_labels[order(group_labels$n, -group_labels$pi), ]

# Bestimme die mittleren Positionen für jede n-pi Kombination
group_centers <- sapply(1:nrow(group_labels), function(i) {
  n_val <- group_labels$n[i]
  pi_val <- group_labels$pi[i]
  pos <- subset_pi$x_pos[subset_pi$n == n_val & subset_pi$pi == pi_val]
  mean(pos)
})

# Erstelle Beschriftungen
x_labels <- paste0("n = ", group_labels$n, "\nπ = ", group_labels$pi)

subset_pi$x_pos <- as.numeric(as.character(subset_pi$x_pos))

t_IOV <- subset_pi$true_IOV[1]          # Correct way: extract first value
t_Skew <- subset_pi$true_Skew[1]        # Correct way: extract first value
t_C1 <- subset_pi$true_C1[1]        # Correct way: extract first value
m_val <- subset_pi$m[1]                 # Correct way: extract first value
p_val <- subset_pi$p[1]                 # Correct way: extract first value
r_val <- subset_pi$r[1]                 # Correct way: extract first value
subset_pi



ggplot(subset_pi, aes(x = x_pos, y = mean_IOV, 
                        color = type, 
                        shape = factor(pi))) +
  # Horizontale Linie für wahren IOV-Wert
  geom_hline(yintercept = t_IOV, 
            color = "darkgreen", 
            linetype = "dashed",
            linewidth = 1,
            alpha = 0.7) +
  # Punkte für Mittelwerte
  geom_point(size = 3.5, position = position_dodge(width = 0.2)) +
  
  # Fehlerbalken
  geom_errorbar(aes(ymin = lower_IOV, ymax = upper_IOV),
                width = 0.15, 
                linewidth = 0.8,
                position = position_dodge(width = 0.2)) +
  
  # Verbindungslinien zwischen Asymptotic und Simulated für gleiche n,pi
  geom_line(data = subset_pi, 
            aes(group = interaction(n, pi)),
            color = "gray50", 
            linetype = "dashed",
            alpha = 0.5,
            position = position_dodge(width = 0.2)) +
  
  # X-Achse anpassen
  scale_x_continuous(
    breaks = group_centers,
    labels = x_labels,
    expand = expansion(mult = 0.1)
  ) +
  
  # Farben und Symbole
  scale_color_manual(
    values = c("Asymptotic" = "#E41A1C", "Simulation" = "#377EB8"),
    name = "Method"
  ) +
  
  scale_shape_manual(
    values = c("1" = 16, "0.75" = 17),  # 16 = Kreis, 17 = Dreieck
    name = "π value"
  ) +
  
  # Labels und Titel
  labs(
    title = "Comparison of Asymptotic vs Simulated IOV Results",
    subtitle = sprintf("For BinAR(1) process with m = %d, p = %.2f, r = %.2f", m_val, p_val, r_val),
    x = "Scenario (Sample Size n and Missing Probability π)",
    y = "IOV Value"
  ) +
  
  # Y-Achse begrenzen auf 0.35-0.55
  coord_cartesian(ylim = c(t_IOV - 0.125,t_IOV + 0.075)) +
  
  # Theme
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.spacing.y = unit(0.2, "cm"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
  ) +
  
  # Optional: Hintergrund für Gruppen
  annotate("rect", 
          xmin = c(0.5, 2.5, 4.5, 6.5, 8.5),
          xmax = c(2.5, 4.5, 6.5, 8.5, 10.5),
          ymin = -Inf, ymax = Inf,
          alpha = 0.05, fill = "gray90")


# 8. Skew Plot erstellen
ggplot(subset_pi, aes(x = x_pos, y = mean_Skew, 
                        color = type, 
                        shape = factor(pi))) +
  # Horizontale Linie für wahren IOV-Wert
  geom_hline(yintercept = t_Skew, 
            color = "darkgreen", 
            linetype = "dashed",
            linewidth = 1,
            alpha = 0.7) +
  # Punkte für Mittelwerte
  geom_point(size = 3.5, position = position_dodge(width = 0.2)) +
  
  # Fehlerbalken
  geom_errorbar(aes(ymin = lower_Skew, ymax = upper_Skew),
                width = 0.15, 
                linewidth = 0.8,
                position = position_dodge(width = 0.2)) +
  
  # Verbindungslinien zwischen Asymptotic und Simulated für gleiche n,pi
  geom_line(data =  subset_pi, 
            aes(group = interaction(n, pi)),
            color = "gray50", 
            linetype = "dashed",
            alpha = 0.5,
            position = position_dodge(width = 0.2)) +
  
  # X-Achse anpassen
  scale_x_continuous(
    breaks = group_centers,
    labels = x_labels,
    expand = expansion(mult = 0.1)
  ) +
  
  # Farben und Symbole
  scale_color_manual(
    values = c("Asymptotic" = "#E41A1C", "Simulation" = "#377EB8"),
    name = "Method"
  ) +
  
  scale_shape_manual(
    values = c("1" = 16, "0.75" = 17),  # 16 = Kreis, 17 = Dreieck
    name = "π value"
  ) +
  
  # Labels und Titel
  labs(
    title = "Comparison of Asymptotic vs Simulated Skew Results",
    subtitle = sprintf("For BinAR(1) process with m = %d, p = %.2f, r = %.2f", m_val, p_val, r_val),
    x = "Scenario (Sample Size n and Missing Probability π)",
    y = "Skew Value"
  ) +
  
  # Y-Achse begrenzen auf 0.35-0.55
  coord_cartesian(ylim = c(t_Skew -0.15, t_Skew+0.15)) +
  
  # Theme
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.spacing.y = unit(0.2, "cm"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
  ) +
  
  # Optional: Hintergrund für Gruppen
  annotate("rect", 
          xmin = c(0.5, 2.5, 4.5, 6.5, 8.5),
          xmax = c(2.5, 4.5, 6.5, 8.5, 10.5),
          ymin = -Inf, ymax = Inf,
          alpha = 0.05, fill = "gray90")

# 9. Cohens(K) Plot erstellen
ggplot(subset_pi, aes(x = x_pos, y = mean_C, 
                        color = type, 
                        shape = factor(pi))) +
  # Horizontale Linie für wahren IOV-Wert
  geom_hline(yintercept = t_C1, 
            color = "darkgreen", 
            linetype = "dashed",
            linewidth = 1,
            alpha = 0.7) +
  # Punkte für Mittelwerte
  geom_point(size = 3.5, position = position_dodge(width = 0.2)) +
  
  # Fehlerbalken
  geom_errorbar(aes(ymin = lower_C, ymax = upper_C),
                width = 0.15, 
                linewidth = 0.8,
                position = position_dodge(width = 0.2)) +
  
  # Verbindungslinien zwischen Asymptotic und Simulated für gleiche n,pi
  geom_line(data = subset_pi, 
            aes(group = interaction(n, pi)),
            color = "gray50", 
            linetype = "dashed",
            alpha = 0.5,
            position = position_dodge(width = 0.2)) +
  
  # X-Achse anpassen
  scale_x_continuous(
    breaks = group_centers,
    labels = x_labels,
    expand = expansion(mult = 0.1)
  ) +
  
  # Farben und Symbole
  scale_color_manual(
    values = c("Asymptotic" = "#E41A1C", "Simulation" = "#377EB8"),
    name = "Method"
  ) +
  
  scale_shape_manual(
    values = c("1" = 16, "0.75" = 17),  # 16 = Kreis, 17 = Dreieck
    name = "π value"
  ) +
  
  # Labels und Titel
  labs(
    title = "Comparison of Asymptotic vs Simulated Cohens k lag(1) Results",
    subtitle = sprintf("For BinAR(1) process with m = %d, p = %.2f, r = %.2f", m_val, p_val, r_val),
    x = "Scenario (Sample Size n and Missing Probability π)",
    y = "Cohens K lag(1) Value"
  ) +
  
  # Y-Achse begrenzen auf 0.35-0.55
  coord_cartesian(ylim = c(-0.3, 0.3)) +
  
  # Theme
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.spacing.y = unit(0.2, "cm"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
  ) +
  
  # Optional: Hintergrund für Gruppen
  annotate("rect", 
          xmin = c(0.5, 2.5, 4.5, 6.5, 8.5),
          xmax = c(2.5, 4.5, 6.5, 8.5, 10.5),
          ymin = -Inf, ymax = Inf,
          alpha = 0.05, fill = "gray90")
        

# ------------------------------------------------------------
# Aymptotic vs Simulated - With serial dependence in O_t
# ------------------------------------------------------------


subset_pi_h <- combined_df %>% filter(pi_h == 0.75 | pi_h == 0.2, pi == 0.75)

unique_groups <- unique(subset_pi_h[, c("n", "pi_h", "type")])
unique_groups <- unique_groups[order(unique_groups$n, unique_groups$pi_h, 
                                    factor(unique_groups$type, levels = c("Simulated","Asymptotic" ))), ]
x_positions <- setNames(1:nrow(unique_groups), 
                      paste(unique_groups$n, unique_groups$pi_h, unique_groups$type, sep = "_"))

subset_pi_h$x_pos <- x_positions[as.character(paste(subset_pi_h$n, subset_pi_h$pi_h, subset_pi_h$type, sep = "_"))]

# 5. Erstelle x-Achsen-Beschriftungen
# Gruppiere nach n und pi
group_labels <- unique(subset_pi_h[, c("n", "pi_h")])
group_labels <- group_labels[order(group_labels$n, group_labels$pi_h), ]


# Bestimme die mittleren Positionen für jede n-pi Kombination
group_centers <- sapply(1:nrow(group_labels), function(i) {
  n_val <- group_labels$n[i]
  pi_val <- group_labels$pi_h[i]
  pos <- subset_pi_h$x_pos[subset_pi_h$n == n_val & subset_pi_h$pi_h == pi_val]
  mean(pos)
})

# Erstelle Beschriftungen
x_labels <- paste0("n = ", group_labels$n, "\nπ(h) = ", group_labels$pi_h)

subset_pi_h$x_pos <- as.numeric(as.character(subset_pi_h$x_pos))

t_IOV <- subset_pi_h$true_IOV[1]          # Correct way: extract first value
t_Skew <- subset_pi_h$true_Skew[1]        # Correct way: extract first value
t_C1 <- subset_pi_h$true_C1[1]        # Correct way: extract first value
m_val <- subset_pi_h$m[1]                 # Correct way: extract first value
p_val <- subset_pi_h$p[1]                 # Correct way: extract first value
r_val <- subset_pi_h$r[1]                 # Correct way: extract first value
subset_pi_h
ggplot(subset_pi_h, aes(x = x_pos, y = mean_IOV, 
                        color = type, 
                        shape = factor(pi_h))) +
  # Horizontale Linie für wahren IOV-Wert
  geom_hline(yintercept = t_IOV, 
            color = "darkgreen", 
            linetype = "dashed",
            linewidth = 1,
            alpha = 0.7) +
  # Punkte für Mittelwerte
  geom_point(size = 3.5, position = position_dodge(width = 0.2)) +
  
  # Fehlerbalken
  geom_errorbar(aes(ymin = lower_IOV, ymax = upper_IOV),
                width = 0.15, 
                linewidth = 0.8,
                position = position_dodge(width = 0.2)) +
  
  # Verbindungslinien zwischen Asymptotic und Simulated für gleiche n,pi
  geom_line(data = subset_pi_h, 
            aes(group = interaction(n, pi_h)),
            color = "gray50", 
            linetype = "dashed",
            alpha = 0.5,
            position = position_dodge(width = 0.2)) +
  
  # X-Achse anpassen
  scale_x_continuous(
    breaks = group_centers,
    labels = x_labels,
    expand = expansion(mult = 0.1)
  ) +
  
  # Farben und Symbole
  scale_color_manual(
    values = c("Asymptotic" = "#E41A1C", "Simulation" = "#377EB8"),
    name = "Method"
  ) +
  
  scale_shape_manual(
    values = c("0.2" = 16, "0.75" = 17),  # 16 = Kreis, 17 = Dreieck
    name = "π(h) value"
  ) +
  
  # Labels und Titel
  labs(
    title = "Comparison of Asymptotic vs Simulated IOV Results \n for serially dependent Missingness",
    subtitle = sprintf("For BinAR(1) process with m = %d, p = %.2f, r = %.2f, π = 0.75",  m_val, p_val, r_val),
    x = "Scenario (Sample Size n and Missing π(h))",
    y = "IOV Value"
  ) +
  
  # Y-Achse begrenzen auf 0.35-0.55
  coord_cartesian(ylim = c(t_IOV - 0.125,t_IOV + 0.075)) +
  
  # Theme
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.spacing.y = unit(0.2, "cm"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
  ) +
  
  # Optional: Hintergrund für Gruppen
  annotate("rect", 
          xmin = c(0.5, 2.5, 4.5, 6.5, 8.5),
          xmax = c(2.5, 4.5, 6.5, 8.5, 10.5),
          ymin = -Inf, ymax = Inf,
          alpha = 0.05, fill = "gray90")


# ==============================================================================
# Wahre Marginalverteilungen für beide Szenarien
# ==============================================================================

# Szenario A: m=3, p=0.20
scen_A <- data.frame(
  scenario  = "Scenario A\nm = 3, p = 0.20, r = 0.35",
  category  = 0:3,
  pmf       = dbinom(0:3, size = 3, prob = 0.20),
  cdf       = pbinom(0:3, size = 3, prob = 0.20)
)

# Szenario B: m=10, p=0.45
scen_B <- data.frame(
  scenario  = "Scenario B\nm = 10, p = 0.45, r = 0.50",
  category  = 0:10,
  pmf       = dbinom(0:10, size = 10, prob = 0.45),
  cdf       = pbinom(0:10, size = 10, prob = 0.45)
)

plot_df <- bind_rows(scen_A, scen_B) %>%
  mutate(scenario = factor(scenario,
                           levels = c(
                             "Scenario A\nm = 3, p = 0.20, r = 0.35",
                             "Scenario B\nm = 10, p = 0.45, r = 0.50"
                           )),
         category = factor(category))


# PMF
p_pmf <- ggplot(plot_df,
                aes(x = category, y = pmf)) +
  geom_col(fill = "grey70", color = "grey40", width = 0.6) +
  geom_text(aes(label = round(pmf, 3)),
            vjust = -0.4, size = 2.8, color = "grey30") +
  facet_wrap(~ scenario, scales = "free_x") +
  scale_y_continuous(limits = c(0, 0.60),
                     name   = "Probability") +
  scale_x_discrete(name = "Category") +
  labs(title = "Marginal PMF") +
  theme_minimal() +
  theme(
    strip.text       = element_text(size = 10, face = "bold"),
    panel.grid.major.x = element_blank(),
    plot.title       = element_text(size = 11)
  )

# CDF
p_cdf <- ggplot(plot_df,
                aes(x = category, y = cdf)) +
  geom_col(fill = "grey70", color = "grey40", width = 0.6) +
  geom_text(aes(label = round(cdf, 3)),
            vjust = -0.4, size = 2.8, color = "grey30") +
  geom_hline(yintercept = 0.5,
             linetype = "dashed", color = "steelblue",
             linewidth = 0.6) +
  annotate("text", x = 0.6, y = 0.52,
           label = "0.5", color = "steelblue", size = 3) +
  facet_wrap(~ scenario, scales = "free_x") +
  scale_y_continuous(limits = c(0, 1.05),
                     name   = "Cumulative Probability") +
  scale_x_discrete(name = "Category") +
  labs(title = "Marginal CDF") +
  theme_minimal() +
  theme(
    strip.text         = element_text(size = 10, face = "bold"),
    panel.grid.major.x = element_blank(),
    plot.title         = element_text(size = 11)
  )

# Kombiniert
library(patchwork)

p_pmf / p_cdf +
  plot_annotation(
    title    = "Marginal distributions of the two simulation scenarios",
    subtitle = "Top: PMF    Bottom: CDF    Dashed line: 0.5",
    theme    = theme(
      plot.title    = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "grey40")
    )
  )

load("Masterarbeit.RData")