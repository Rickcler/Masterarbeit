library(ggplot2)

combined_df <- rbind(asymp_df, sim_df)

# 4. Erstelle eine kombinierte x-Position
# Wir wollen: n → pi → type (Asymptotic/Simulated)
combined_df$group_id <- interaction(
  combined_df$n, 
  combined_df$pi, 
  combined_df$type,
  sep = "_"
)


# Erstelle eindeutige x-Positionen
unique_groups <- unique(combined_df[, c("n", "pi", "type")])
unique_groups <- unique_groups[order(unique_groups$n, -unique_groups$pi, 
                                     factor(unique_groups$type, levels = c("Simulated","Asymptotic" ))), ]

# Weise x-Positionen zu
x_positions <- setNames(1:nrow(unique_groups), 
                       paste(unique_groups$n, unique_groups$pi, unique_groups$type, sep = "_"))

combined_df$x_pos <- x_positions[as.character(paste(combined_df$n, combined_df$pi, combined_df$type, sep = "_"))]

# 5. Erstelle x-Achsen-Beschriftungen
# Gruppiere nach n und pi
group_labels <- unique(combined_df[, c("n", "pi")])
group_labels <- group_labels[order(group_labels$n, -group_labels$pi), ]

# Bestimme die mittleren Positionen für jede n-pi Kombination
group_centers <- sapply(1:nrow(group_labels), function(i) {
  n_val <- group_labels$n[i]
  pi_val <- group_labels$pi[i]
  pos <- combined_df$x_pos[combined_df$n == n_val & combined_df$pi == pi_val]
  mean(pos)
})

# Erstelle Beschriftungen
x_labels <- paste0("n = ", group_labels$n, "\nπ = ", group_labels$pi)

# Marginal CDF berechnen
marginal_cdf <- pbinom(0:2, 3, 0.2)  # Nur i=1,...,m (nicht 0 oder m)
# Für m=3: i=1,2,3

# Wahrer IOV
true_IOV <- (4/3) * sum(marginal_cdf * (1 - marginal_cdf))
true_skew <- (2/m) * sum(marginal_cdf- 1)
# 6.IOV Plot erstellen
ggplot(combined_df, aes(x = x_pos, y = mean_IOV, 
                        color = type, 
                        shape = factor(pi))) +
  # Horizontale Linie für wahren IOV-Wert
  geom_hline(yintercept = true_IOV, 
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
  geom_line(data = combined_df, 
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
    values = c("Asymptotic" = "#E41A1C", "Simulated" = "#377EB8"),
    name = "Method"
  ) +
  
  scale_shape_manual(
    values = c("1" = 16, "0.75" = 17),  # 16 = Kreis, 17 = Dreieck
    name = "π value"
  ) +
  
  # Labels und Titel
  labs(
    title = "Comparison of Asymptotic vs Simulated IOV Results",
    subtitle = "For BinAR(1) process with m=3, p=0.2, r=0.35",
    x = "Scenario (Sample Size n and Missing Probability π)",
    y = "IOV Value"
  ) +
  
  # Y-Achse begrenzen auf 0.35-0.55
  coord_cartesian(ylim = c(0.35, 0.55)) +
  
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
ggplot(combined_df, aes(x = x_pos, y = mean_Skew, 
                        color = type, 
                        shape = factor(pi))) +
  # Horizontale Linie für wahren IOV-Wert
  geom_hline(yintercept = true_skew, 
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
  geom_line(data = combined_df, 
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
    values = c("Asymptotic" = "#E41A1C", "Simulated" = "#377EB8"),
    name = "Method"
  ) +
  
  scale_shape_manual(
    values = c("1" = 16, "0.75" = 17),  # 16 = Kreis, 17 = Dreieck
    name = "π value"
  ) +
  
  # Labels und Titel
  labs(
    title = "Comparison of Asymptotic vs Simulated Skew Results",
    subtitle = "For BinAR(1) process with m=3, p=0.2, r=0.35",
    x = "Scenario (Sample Size n and Missing Probability π)",
    y = "IOV Value"
  ) +
  
  # Y-Achse begrenzen auf 0.35-0.55
  coord_cartesian(ylim = c(-0.55, -0.25)) +
  
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
