library(dplyr)
library(plyr)
library(tidyr)



#' Generate BinAR(1) process
generate_binar1 <- function(n, m, p, r) {
  # n: Länge der Zeitreihe
  # m: Parameter der Binomialverteilung
  # p: Erwartungswert der stationären Verteilung = E[I_t]
  # r: Thinning-Parameter (muss in [0, 1] sein)
  alpha <- (p * (1-r)) + r
  beta <- p * (1-r)
  # Initialisiere
  I <- numeric(n)
  
  # Startwert aus stationärer Verteilung: Bin(m, p)
  I[1] <- rbinom(1, m, p)
  
  # Erzeuge den Prozess
  for (t in 2:n) {
    
    I[t] <- rbinom(1, I[t-1], alpha) + rbinom(1, m - I[t-1], beta)  
  }
  
  return(I)
}
#'' Generate Amplitude-modulating process
generate_O <- function(n, pi) {
  rbinom(n, 1, pi)
}
# Test with small example
set.seed(123)
test_counts <- generate_binar1(n = 10000, m = 3, p = 0.2, r = 0.35)


#-- Function to compute the empirical marginal cdf 
marginal_probs_e <- function(m, data, n) {
  marg_probs <- numeric(m)
  for (i in 0:m-1) {
    marg_probs[i+1] <- sum(data <= i)*(1/n)
  } 
  return(marg_probs)
}

#-- Function to compute the empirical lag(1) cdf
biv_probs_e <- function(m, data, n, h = 1) {
  biv_probs <- numeric(m)
  for (i in 0:m-1) {
    biv_probs[i+1] <- sum((data[1:(n-h)] <= i)*(data[(h+1):n] <= i))*(1/(n-h))
  } 
  return(biv_probs)
}

#' Simulation function
simulation <- function(n, m, p, r, pi, n_reps = 100) {
  results <- matrix(NA, nrow = n_reps, ncol = 4, dimnames = list(NULL, c("IOV", "Skew", "lag1_Cohen", "lag2_Cohen")))
  for (rep in 1:n_reps) {
  # Generate BinAR(1) process
  count_process <- generate_binar1(n, m, p, r)
  Missing_process <- generate_O(n, pi)
  observed_counts <- count_process[Missing_process ==1]
  CDF <- marginal_probs_e(m, observed_counts, length(observed_counts))
  results[rep, 1] <- (4/m)*(sum(CDF *(1-CDF)))
  results[rep, 2] <- (2/m)*sum(CDF-1)
  results[rep, 3] <- sum(biv_probs_e(m, observed_counts, length(observed_counts), h = 1) - CDF^2)/sum(CDF*(1-CDF))
  results[rep, 4] <- sum(biv_probs_e(m, observed_counts, length(observed_counts), h = 2) - CDF^2)/sum(CDF*(1-CDF))
  }
  return(rbind(colMeans(results, na.rm= T), apply(results, 2, sd)))
}


simulation(n = 100, m = 3, p = 0.2, r = 0.35, pi = 0.8, n_reps = 10)

# Alle Szenarien generieren
scenarios <- expand.grid(
  n = c(50, 100, 250, 500, 1000),
  m = c(3, 10),
  p = c(0.20, 0.45),
  r = c(0, 0.35, 0.50),
  pi = c(1, 0.95, 0.9, 0.85, 0.8, 0.75, 0.5)
)
scenarios <- scenarios[
  (scenarios$m == 3 & scenarios$p == 0.20 & scenarios$r %in% c(0, 0.35)) |
  (scenarios$m == 10 & scenarios$p == 0.45 & scenarios$r %in% c(0, 0.50)),
]

# Index reset
rownames(scenarios) <- NULL

results <- apply(scenarios, 1, function(row) {
  simulation(
    n = row["n"], m = row["m"], p = row["p"], 
    r = row["r"], pi = row["pi"], n_reps = 1000
  )
})

target_scenarios <-  scenarios[
  scenarios$m == 3 & 
  scenarios$p == 0.20 & 
  scenarios$r == 0.35 &
  scenarios$pi %in% c(1, 0.75) &
  scenarios$n %in% c(50, 100, 250, 500, 1000),
]

sim_indices <- which(
  scenarios$m == 3 & 
  scenarios$p == 0.20 & 
  scenarios$r == 0.35 &
  scenarios$pi %in% c(1, 0.75) &
  scenarios$n %in% c(50, 100, 250, 500, 1000)
)

IOV_sim <- results[1:2,sim_indices ]

# Sortiere die Szenarien wie in group_info
sim_scenarios <- scenarios[sim_indices, ]


# Extrahiere die entsprechenden simulierten Ergebnisse
sim_means <- as.numeric(IOV_sim[1, ])  # Erste Zeile: Mittelwerte
sim_sds <- as.numeric(IOV_sim[2, ])    # Zweite Zeile: Standardabweichungen

sim_df <- data.frame(
  n = sim_scenarios$n,
  pi = sim_scenarios$pi,
  m = sim_scenarios$m,
  p = sim_scenarios$p,
  r = sim_scenarios$r,
  type = "Simulated",
  mean = sim_means,
  sd = sim_sds,
  lower = sim_means - sim_sds,
  upper = sim_means + sim_sds
)



selected_indices <- which(scenarios$m == 3 & 
                         scenarios$p == 0.20 & 
                         scenarios$r == 0.35 &
                         scenarios$pi %in% c(1, 0.75) &
                         scenarios$n %in% c(50, 100, 250, 500, 1000))


group_info <- data.frame(
  index = selected_indices,
  n = scenarios$n[selected_indices],
  pi = scenarios$pi[selected_indices]
)
group_info <- group_info[order(group_info$n, -group_info$pi), ]
plot_1 <- plot_1[, as.character(group_info$index)]

x_labels <- paste0("n=", group_info$n, "\nπ=", group_info$pi)
# Correct syntax - remove the duplicate 'levels' argument
# Create the plot data
plot_data <- data.frame(
  scenario = factor(1:ncol(plot_1)),
  label = factor(x_labels, levels = x_labels),  # Preserve order
  n = group_info$n,
  pi = group_info$pi,
  mean = as.numeric(plot_1[1, ]),
  sd = as.numeric(plot_1[2, ])
)




  # Add spacing between different n values
plot_data$group_pos <- as.numeric(factor(plot_data$n, 
                                         levels = sort(unique(plot_data$n))))

# Adjust x-position based on n and pi
plot_data$x_pos <- plot_data$group_pos + 
                   ifelse(plot_data$pi == 1, -0.2, 0.2)

# Create custom x-axis labels (just show n once per group)
n_labels <- sort(unique(plot_data$n))
x_breaks <- seq_along(n_labels)
x_labels <- paste0("n = ", n_labels)

ggplot(plot_data, aes(x = x_pos, y = mean, color = factor(pi))) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                width = 0.1, linewidth = 1) +
  scale_x_continuous(breaks = x_breaks,
                     labels = x_labels,
                     minor_breaks = NULL) +
  labs(title = "Simulation Results",
       subtitle = "π = 1 (left marker) vs π = 0.75 (right marker) for each n",
       x = "Sample Size (n)",
       y = "Value",
       color = "π") +
  scale_color_manual(values = c("1" = "steelblue", "0.75" = "coral"),
                     labels = c("1" = "π = 1", "0.75" = "π = 0.75")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "bottom")
