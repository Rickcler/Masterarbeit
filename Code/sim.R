library(dplyr)
library(plyr)
library(tidyr)
library(ggplot2)


#' Generate BinAR(1) process

generate_binar1 <- function(n, m, p, r) {
  # Initialize
  I <- numeric(n)
  
  # Initial value from stationary distribution: Bin(m, p)
  I[1] <- rbinom(1, m, p)
  
  # Generate the process
  for (t in 2:n) {
    # BinAR(1) evolution: I_t = α∘I_{t-1} + ε_t
    # where α∘ is binomial thinning with parameter r
    # and ε_t ~ Bin(m, p(1-r)/(1-p*r))
    alpha_thinned <- rbinom(1, I[t-1], r)
    epsilon <- rbinom(1, m, p*(1-r)/(1-p*r)) # Ensures that the marginal distribution is Bin(m, p)
    I[t] <- alpha_thinned + epsilon
  }
  
  return(I)
}

# Test with small example
set.seed(123)
test_counts <- generate_binar1(n = 10000, m = 3, p = 0.2, r = 0.35)
print(test_counts)

#' Generate amplitude-modulating process (missing observations)

generate_O <- function(n, pi) {
  rbinom(n, 1, pi)
}

marginal_probs_e <- function(m, data, n) {
  marg_probs <- numeric(m)
  for (i in 0:m-1) {
    marg_probs[i+1] <- sum(data <= i)*(1/n)
  } 
  return(marg_probs)
}
marginal_probs_e(m = 3, data = test_counts, n = length(test_counts))

biv_probs_e <- function(m, data, n, h = 1) {
  biv_probs <- numeric(m)
  for (i in 0:m-1) {
    biv_probs[i+1] <- sum((data[1:(n-h)] <= i)*(data[(h+1):n] <= i))*(1/(n-h))
  } 
  return(biv_probs)
}

biv_probs_e(m = 3, data = test_counts, n = length(test_counts), h = 1)


marginal_probs_T <- function(m, p) {
  cumsum(dbinom(0:(m-1), m, p))
}
marginal_probs_T(m = 3, p = 0.2)

joint_cdf_lag1 <- function(m, p, r) {
  # p* for innovation term
  p_star <- p * (1 - r) / (1 - p * r)
  
  joint_pmf <- matrix(0, nrow = m, ncol = m)

  marginal <- dbinom(0:(m-1), m, p)

  for (a in 0: (m-1)) {
    for (b in 0:(m-1)) {
      trans_prob <- 0
      for (k in 0:min(a, b)) {
        thinning_prob <- choose(a,k) * r^k*(1-r)^(a-k)

        innov_prob <- choose(m - a, b - k) * 
          p_star^(b - k) * 
          (1 - p_star)^(m - a - b + k)

        trans_prob <- trans_prob + thinning_prob * innov_prob
      }
      joint_pmf[a + 1, b + 1] <- marginal[a + 1] * trans_prob
    }
  }
  joint_cdf <- matrix(0, nrow = m, ncol = m)
  for (a in 0:(m-1)) {
    for (b in 0:(m-1)) {
      joint_cdf[a + 1, b + 1] <- sum(joint_pmf[1:(a + 1), 1:(b + 1)])
    }
  }
  return(joint_cdf)
}

joint_cdf_lag1(m = 3, p = 0.2, r = 0.35)

#' Simulation function
simulation <- function(n, m, p, r, pi, n_reps = 100) {
  results <- matrix(NA, nrow = n_reps, ncol = 4, dimnames = list(NULL, c("IOV", "Skew", "lag1_Cohen", "lag2_Cohen")))
  for (rep in 1:n_reps) {
  # Generate BinAR(1) process
  count_process <- generate_binar1(n, m, p, r)
  Missing_process <- generate_O(n, pi)
  observed_counts <- count_process[Missing_process ==1]
  CDF <- marginal_probs(m, observed_counts, length(observed_counts))
  results[rep, 1] <- (4/m)*(sum(CDF *(1-CDF)))
  results[rep, 2] <- (2/m)*sum(CDF-1)
  results[rep, 3] <- sum(biv_probs(m, observed_counts, length(observed_counts), h = 1) - CDF^2)/sum(CDF*(1-CDF))
  results[rep, 4] <- sum(biv_probs(m, observed_counts, length(observed_counts), h = 2) - CDF^2)/sum(CDF*(1-CDF))
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
scenarios
# Filter für Paper-Szenarien
scenarios <- scenarios[
  (scenarios$m == 3 & scenarios$p == 0.20 & scenarios$r %in% c(0, 0.35)) |
  (scenarios$m == 10 & scenarios$p == 0.45 & scenarios$r %in% c(0, 0.50)),
]


# Convert scenarios to a list of lists
scenarios_list <- split(scenarios, seq(nrow(scenarios)))
results <- lapply(scenarios_list, function(params) {
  simulation(
    n = params$n, 
    m = params$m, 
    p = params$p, 
    r = params$r, 
    pi = params$pi, 
    n_reps = 100
  )
})
results <- apply(scenarios, 1, function(row) {
  simulation(
    n = row["n"], m = row["m"], p = row["p"], 
    r = row["r"], pi = row["pi"], n_reps = 100
  )
})

target_scenarios <- scenarios[
  scenarios$m == 3 & 
  scenarios$p == 0.20 & 
  scenarios$r == 0.35 &
  scenarios$pi %in% c(1, 0.75) &
  scenarios$n %in% c(50, 100, 250, 500, 1000)
]
plot_1 <- results[,scenarios$m == 3 & 
  scenarios$p == 0.20 & 
  scenarios$r == 0.35 &
  scenarios$pi %in% c(1, 0.75) &
  scenarios$n %in% c(50, 100, 250, 500, 1000)]


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

get_marginal_cdf  <- function(m,p) {
  cdf <- numeric(m-1)
  for (i in 0:(m-1)) {
    cdf[i+1] <- pbinom(i, m, p)
  }
  return(cdf)
}

get_marginal_cdf(3, 0.2)