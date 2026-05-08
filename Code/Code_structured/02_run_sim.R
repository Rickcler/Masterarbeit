# ==============================================================================
# 02_run_sim.R
# Szenarien definieren und Simulationen ausführen
# Erzeugt: sim_df, cdf_df, results_MAR, results_MAR_inv, sim_data, sim_data_2
# ==============================================================================

source("00_setup.R")
source("01_functions.R")
load("Masterarbeit.RData")  # für UNIQUE_N
# ------------------------------------------------------------------------------
# Szenarien: MCAR / serially dependent missingness
# ------------------------------------------------------------------------------

scenarios <- expand.grid(
  n    = UNIQUE_N,
  m    = c(3, 10),
  p    = c(0.20, 0.45),
  r    = c(0, 0.35, 0.50),
  pi   = c(1, 0.75),
  pi_h = c(0, 0.2, 0.75)
)

scenarios <- scenarios[
  (scenarios$m == 3  & scenarios$p == 0.20 & scenarios$r == 0.35 & scenarios$pi_h == 0)   |
  (scenarios$m == 3  & scenarios$p == 0.20 & scenarios$r == 0.35 & scenarios$pi_h == 0.2)  |
  (scenarios$m == 3  & scenarios$p == 0.20 & scenarios$r == 0.35 & scenarios$pi_h == 0.75) |
  (scenarios$m == 10 & scenarios$p == 0.45 & scenarios$r == 0.50 & scenarios$pi_h == 0),
]
rownames(scenarios) <- NULL

# ------------------------------------------------------------------------------
# Szenarien: MAR
# ------------------------------------------------------------------------------

scenarios_MAR <- expand.grid(
  n       = UNIQUE_N,
  m       = c(3, 10),
  p       = c(0.20, 0.45),
  r       = c(0.35, 0.50),
  pi_low  = 0.4,
  pi_high = 0.9
)

scenarios_MAR <- scenarios_MAR[
  (scenarios_MAR$m == 3  & scenarios_MAR$p == 0.20 & scenarios_MAR$r == 0.35) |
  (scenarios_MAR$m == 10 & scenarios_MAR$p == 0.45 & scenarios_MAR$r == 0.50),
]
rownames(scenarios_MAR) <- NULL

# ------------------------------------------------------------------------------
# Hauptsimulationen
# ------------------------------------------------------------------------------

set.seed(SEED)
results <- apply(scenarios, 1, function(row) {
  simulation(
    n    = as.numeric(row["n"]),
    m    = as.numeric(row["m"]),
    p    = as.numeric(row["p"]),
    r    = as.numeric(row["r"]),
    pi   = as.numeric(row["pi"]),
    pi_h = as.numeric(row["pi_h"])
  )
})

# ------------------------------------------------------------------------------
# MAR-Simulationen
# ------------------------------------------------------------------------------

set.seed(SEED)
results_MAR <- apply(scenarios_MAR, 1, function(row) {
  simulation_MAR(
    n       = as.numeric(row["n"]),
    m       = as.numeric(row["m"]),
    p       = as.numeric(row["p"]),
    r       = as.numeric(row["r"]),
    pi_low  = as.numeric(row["pi_low"]),
    pi_high = as.numeric(row["pi_high"])
  )
})

set.seed(SEED)
results_MAR_inv <- apply(scenarios_MAR, 1, function(row) {
  simulation_MAR_inv(
    n       = as.numeric(row["n"]),
    m       = as.numeric(row["m"]),
    p       = as.numeric(row["p"]),
    r       = as.numeric(row["r"]),
    pi_low  = as.numeric(row["pi_low"]),
    pi_high = as.numeric(row["pi_high"])
  )
})

# ------------------------------------------------------------------------------
# Ergebnisse in Data Frames
# ------------------------------------------------------------------------------

# Zusammenfassung je Szenario
sim_df <- do.call(rbind, lapply(seq_along(results), function(i) {
  s <- results[[i]]$summary
  data.frame(
    scenarios[i, ],
    type       = "Simulation",
    mean_IOV   = s["mean", "IOV"],
    sd_IOV     = s["sd",   "IOV"],
    lower_IOV  = s["mean", "IOV"] - s["sd", "IOV"],
    upper_IOV  = s["mean", "IOV"] + s["sd", "IOV"],
    mean_Skew  = s["mean", "Skew"],
    sd_Skew    = s["sd",   "Skew"],
    lower_Skew = s["mean", "Skew"] - s["sd", "Skew"],
    upper_Skew = s["mean", "Skew"] + s["sd", "Skew"],
    mean_C     = s["mean", "lag1_Cohen"],
    sd_C       = s["sd",   "lag1_Cohen"],
    lower_C    = s["mean", "lag1_Cohen"] - s["sd", "lag1_Cohen"],
    upper_C    = s["mean", "lag1_Cohen"] + s["sd", "lag1_Cohen"],
    mean_C_bc  = s["mean", "lag1_Cohen_bc"],
    sd_C_bc    = s["sd", "lag1_Cohen_bc"],
    lower_C_bc = s["mean", "lag1_Cohen_bc"] - s["sd", "lag1_Cohen_bc"],
    upper_C_bc = s["mean", "lag1_Cohen_bc"] + s["sd", "lag1_Cohen_bc"]
  )
}))

# Geschätzte CDFs (nur m = 3)
cdf_df <- do.call(rbind, lapply(which(scenarios$m == 3), function(i) {
  cdf <- results[[i]]$cdf
  data.frame(
    scenario  = i,
    statistic = c("mean", "sd"),
    scenarios[i, ],
    cdf
  )
}))

# ------------------------------------------------------------------------------
# Bias / CLT Simulationen (für 02_run_bias.R genutzt)
# ------------------------------------------------------------------------------

set.seed(SEED)
sim_50   <- simulation_raw(50,   10, 0.3, 0.2, 0.75, n_reps = 2000)
sim_1000 <- simulation_raw(1000, 10, 0.3, 0.2, 0.75, n_reps = 2000)
sim_5000 <- simulation_raw(5000, 10, 0.3, 0.2, 0.75, n_reps = 2000)

sim_data          <- rbind(sim_50, sim_1000, sim_5000)
true_val          <- true_IOV(m = 10, p = 0.3)
sim_data$true_IOV    <- true_val
sim_data$diff        <- sim_data$IOV - true_val
sim_data$scaled_CLT  <- sqrt(sim_data$n) * sim_data$diff
sim_data$scaled_bias <- sim_data$n       * sim_data$diff

# Feinere n-Gitter für Bias-Kurve
set.seed(SEED)
sim_data_2 <- do.call(rbind, lapply(
  c(2, 3, 5, 10, 15, 25, 37, 50, 75, 100, 150, 200, 375, 500, 750, 1000),
  function(n_val) simulation_raw(n_val, 10, 0.3, 0.2, 0.75, n_reps  = 50000)
))

sim_data_2$true_IOV    <- true_val
sim_data_2$diff        <- sim_data_2$IOV - true_val
sim_data_2$scaled_CLT  <- sqrt(sim_data_2$n) * sim_data_2$diff
sim_data_2$scaled_bias <- sim_data_2$n       * sim_data_2$diff

# ------------------------------------------------------------------------------
# Speichern
# ------------------------------------------------------------------------------
Sigma_raw
save.image("Masterarbeit.RData")


