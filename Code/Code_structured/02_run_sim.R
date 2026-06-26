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
  r_pi = c(0, 0.2, 0.75)
)

scenarios <- scenarios[
  (scenarios$m == 3  & scenarios$p == 0.20 & scenarios$r == 0.35 & scenarios$r_pi == 0)   |
  (scenarios$m == 3  & scenarios$p == 0.20 & scenarios$r == 0.35 & scenarios$r_pi == 0.2)  |
  (scenarios$m == 3  & scenarios$p == 0.20 & scenarios$r == 0.35 & scenarios$r_pi == 0.75) |
  (scenarios$m == 10 & scenarios$p == 0.45 & scenarios$r == 0.50 & scenarios$r_pi == 0) |
  (scenarios$m == 10 & scenarios$p == 0.45 & scenarios$r == 0.50 & scenarios$r_pi == 0.2)  |
  (scenarios$m == 10 & scenarios$p == 0.45 & scenarios$r == 0.50 & scenarios$r_pi == 0.75), 
]
rownames(scenarios) <- NULL

# ------------------------------------------------------------------------------
# Szenarien: MAR & MNAR
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


# Scenarios
scenarios_mnar <- expand.grid(
  r = c(0.00, 0.15, 0.35, 0.6),
  mechanism = c("increasing", "decreasing"),
  stringsAsFactors = FALSE
)
# Simulation
n_grid_mnar <- c(50, 100, 250, 500, 1000, 2000)

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
    r_pi = as.numeric(row["r_pi"])
  )
})

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


# ------------------------------------------------------------------------------
# MAR-Simulationen
# ------------------------------------------------------------------------------
n_grid  <- c(50, 100, 250, 500, 1000, 2000) 

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

# Figure 4.12 left plot
set.seed(42)

kappa_inc_list <- lapply(n_grid, function(n) {
  vals <- simulation_kappa_HA_MNAR(
    n, m_val, p_val, r_val, pi_val,
    h = h_val,
    n_reps = 1000,
    mechanism = "increasing"
  )
  data.frame(n = n, kappa_hat = vals, mech = "Increasing")
})

kappa_dec_list <- lapply(n_grid, function(n) {
  vals <- simulation_kappa_HA_MNAR(
    n, m_val, p_val, r_val, pi_val,
    h = h_val,
    n_reps = 1000,
    mechanism = "decreasing"
  )
  data.frame(n = n, kappa_hat = vals, mech = "Decreasing")
})

# Figure 4.12 right plot
set.seed(42)
rej_list_mnar <- lapply(1:nrow(scenarios_mnar), function(s) {

  sc <- scenarios_mnar[s, ]

  rates <- sapply(n_grid_mnar, function(n) {

    rejection_rate_MNAR(
      n          = n,
      m          = 3,
      p          = 0.20,
      r          = sc$r,
      pi_low     = 0.5,
      pi_high    = 0.9,
      h          = 1,
      alpha      = 0.05,
      n_reps     = 1000,
      mechanism  = sc$mechanism
    )
  })

  data.frame(
    n = n_grid_mnar,
    rate = rates,
    r = sc$r,
    mechanism = sc$mechanism
  )
})


rej_df_mnar <- do.call(rbind, rej_list_mnar)



# ------------------------------------------------------------------------------
# Bias / CLT Simulationen 
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
set.seed(42)
sim_data_2 <- do.call(rbind, lapply(
  c(2, 3, 5, 10, 15, 25, 37, 50, 75, 100, 150, 200, 375, 500, 750, 1000),
  function(n_val) simulation_raw(n_val, 10, 0.3, 0.2, 0.75, n_reps  = 50000)
))

sim_data_2$true_IOV    <- true_val
sim_data_2$diff        <- sim_data_2$IOV - true_val
sim_data_2$scaled_CLT  <- sqrt(sim_data_2$n) * sim_data_2$diff
sim_data_2$scaled_bias <- sim_data_2$n       * sim_data_2$diff

# ------------------------------------------------------------------------------
# Kappa H_0 Simulationen 
# ------------------------------------------------------------------------------

n_grid_H0  <- c(25, 37, 50, 75, 100,  200, 375, 500, 750, 1000, 2500, 5000)
alpha_test <- 0.05

# Szenarien A und B mit r = 0 (H_0: i.i.d.)
scenarios_H0 <- data.frame(
  m     = c(3,    10),
  p     = c(0.20, 0.45),
  pi    = c(0.75, 0.75),
  pi_h  = c(0,    0),

  label = c("Scenario A  (m = 3, p = 0.20)",
            "Scenario B  (m = 10, p = 0.45)")
)


set.seed(42)
rej_H0_list <- lapply(1:nrow(scenarios_H0), function(s) {
  sc <- scenarios_H0[s, ]

  do.call(rbind, lapply(n_grid_H0, function(n) {

    # Asymptotische Varianz unter H_0 für kritischen Wert
    asymp <- Cohens_asymptotic_iid(
      n,
      sc$pi,
      sc$m,
      pbinom(0:(sc$m - 1), sc$m, sc$p)
    )

    crit <- qnorm(1 - alpha_test / 2) * asymp["sd"]


    # Wir brauchen rohe Werte — simulation() gibt nur mean/sd zurück
    # Daher raw-Replikation inline
    kappa_raw <- replicate(5000, {
      iid_proc <- generate_iid(n, sc$m, sc$p)
      O_proc   <- generate_binar1(n, 1, sc$pi, sc$pi_h)
      obs_iid  <- iid_proc[O_proc == 1]
      CDF_iid  <- marginal_probs_e(sc$m, obs_iid, length(obs_iid))
      pi_est   <- mean(O_proc)
      biv      <- biv_probs_e(sc$m, iid_proc, O_proc, n, h = 1)
      denom    <- sum(CDF_iid * (1 - CDF_iid))
      if (denom > 0 && pi_est > 0) {
        kappa    <- sum(biv - CDF_iid^2) / denom
        kappa_bc <- kappa + 1 / (pi_est * n)
        c(kappa, kappa_bc)
      } else {
        c(NA, NA)
      }
    })

    data.frame(
      n          = n,
      scenario   = sc$label,
      m          = sc$m,
      p          = sc$p,
      pi         = sc$pi,
      rej_uncorr = rejection_from_vals(kappa_raw[1, ], crit),
      rej_bc     = rejection_from_vals(kappa_raw[2, ], crit)
    )
  }))
})

rej_H0_df <- do.call(rbind, rej_H0_list) %>%
  pivot_longer(
    cols      = c(rej_uncorr, rej_bc),
    names_to  = "estimator",
    values_to = "rejection_rate"
  ) %>%
  mutate(
    estimator = factor(
      estimator,
      levels = c("rej_uncorr", "rej_bc"),
      labels = c("Uncorrected", "Bias-corrected")
    ),
    scenario = factor(scenario, levels = scenarios_H0$label)
  )


# ------------------------------------------------------------------------------
# Kappa H_A Simulationen 
# ------------------------------------------------------------------------------

### Confidence Intervals 

# Szenarien: verschiedene n, festes r > 0
n_grid  <- c(50, 100, 250, 500, 1000, 2000)
m_val   <- 3
p_val   <- 0.20
r_val   <- 0.35
pi_val  <- 0.75
h_val   <- 1

true_kappa <- true_kappa_HA(m_val, p_val, r_val, h = h_val)

set.seed(SEED)
kappa_HA_list <- lapply(n_grid, function(n) {
  vals <- simulation_kappa_HA(n, m_val, p_val, r_val, pi_val,
                               h = h_val, n_reps = 1000)
  data.frame(
    n         = n,
    kappa_hat = vals
  )
})


kappa_HA_df <- do.call(rbind, kappa_HA_list) %>%
  mutate(n = factor(n, levels = n_grid))

# Zusammenfassung: Mittelwert und SD pro n
kappa_summary <- kappa_HA_df %>%
  group_by(n) %>%
  summarise(
    mean_kappa = mean(kappa_hat, na.rm = TRUE),
    sd_kappa   = sd(kappa_hat,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(n_num = as.numeric(as.character(n)))

# 95%-Konfidenzband unter H_0 für jedes n
ci_df <- data.frame(
  n_num = n_grid
) %>%
  rowwise() %>%
  mutate(
    var_kappa = asymp_var_kappa_H0(m_val, p_val, pi_val) / n_num,
    sd_H0     = sqrt(var_kappa),
    ci_lower  = -1.96 * sd_H0,
    ci_upper  =  1.96 * sd_H0
  ) %>%
  ungroup()
 
### Rejection Rates

scenarios_rej <- expand.grid(
  r  = c(0.15, 0.35, 0.60),
  pi = c(1.00, 0.75)
) %>%
  mutate(
    label = paste0(
      "r = ", r, ", π = ", pi
    ),
    # Linientyp nach pi
    lty = ifelse(pi == 1.00, "solid", "dashed"),
    # Farbe nach r
    col = case_when(
      r == 0.15 ~ "#5B8DB8",
      r == 0.35 ~ "#2C3E6B",
      r == 0.60 ~ "#C0392B"
    )
  )

n_grid_rej <- c(50, 100, 250, 500, 1000, 2000)
m_val      <- 3
p_val      <- 0.20
h_val      <- 1

set.seed(42)
rej_list <- lapply(1:nrow(scenarios_rej), function(s) {
  sc <- scenarios_rej[s, ]
  message("Scenario: r = ", sc$r, ", pi = ", sc$pi)

  rates <- sapply(n_grid_rej, function(n) {
    rejection_rate(
      n      = n,
      m      = m_val,
      p      = p_val,
      r      = sc$r,
      pi     = sc$pi,
      h      = h_val,
      alpha  = 0.05,
      n_reps = 1000
    )
  })

  data.frame(
    n     = n_grid_rej,
    rate  = rates,
    r     = sc$r,
    pi    = sc$pi,
    label = sc$label,
    col   = sc$col,
    lty   = sc$lty
  )
})

rej_df <- do.call(rbind, rej_list) %>%
  mutate(
    label = factor(label, levels = scenarios_rej$label),
    r_fac = factor(paste0("r = ", r),
                   levels = c("r = 0.15", "r = 0.35", "r = 0.6")),
    pi_fac = factor(paste0("π = ", pi),
                    levels = c("π = 1", "π = 0.75"))
  )

# ------------------------------------------------------------------------------
# Speichern
# ------------------------------------------------------------------------------
save.image("Masterarbeit.RData")


