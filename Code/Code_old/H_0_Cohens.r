# ==============================================================================
# Rejection Rate unter H_0 — Level Control
# Nutzt simulation() direkt mit r = 0
# ==============================================================================

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

