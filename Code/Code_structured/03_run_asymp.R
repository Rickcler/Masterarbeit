# ==============================================================================
# 03_run_asymp.R
# Asymptotische Ergebnisse berechnen
# Erzeugt: asymp_df, Sigma_raw (für Bias-Kurve)
# Voraussetzung: 01_functions.R muss geladen sein
# ==============================================================================

source("00_setup.R")
source("01_functions.R")

# ------------------------------------------------------------------------------
# Parameterkoeffizientengitter
# ------------------------------------------------------------------------------

unique_coeffs <- data.frame(rbind(
  c(m = 3,  p = 0.20, r = 0.35, pi = 1,    pi_h = 0),
  c(m = 3,  p = 0.20, r = 0.35, pi = 0.75, pi_h = 0),
  c(m = 3,  p = 0.20, r = 0.35, pi = 1,    pi_h = 0.2),
  c(m = 3,  p = 0.20, r = 0.35, pi = 0.75, pi_h = 0.2),
  c(m = 3,  p = 0.20, r = 0.35, pi = 1,    pi_h = 0.75),
  c(m = 3,  p = 0.20, r = 0.35, pi = 0.75, pi_h = 0.75),
  c(m = 10, p = 0.45, r = 0.50, pi = 1,    pi_h = 0),
  c(m = 10, p = 0.45, r = 0.50, pi = 0.75, pi_h = 0)
))

# ------------------------------------------------------------------------------
# Asymptotische Ergebnisse über alle Parameterkombinationen und n-Werte
# ------------------------------------------------------------------------------

asymp_df <- do.call(rbind, lapply(seq_len(nrow(unique_coeffs)), function(idx) {
  params  <- unique_coeffs[idx, ]
  m_val   <- as.numeric(params["m"])
  p_val   <- as.numeric(params["p"])
  r_val   <- as.numeric(params["r"])
  pi_val  <- as.numeric(params["pi"])
  pi_h_val <- as.numeric(params["pi_h"])

  marginal <- pbinom(0:(m_val - 1), m_val, p_val)
  Sigma    <- Sigma_Star(m_val, p_val, r_val, pi_val, pi_h_val)

  do.call(rbind, lapply(UNIQUE_N, function(n_val) {
    iov_res  <- IOV_asymptotic(Sigma, marginal, n_val, m_val)
    skew_res <- Skew_asymptotic(Sigma, marginal, n_val, m_val)
    C_res    <- Cohens_asymptotic_iid(n_val, pi_val, m_val, marginal)

    data.frame(
      n          = n_val,
      m          = m_val,
      p          = p_val,
      r          = r_val,
      pi         = pi_val,
      pi_h       = pi_h_val,
      type       = "Asymptotic",
      mean_IOV   = iov_res["expectation"],
      sd_IOV     = iov_res["sd"],
      lower_IOV  = iov_res["expectation"] - iov_res["sd"],
      upper_IOV  = iov_res["expectation"] + iov_res["sd"],
      mean_Skew  = skew_res["expectation"],
      sd_Skew    = skew_res["sd"],
      lower_Skew = skew_res["expectation"] - skew_res["sd"],
      upper_Skew = skew_res["expectation"] + skew_res["sd"],
      mean_C     = C_res["expectation"],
      sd_C       = C_res["sd"],
      lower_C    = C_res["expectation"] - C_res["sd"],
      upper_C    = C_res["expectation"] + C_res["sd"],
      row.names  = NULL
    )
  }))
}))

# ------------------------------------------------------------------------------
# Sigma für Bias-Kurve (m=10, p=0.3, r=0.2, pi=0.75, MCAR)
# Wird in 04_visuals.R (Bias-Plots) benötigt
# ------------------------------------------------------------------------------

Sigma_raw <- Sigma_Star(m = 10, p = 0.3, r = 0.2, pi = 0.75, pi_h = 0)

# ------------------------------------------------------------------------------
# Speichern
# ------------------------------------------------------------------------------

save.image("Masterarbeit.RData")
