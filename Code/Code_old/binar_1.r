# ----------------------------------------------------------------------
# BinAR(1)-Prozess: Theoretisches ordinales Kappa nach Cohen
# ----------------------------------------------------------------------

# Parameter
m <- 10         # Anzahl der Versuche (Zustandsraum {0,...,m})
p <- 0.45        # Erfolgswahrscheinlichkeit der marginalen Binomialverteilung
r_vals <- c(-0.5, -0.2, 0.2, 0.4, 0.6, 0.8, 0.95)   # verschiedene Autokorrelationsparameter
max_h <- 100    # maximale Verzögerung

# ----------------------------------------------------------------------
# Hilfsfunktionen (aus der Aufgabenstellung)
# ----------------------------------------------------------------------

# Gemeinsame bedingte Wahrscheinlichkeit P(Y_{t+h}=j | Y_t=i)
lag_h_conditional <- function(m, p, r, h = 1) {
  beta  <- p * (1 - r^h)
  alpha <- beta + r^h
  res <- matrix(0, nrow = m + 1, ncol = m + 1)
  for (i in 0:m) {
    for (j in 0:m) {
      smaller <- min(i, j)
      bigger  <- max(0, i + j - m)
      res[i + 1, j + 1] <- sum(sapply(bigger:smaller, function(n) {
        choose(i, n) * choose(m - i, j - n) *
          alpha^n * (1 - alpha)^(i - n) *
          beta^(j - n) * (1 - beta)^(m - i - j + n)
      }))
    }
  }
  return(res)
}

# Gemeinsame Wahrscheinlichkeitsfunktion P(Y_t=i, Y_{t+h}=j)
lag_h_joint_pmf <- function(m, p, r, h = 1) {
  lag_h_conditional(m, p, r, h = h) * dbinom(0:m, m, p)
}


#-- Get CDF from PMF matrix
pmf_to_cdf <- function(pmf) {
  cdf <- apply(pmf, 2, cumsum)
  cdf <- t(apply(cdf, 1, cumsum))
  return(cdf)
}


lag_h_joint_cdf <- function(m, p, r, h = 1) {
  return(pmf_to_cdf(lag_h_joint_pmf(m,p, r, h=h)))
}

# ----------------------------------------------------------------------
# Berechnung von Cohens ordinalem Kappa für gegebenes h
# κ_ord(h) = (∑_{j=0}^{m-1} (p_{jj}(h) - p_j^2)) / (∑_{i=0}^{m-1} p_i (1-p_i))
# ----------------------------------------------------------------------
kappa_ord <- function(m, p, r, h) {
  # gemeinsame CDF für Lag h
  joint_cdf <- lag_h_joint_cdf(m, p, r, h)
  # Randwahrscheinlichkeiten (stationär, daher gleich)
  marg_cdf <- pbinom(0:m, m, p)
  # Summe über j = 0 ... m-1
  idx <- 1:m   # entspricht Werten 0,1,...,m-1
  f_jj <- sapply(idx, function(j) joint_cdf[j, j])   # Diagonalelemente
  f_j  <- marg_cdf[idx]                              # Randwahrscheinlichkeiten
  # Zähler und Nenner
  numerator   <- sum(f_jj - f_j^2)
  denominator <- sum(f_j * (1 - f_j))
  return(numerator / denominator)
}

# ----------------------------------------------------------------------
# Berechnung der κ-Werte für alle h und alle r
# ----------------------------------------------------------------------
kappa_matrix <- matrix(NA, nrow = max_h, ncol = length(r_vals))
colnames(kappa_matrix) <- paste0("r = ", r_vals)

for (i in seq_along(r_vals)) {
  r <- r_vals[i]
  for (h in 1:max_h) {
    kappa_matrix[h, i] <- kappa_ord(m, p, r, h)
  }
}

# ----------------------------------------------------------------------
# Grafik: Punkte mit verbindenden Linien, breites Format
# ----------------------------------------------------------------------
# Farben für die verschiedenen r
library(ggplot2)
library(tidyr)   # für pivot_longer

# 1. Daten vorbereiten: aus Matrix einen langen data.frame machen
#    Annahme: h_werte ist Vektor mit den Zeilenwerten
h_werte <- 1:max_h 
df <- data.frame(h = h_werte, kappa_matrix, check.names = FALSE)   # Spalten: h, r1, r2, ...


df_long <- pivot_longer(df, 
                        cols = -h,          # alle Spalten außer 'h' werden umgeformt
                        names_to = "r", 
                        values_to = "kappa")

# 2. Plot
print(ggplot(df_long, aes(x = h, y = kappa, color = r)) +
  geom_line() +          # Linien
  geom_point() +         # Punkte (optional)
  labs(x = "h", 
       y = "\u03BA(h)",
       color = "r-Values",
       title = "\u03BA(h) in relation to h",
       subtitle = "BinAR(1) with p = 0.45, m = 10 over multiple r-Values") +
  theme_minimal() +
  coord_fixed(ratio = 40))

