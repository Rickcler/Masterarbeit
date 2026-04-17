# ==============================================================================
# 00_setup.R
# Pakete und globale Einstellungen
# ==============================================================================

# Pakete (einmalig installieren falls nötig)
# install.packages(c("dplyr", "tidyr", "ggplot2", "patchwork"))

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# Globale Simulationsparameter
SEED    <- 123
N_REPS  <- 1000
UNIQUE_N <- c(50, 100, 250, 500, 1000)

