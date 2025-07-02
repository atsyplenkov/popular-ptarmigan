# Initial mean annual scenario for Gizhgit was estimated by
# M. Uspenskii file /analysis/from/rusle_1989.tif
# Here, just some figures for the article are estimated

library(dplyr)
library(terra)

# 1 July 2025 ------------------------------------------------------------
sdr <- terra::rast("R/04_connectivity/data/sdr_corr_2025.tiff")
rusle_1989 <-
  terra::rast("../analysis/from/rusle_1989.tif") |>
  terra::resample(sdr, "near")

# Mean annual situation --------------------------------------------------
# Soil loss stats
rusle_1989 |>
  terra::global(\(x) {
    dplyr::tibble(
      mean = terra::mean(x, na.rm = TRUE),
      sd = stats::sd(x, na.rm = TRUE),
      sum = sum(x, na.rm = TRUE) / 10000
    )
  })

# Sediment export by RUSLE in t/yr
rusle_se <- sdr * rusle_1989 / 10000
terra::global(rusle_se, sum, na.rm = T)

# in t/km2/yr
terra::global(rusle_se, sum, na.rm = T) / 1.84

# Various sediment export, considering rainfall variability --------------
MEAN_R <- 807
LWR_R <- 689
UP_R <- 924

ckls <- rusle_1989 / MEAN_R

rusle_1989_lwr <- ckls * LWR_R
names(rusle_1989_lwr) <- "1989 lower"
rusle_1989_up <- ckls * UP_R
names(rusle_1989_up) <- "1989 upper"

rusle_se_ci <-
  #fmt: skip
  c(rusle_1989_lwr, rusle_1989_up) *
  sdr / 10000

rusle_se_ci |>
  terra::global(\(x) {
    dplyr::tibble(
      sum = sum(x, na.rm = TRUE),
      sy = sum / 1.84
    )
  })
