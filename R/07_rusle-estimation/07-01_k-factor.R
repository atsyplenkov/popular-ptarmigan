library(tidyverse)
library(terra)
library(sf)
library(here)

# 1) Load AOI boundary -------------------------------------------------
ws <- sf::st_read(
  "data/vector/watershed/giz_2022_ws-dem.shp",
  quiet = TRUE
)

# 2) Sand -----------------------------------------------------------------
# Mean
l_sand_mean <- list.files(
  here::here("data", "raster", "soilgrid", "mean"),
  pattern = "sand",
  full.names = T
)

#fmt: skip
sand_mean <-
  l_sand_mean |>
  terra::rast() |>
  terra::mean() / 10

# 0.05
l_sand_lwr <- list.files(
  here::here("data", "raster", "soilgrid", "lwr"),
  pattern = "sand",
  full.names = T
)

#fmt:skip
sand_lwr <-
  l_sand_lwr |> 
  terra::rast() |>
  terra::mean() / 10

# 0.95
l_sand_upr <- list.files(
  here::here("data", "raster", "soilgrid", "upr"),
  pattern = "sand",
  full.names = T
)

#fmt:skip
sand_upr <-
  l_sand_upr |> 
  terra::rast() |> 
  terra::mean() / 10


# 3) Clay -----------------------------------------------------------------
# Mean
l_clay_mean <- list.files(
  here::here("data", "raster", "soilgrid", "mean"),
  pattern = "clay",
  full.names = T
)

#fmt:skip
clay_mean <-
  l_clay_mean |>
  terra::rast() |>
  terra::mean() / 10

# 0.05
l_clay_lwr <- list.files(
  here::here("data", "raster", "soilgrid", "lwr"),
  pattern = "clay",
  full.names = T
)

#fmt:skip
clay_lwr <-
  l_clay_lwr |>
  terra::rast() |>
  terra::mean() / 10

# 0.95
l_clay_upr <- list.files(
  here::here("data", "raster", "soilgrid", "upr"),
  pattern = "clay",
  full.names = T
)

#fmt:skip
clay_upr <-
  l_clay_upr |>
  terra::rast() |>
  terra::mean() / 10

# 4) Silt -----------------------------------------------------------------
# Mean
l_silt_mean <- list.files(
  here::here("data", "raster", "soilgrid", "mean"),
  pattern = "silt",
  full.names = T
)

#fmt:skip
silt_mean <-
  l_silt_mean |>
  terra::rast() |>
  terra::mean() / 10

# 0.05
l_silt_lwr <- list.files(
  here::here("data", "raster", "soilgrid", "lwr"),
  pattern = "silt",
  full.names = T
)

#fmt:skip
silt_lwr <-
  l_silt_lwr |>
  terra::rast() |>
  terra::mean() / 10

# 0.95
l_silt_upr <- list.files(
  here::here("data", "raster", "soilgrid", "upr"),
  pattern = "silt",
  full.names = T
)

#fmt:skip
silt_upr <-
  l_silt_upr |>
  terra::rast() |>
  terra::mean() / 10

# Organic Matter
# Mean
l_soc_mean <- list.files(
  here::here("data", "raster", "soilgrid", "mean"),
  pattern = "soc",
  full.names = T
)

#fmt:skip
soc_mean <-
  l_soc_mean |>
  terra::rast() |>
  terra::mean() / 100

# 0.05
l_soc_lwr <- list.files(
  here::here("data", "raster", "soilgrid", "lwr"),
  pattern = "soc",
  full.names = T
)

#fmt:skip
soc_lwr <-
  l_soc_lwr |>
  terra::rast() |>
  terra::mean() / 100

# 0.95
l_soc_upr <- list.files(
  here::here("data", "raster", "soilgrid", "upr"),
  pattern = "soc",
  full.names = T
)

#fmt:skip
soc_upr <-
  l_soc_upr |>
  terra::rast() |>
  terra::mean() / 100

# 5) K factor equation ----------------------------------------------------
# K factor equation from Williams and Renard (1983) as cited in Chen et al. (2011)
# in ton hours per megajoules per millimetre

k_factor <- function(sand, clay, silt, soc) {
  # x[1] - sand, %
  # x[2] - clay, %
  # x[3] - silt, %
  # x[4] - soc, %

  Sa <- sand
  Si <- silt
  Cl <- clay
  sn <- 1 - (sand / 100)
  C <- soc

  fsand <- 0.3 * exp(0.0256 * Sa * (1 - Si / 100))
  fclsi <- (Si / (Cl + Si))^0.3
  forg <- 1 - ((0.25 * C) / (C + exp(3.72 - 2.95 * C)))
  fhisand <- 1 - ((0.7 * sn) / (sn + exp(-5.51 + 22.9 * sn)))

  kfactor <- 0.1317 * (0.2 + fsand * fclsi * forg * fhisand)

  return(kfactor)
}

# Mean
k_mean <- k_factor(sand_mean, clay_mean, silt_mean, soc_mean)

# 0.05
k_lwr <- k_factor(sand_lwr, clay_lwr, silt_lwr, soc_lwr)

# 0.95
k_upr <- k_factor(sand_upr, clay_upr, silt_upr, soc_upr)

# 6) Uncertainty ----------------------------------------------------------
sand_iqr <- sand_upr - sand_lwr
sand_sd <- sand_iqr / 1.35
sand_iqr / sand_mean


terra::global(k_upr, mean, na.rm = T)

# Save --------------------------------------------------------------------
terra::writeRaster(
  k_mean,
  filename = "R/07_rusle-estimation/data/giz_k_mean.tif",
  overwrite = T
)

terra::writeRaster(
  k_lwr,
  filename = "R/07_rusle-estimation/data/giz_k_lwr.tif",
  overwrite = T
)

terra::writeRaster(
  k_upr,
  filename = "R/07_rusle-estimation/data/giz_k_upr.tif",
  overwrite = T
)
