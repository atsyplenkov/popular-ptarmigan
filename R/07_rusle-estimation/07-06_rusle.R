library(tidyverse)
library(terra)
library(Rsagacmd)
library(rusleR)
library(sf)
library(here)

# initiate a saga object
saga <- Rsagacmd::saga_gis(raster_backend = "terra")

# 1) Load input data ------------------------------------------------------
# DEM
dem <-
  terra::rast(here::here(
    "data/raster/basin_dem/giz_2022_dem_r5-f-lake.tif"
  ))

# K factor
k_mean <-
  terra::rast(here::here("R/07_rusle-estimation/data/giz_k_mean.tif"))

# R
load(here::here("R/07_rusle-estimation/data/giz_r.Rdata"))

# LULC
lulc <-
  sf::st_read(
    here::here(
      "data/vector/lulc-ortho/giz_2022_lulc-ortho_25may.shp"
    ),
    quiet = TRUE
  )

# Watershed border
ws <-
  terra::vect(here::here(
    "data",
    "vector",
    "watershed",
    "giz_2022_ws-dem.shp"
  ))

# Lake shoreline (polygon)
lake <-
  terra::vect(here::here("data/vector/shoreline_dem/gizh_bl_dem.shp"))

# Clip dem to ws border
dem_ws <-
  dem |>
  terra::crop(ws, mask = T)

# SDR
sdr <-
  #FIXME:
  # REPLACE SDR path
  terra::rast(here::here(
    "R",
    "04_connectivity",
    "data",
    "sdr_rockfall_corr.tiff"
  )) |>
  terra::crop(ws, mask = T) |>
  terra::mask(lake, inverse = T)

# 2) Compute LS -----------------------------------------------------------
ls_100 <-
  rusleR::ls_alpine(dem_ws, threshold = 120)

# 3) C-factor -------------------------------------------------------------
c_val <-
  lulc |>
  dplyr::mutate(
    lulc = dplyr::case_when(
      id == 1 ~ "shrub",
      id == 2 ~ "patchy shrub",
      id == 3 ~ "forest",
      id == 4 ~ "natural grassland"
    )
  ) |>
  dplyr::mutate(
    C = dplyr::case_when(
      id == 1 ~ 0.001, # shrub
      id == 2 ~ 0.003, # patchy shrub
      id == 3 ~ 0.001, # forest
      id == 4 ~ 0.01, # grassland/pasture
      TRUE ~ NA_real_
    )
  )

# Rasterize
c_rast <- terra::rast(
  terra::vect(c_val),
  ncols = terra::ncol(dem_ws),
  nrows = terra::nrow(dem_ws)
)

c_rast <- terra::rasterize(terra::vect(c_val), c_rast, "C")

# 4) Resample -------------------------------------------------------------
# C-factor
c_res <-
  c_rast |>
  terra::resample(dem_ws, "near") |>
  terra::crop(ws, mask = T)

# K-factor
k_res <-
  k_mean |>
  terra::resample(dem_ws, "near") |>
  terra::crop(ws, mask = T)

# 5) R-factor -------------------------------------------------------------
R_1989 <-
  giz_erosivity |>
  dplyr::filter(period == "1989-2020") |>
  dplyr::pull(fit)

R_vars <-
  giz_erosivity |>
  dplyr::filter(period == "1989-2020") |>
  dplyr::select(fit:upr) |>
  as.list()

# 6) RUSLE ----------------------------------------------------------------
rusle_1989 <-
  ls_100 * c_res * k_res * R_1989

terra::names(rusle_1989) <- "RUSLE_1989-2020"

# Stats
terra::global(rusle_1989, mean, na.rm = T) |>
  round(1)

terra::global(rusle_1989, sd, na.rm = T) |>
  round(1)

terra::global(rusle_1989, median, na.rm = T)

terra::global(rusle_1989 / 10000, sum, na.rm = T) |>
  round(1)

# RUSLE considering the uncertainty in R
ckls <-
  ls_100 * c_res * k_res

R_vars |>
  purrr::map(~ (. * ckls)) |> # t/ha/yr
  purrr::map(~ (. * sdr / 10000)) |> # t/yr
  purrr::map_df(~ terra::global(., sum, na.rm = T)) |> # t/yr
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), round)) |> # t/yr
  magrittr::divide_by(1.84) |> # t/km2/yr
  # magrittr::multiply_by(2.5) |>  # total sed export
  round(1)

# 7) rusle in mm ----------------------------------------------------------
#FIXME:
# There was some discussion re: soil density. Double-check

# Get mean topsoil bulk density
dens <-
  rusleR::get_soilgrids(ws, layer = "bdod")

# Mean topsoil density in AOI
dens_mean <-
  dens |>
  terra::mean() |>
  magrittr::divide_by(100) |>
  terra::global("mean", na.rm = T) # t/m³
#> 1.216963

# Mean soil loss rate in mm/yr
R_vars |>
  purrr::map(~ (. * ckls)) |>
  purrr::map_dfr(~ terra::global(.x, mean, na.rm = T)) |>
  dplyr::as_tibble() |>
  dplyr::rename(rate = 1) |>
  dplyr::mutate(type = c("mean", "lwr", "upr"), .before = rate) |>
  dplyr::mutate(rate_mm = (rate / dens_mean$mean) / 10) # mm/yr

# 95% interval of soil erosion rates
rusle_1989 |>
  terra::values(na.rm = T) |>
  ggdist::mean_qi() |>
  dplyr::mutate(dplyr::across(
    dplyr::contains("y"),
    ~ (. / dens_mean$mean / 10)
  )) |> # mm/yr
  dplyr::mutate(dplyr::across(
    dplyr::where(is.numeric),
    ~ round(., 2)
  ))

# 7) ----------------------------------------------------------------------
# Sediment load in t / yr
sl <-
  sdr * rusle_1989 / 10000

terra::global(sl, sum, na.rm = T)

# 8) ----------------------------------------------------------------------
terra::writeRaster(
  rusle_1989,
  "R/07_rusle-estimation/data/rusle_1989-2020.tif",
  overwrite = TRUE
)

terra::writeRaster(ls_100, "R/07_rusle-estimation/data/ls_120.tif")
