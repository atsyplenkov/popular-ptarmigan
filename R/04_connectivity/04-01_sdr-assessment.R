library(dplyr)
library(terra)
library(here)
library(sf)
library(atslib)
library(rusleR)
library(magrittr)
library(ggdist)

# 1) load data ------------------------------------------------------------
# Watershed
ws <-
  vect(here("data", "vector", "watershed", "giz_2022_ws-dem.shp"))

# Lake shoreline (polygon)
lake <-
  vect(here("data", "vector", "shoreline_dem", "gizh_bl_dem.shp"))

# Connectivity index
ic <-
  rast(here("data", "raster", "ic_with_targets.tif")) %>%
  terra::crop(ws, mask = T) %>%
  terra::mask(lake, inverse = T)

# 2) Estimate SDRmax ------------------------------------------------------
# SDRmax is defined as a fraction of topsoil particles finer than coarse
# sand (Vigiak et al., 2012)
cfvo_q50 <-
  get_soilgrids(ws, layer = "cfvo", pred.quantile = "Q0.5")

mean_cfvo <-
  cfvo_q50 %>%
  mean() %>%
  global(mean, na.rm = T) %>%
  divide_by(10) %>%
  pull(1)

# 3) Calculate SDR --------------------------------------------------------
# See Vigiak et al., 2012
sdr_max <- round(1 - mean_cfvo / 100, 1) # 0.9
k <- 2.0
ic_0 <- 0.5

sdr_k1_max <- sdr_max / (1 + exp((ic_0 - ic) / k))

# Overall statistics
global(
  sdr_k1_max,
  fun = function(x) {
    q50 <- quantile(x, c(0.5), na.rm = T)
    m <- mean(x, na.rm = T)

    tibble(mean = m, med = q50) %>%
      round(3)
  }
)

# Mean and CI
sdr_k1_max |>
  terra::values() %>%
  ggdist::mean_qi(na.rm = T) %>%
  mutate(across(where(is.numeric), ~ round(., 3)))


# Get SDR values for the slopes adjacent to the lake ---------------------
lake_slopes <-
  stringr::str_squish(
    "POLYGON ((42.99472280063378804 43.47123297249415685, 
    42.99499318515330515 43.47116333752526174, 42.99489095926940507 
    43.4707803561009527, 42.9955676002677194 43.47039244106963451, 
    42.9961964383145272 43.47029189949891048, 42.99656502848106499 
    43.46995436138414703, 42.99696684236003819 43.47015664624734654, 
    42.99807711135842681 43.47016679990915122, 42.99845722163060202 
    43.46986664652747123, 42.99935642326463636 43.46969644331153404, 
    42.99970671654739363 43.46952592571579999, 43.0001787651794487 
    43.46953418557834681, 43.00054089406462055 43.46939176431815355, 
    43.00067000654629368 43.46934753637965798, 43.0012556512216193 
    43.46939496781121903, 43.00183017327536561 43.46900522974006265, 
    43.00219841849219904 43.4686769645458071, 43.00268790798891416 
    43.46815557393486529, 43.00313636183891219 43.46771713981128471, 
    43.00325602581532536 43.46524615317182594, 42.99832605636241567 
    43.46494608879676491, 42.99496772677294132 43.46537068726958353, 
    42.99338873819837659 43.4679555141812628, 42.99472280063378804 
    43.47123297249415685))"
  ) |> 
  wk::wkt(crs = 4326) |> 
  sf::st_as_sf() |> 
  sf::st_transform(32638)

library(exactextractr)

lake_slopes_sdr <- 
  exact_extract(sdr_k1_max, lake_slopes, "mean")

lake_slopes

# # 4) Correct the SDR values -----------------------------------------------
# # calculate upper and lower values of SDR
# # with 1% and 99% probability
# sdr_quantiles <-
#   global(
#     sdr_k1_max,
#     fun = function(x) {
#       quantile(x, probs = c(0.005, 0.995), na.rm = T)
#     }
#   )

# # Keep only sdr values between those quantiles
# sdr_corr <-
#   sdr_k1_max %>%
#   clamp(sdr_quantiles[, 1], sdr_quantiles[, 2], values = F)

# # 5) Descriptive statistics -----------------------------------------------
# # Means and CI
# sdr_corr %>%
#   terra::values() %>%
#   ggdist::mean_qi(na.rm = T) %>%
#   mutate(across(where(is.numeric), ~ round(., 3)))

# # Max and Min
# sdr_corr %>%
#   global(
#     fun = function(x) {
#       quantile(x, probs = c(0, 1), na.rm = T)
#     }
#   ) %>%
#   mutate(across(where(is.numeric), ~ round(., 2)))

# save --------------------------------------------------------------------
writeRaster(
  x = sdr_k1_max,
  filename = here(
    "R",
    "04_connectivity",
    "data",
    "sdr_corr_2025.tiff"
  ),
  gdal = c("COPY_SRC_OVERVIEWS=YES", "COMPRESS=LZW", "TILED=YES"),
  overwrite = TRUE
)

# test -------------------------------------------------------------------

rusle_maxim <- terra::rast("../analysis/from/rusle_1989.tif")
rusle_1989 <- terra::rast(
  "R/07_rusle-estimation/data/rusle_1989-2020.tif"
)

global(rusle_1989, "sum", na.rm = TRUE) / 10000
global(rusle_maxim, "sd", na.rm = TRUE) / 10000

se_rusle <- resample(sdr_k1_max, rusle_maxim, "near") * rusle_maxim
global(se_rusle, "sum", na.rm = TRUE) / 10000
