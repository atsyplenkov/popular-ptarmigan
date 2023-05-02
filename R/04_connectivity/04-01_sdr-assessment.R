# libraries ---------------------------------------------------------------
library(tidyverse)    # CRAN v2.0.0
library(terra)        # CRAN v1.7-18
library(here)         # CRAN v1.0.1
library(sf)           # CRAN v1.0-12
library(atslib)       # [github::atsyplenkov/atslib]
library(rusleR)       # [github::atsyplenkov/rusleR]
library(magrittr)     # CRAN v2.0.3
library(ggdist)       # CRAN v3.2.1

# 1) load data ------------------------------------------------------------
# Watershed
ws <- 
  vect(here("data", "vector", "watershed", "giz_2022_ws-dem.shp"))

# Lake shoreline (polygon)
lake <- 
  vect(here("data", "vector", "shoreline_dem", "gizh_bl_dem.shp"))

# Connectivity index
ic <-
  rast(here("data", "vector", "connectivity", "ic.tiff")) %>%
  terra::crop(ws,
              mask = T) %>% 
  terra::mask(lake,
              inverse = T)

# 2) Estimate SDRmax ------------------------------------------------------
# SDRmax is defined as a fraction of topsoil particles finer than coarse
# sand (Vigiak et al., 2012)
cfvo_q50 <- 
  get_soilgrids(ws,
                layer = "cfvo",
                pred.quantile = "Q0.5")

mean_cfvo <- 
  cfvo_q50 %>% 
  mean() %>% 
  global(mean, na.rm = T) %>% 
  divide_by(10) %>% 
  pull(1)

# 3) Calculate SDR --------------------------------------------------------
# See Vigiak et al., 2012
sdr_max <- 1 - mean_cfvo/100
k <- 1
ic_0 <- 0.5

sdr_k1_max <- sdr_max / (1 + exp((ic_0 - ic) / k))

# Overall statistics
global(
  sdr_k1_max,
  fun = function(x){
    q50 <- quantile(x, c(0.5), na.rm = T)
    m <- mean(x, na.rm = T)
    
    tibble(mean = m,
           med = q50) %>% 
      round(3)
  }
)

# 4) Correct the SDR values -----------------------------------------------
# calculate upper and lower values of SDR
# with 1% and 99% probability
sdr_quantiles <-
  global(
    sdr_k1_max,
    fun = function(x){
      quantile(x,
               probs = c(0.01, 0.99),
               na.rm = T)
    }
  )

# Keep only sdr values between those quantiles
sdr_corr <-
  sdr_k1_max %>% 
  clamp(sdr_quantiles[,1],
        sdr_quantiles[,2],
        values = F) 

# 5) Descriptive statistics -----------------------------------------------
# Means and CI
sdr_corr %>% 
  terra::values() %>% 
  ggdist::mean_qi(na.rm = T) %>% 
  mutate(across(where(is.numeric),
                ~round(., 3)))

# Max and Min
sdr_corr %>% 
  global(
    fun = function(x){
      quantile(x,
               probs = c(0, 1),
               na.rm = T)
    }
  ) %>% 
  mutate(across(where(is.numeric),
                ~round(., 2)))

# save --------------------------------------------------------------------
writeRaster(
  x = sdr_corr,
  filename =  here("R", "04_connectivity", "data", "sdr_corr.tiff"),
  gdal = c("COPY_SRC_OVERVIEWS=YES", "COMPRESS=LZW", "TILED=YES"),
  overwrite = TRUE
)
