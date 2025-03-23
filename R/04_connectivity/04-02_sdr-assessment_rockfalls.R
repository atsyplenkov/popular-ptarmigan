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

# Erosion processes shapefile
pr <- 
  st_read(
    here("data", "vector", "erosion-map", "processes.shp"),
    quiet = T
  ) %>% 
  st_transform(crs(ws, proj = T)) %>% 
  st_intersection(st_as_sf(ws))

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

sdr_max <- 1 - mean_cfvo/100

# 3) Spatially-variable SDRmax --------------------------------------------
sdr_vect <- 
  pr %>% 
  transmute(
    sdr_max = ifelse(
      Group_1 == 4, 
      1, sdr_max
    )
  ) %>% 
  vect()

# dummy raster object
raster_obj <- 
  rast(
    sdr_vect,
    ncols = ncol(ic),
    nrows = nrow(ic)
  )

# populate raster obj with values
sdr_max_rast <- 
  rasterize(sdr_vect, raster_obj, "sdr_max") %>% 
  resample(ic)

# 4) Calculate SDR --------------------------------------------------------
# See Vigiak et al., 2012
k <- 1
ic_0 <- 0.5

sdr_rockfall <- 
  sdr_max_rast / (1 + exp((ic_0 - ic) / k))

# Overall statistics
global(
  sdr_rockfall,
  fun = function(x){
    q50 <- quantile(x, c(0.5), na.rm = T)
    m <- mean(x, na.rm = T)
    
    tibble(mean = m,
           med = q50) %>% 
      round(3)
  }
)

# 5) Correct the SDR values -----------------------------------------------
# calculate upper and lower values of SDR
# with 1% and 99% probability
sdr_rockfall_quantiles <-
  global(
    sdr_rockfall,
    fun = function(x){
      quantile(x,
               probs = c(0.01, 0.99),
               na.rm = T)
    }
  )

# Keep only sdr values between those quantiles
sdr_rockfall_corr <-
  sdr_rockfall %>% 
  clamp(sdr_rockfall_quantiles[,1],
        sdr_rockfall_quantiles[,2],
        values = F) 

# 5) Descriptive statistics -----------------------------------------------
# Means and CI
sdr_rockfall_corr %>% 
  terra::values() %>% 
  ggdist::mean_qi(na.rm = T) %>% 
  mutate(across(where(is.numeric),
                ~round(., 3)))

# Max and Min
sdr_rockfall_corr %>% 
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
  x = sdr_rockfall_corr,
  filename =  here("R", "04_connectivity", "data", "sdr_rockfall_corr.tiff"),
  gdal = c("COPY_SRC_OVERVIEWS=YES", "COMPRESS=LZW", "TILED=YES"),
  overwrite = TRUE
)
