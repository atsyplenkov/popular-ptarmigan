# libraries ---------------------------------------------------------------
library(tidyverse)
library(sf)
library(terra)
library(here)
library(qs)
library(ggdist)
library(truncnorm)

# 1) load data ------------------------------------------------------------
# Watershed
ws <- 
  vect(here("data", "vector", "watershed", "giz_2022_ws-dem.shp"))

# Lake shoreline (polygon)
lake <- 
  vect(here("data", "vector", "shoreline_dem", "gizh_bl_dem.shp"))

# Basin DEM
dem <- 
  rast(here("data", "raster", "basin_dem", "giz_2022_dem_r5-f-lake.tif"))

# Erosion processes shapefile
pr <- 
  st_read(
    here("data", "vector", "erosion-map", "erosion-map.shp"),
    quiet = T
  ) %>% 
  st_transform(crs(ws, proj = T)) %>% 
  st_intersection(st_as_sf(ws))

# SDR
sdr_corr <- 
  rast(here("R", "04_connectivity", "data", "sdr_rockfall_corr.tiff")) %>% 
  resample(dem, "cubic")

# All possible sedimentation rates, t/yr
sed_rates <- qs::qread(here("R", "03_lake-bathymetry", "data", "sed_rates.qs"))

giz_area <- 1.84

# 2) Mean erosion rates ---------------------------------------------------
# in m / yr
rates <- 
  tibble(
    `gullies` = 0.076, #!!!!
    `soil creep` = 0.0015,
    none = 0,
    `sheet wash` = 0.0015 + 0.0015, # since creep also applies here
    `rockfalls` = 0.013 #!!!!
  ) %>% 
  gather(type, rate)

# Convert to raster
pr_vect <- 
  pr %>% 
  left_join(rates,
            by = "type") %>% 
  vect()

pr_rast <- 
  rasterize(
    pr_vect,
    dem,
    "rate"
  )

# 3) Calculate total loss -------------------------------------------------
# Slope in radians
dem_slope <- 
  terrain(dem,
          unit = "radians")

# Slope cos
dem_cos <- cos(dem_slope)

# Tot loss in m³/yr
loss <- pr_rast / dem_cos

global(loss,
       sum,
       na.rm = T)

# 4) Volumetric sediment export -------------------------------------------
# Actual total loss with SDR
# loss_sdr <-  sdr_k1 * loss
loss_sdr <-  sdr_corr * loss

# Volumetric sediment export m³/yr
global(loss_sdr,
       sum,
       na.rm = T)

# 5) Bulk density rasters -------------------------------------------------
bd_types <- 
  pr %>% 
  transmute(
    bd_type = case_when(
      type == "rockfalls" ~ 1L,
      TRUE ~ 2L
    )
  ) %>% 
  vect() %>% 
  rasterize(dem, "bd_type")

# 6) Bulk density generations ---------------------------------------------
N_mc <- 500

set.seed(123)
bd_argillites <- 
  rtruncnorm(
    N_mc,
    a = 2.2,
    b = 4,
    mean = 2.55,
    sd = 1.396
  )

bd_soils <- 
  rtruncnorm(
    N_mc,
    a = 0.732,
    b = 1.679,
    mean = 1.217,
    sd = 1.603
  )

# 7) Spatial generations of bulk densities --------------------------------
gen_bd_geol <-
  function(i){
    
    # Reclass
    m <- c(0, 1, bd_argillites[i],
           1, 2, bd_soils[i])
    
    rclmat <- matrix(m, ncol=3, byrow=TRUE)
    
    geol_cl <- 
      bd_types %>% 
      classify(rclmat, include.lowest = T) 
    
    names(geol_cl) <- paste0("db_", i)
    
    if (i %% 10 == 0) {
      cat(i, "...")
    }
    
    return(geol_cl)
  }

db_rast <- 
  seq_len(N_mc) %>% 
  map(~gen_bd_geol(i = .x))

# 8) All possible sed export in t/yr --------------------------------------
library(tictoc)

tic()
loss_sdr_mc <- 
  db_rast %>% 
  map(~(.x * loss_sdr)) %>%
  map_df(~global(.x, sum, na.rm = T))
toc()
beepr::beep()


# QI of sed export in t/yr
mean_qi(loss_sdr_mc)

# QI of sed yield in t/km2/yr
mean_qi(loss_sdr_mc / giz_area)

# All possible loss in t/yr (without sdr)
loss_mc <- 
  db_rast[1:25] %>% 
  map(~(.x * loss)) %>%
  map_df(~global(.x, sum, na.rm = T))

# QI of loss in t/yr
mean_qi(loss_mc)

# 9) Save -----------------------------------------------------------------
qsave(
  loss_sdr_mc,
  here("R", "05_erosion-mapping","data", "loss_sdr_t-yr.qs")
)

qsave(
  loss_mc,
  here("R", "05_erosion-mapping","data", "loss_t-yr.qs")
)

