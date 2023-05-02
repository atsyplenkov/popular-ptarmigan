# libraries ---------------------------------------------------------------
library(tidyverse)     # CRAN v2.0.0
library(exactextractr) # [github::isciences/exactextractr] v0.10.0
library(rstatix)       # CRAN v0.7.2
library(terra)         # CRAN v1.7-18
library(whitebox)      # CRAN v2.3.0
library(here)          # CRAN v1.0.1
library(sf)            # CRAN v1.0-12
library(atslib)        # [github::atsyplenkov/atslib]

theme_set(theme_hp(base_font_family = "Verdana"))

# whitebox setup ----------------------------------------------------------
# change it according to your system setup
wbt_init(exe_path = "C:/Soft/WBT/whitebox_tools.exe")

# load files --------------------------------------------------------------
# Crop area, i.e. AOI
mask <- vect(here("data", "vector","dod_area", "dod_area.shp"))

# AOI are in m2
terra::expanse(mask)

# Path to DEMs
dem_2021_path <-
  here("data", "raster", "dod_dems", "GE_2021_dem_cut.tiff")

dem_2022_path <-
  here("data", "raster", "dod_dems", "GE_2022_dem_cut_ref_zcor.tiff")

# UAV dates
dod_2021 <-
  as.Date("2021-05-23")

dod_2022 <-
  as.Date("2022-07-26")

# Difference in years
dod_difftime <- 
  difftime(dod_2022, dod_2021, "days") %>%
  as.numeric() %>% 
  magrittr::divide_by(365)

# Co-registration assessment
# This assessment was done manually in Global Mapper by
# Sergey Kharchenko
err <-
  c(
    -0.0045166016,
    -0.0255126953,
    -0.0144042969,
    0.0335693359,
    -0.0025634766,
    -0.0057373047,
    0.0004882812,
    -0.0059814453,
    0.0069580078,
    -0.0010986328,
    -0.0079345703,
    -0.0212402344,
    -0.0288085938,
    -0.0100097656,
    0.0102539062,
    0.0046386719,
    0.0046386719,
    -0.0069580078,
    0.0372314453,
    0.0107421875,
    0.0054931641,
    -0.0151367188,
    0.0028076172
  )

dod_sd <- 
  sd(err)

dod_mean_err <- 
  mean(err)

# 1) DoD creation ---------------------------------------------------------
# !NB DEMs are already co-registered using Nuuth and Kaab method

# Load 2021 year cropped to AOI
dem2021 <- 
  rast(dem_2021_path) %>% 
  crop(mask, mask = T)

# Load 2022 year cropped to AOI
dem2022 <- 
  rast(dem_2022_path) %>% 
  crop(mask, mask = T)

# Dem of Difference
dod_rast <- 
  dem2022 - dem2021

# Apply level of detection
dod_rast_2sd <- 
  ifel(
    dod_rast > -2*dod_sd & dod_rast < 2*dod_sd,
    0,
    dod_rast
  )

# save DoDs
dod_path <- 
  here("R", "02_dem-of-difference", "data",
       "dod_2021-2022_raw.tiff")

dod_det_path <- 
  here("R", "02_dem-of-difference", "data",
       "dod_2021-2022_2sd.tiff")

writeRaster(
  x = dod_rast_2sd,
  filename = dod_det_path,
  gdal = c("COPY_SRC_OVERVIEWS=YES", "COMPRESS=LZW", "TILED=YES"),
  overwrite = TRUE
)

# 2) Buffer zone calculation ----------------------------------------------
# Fill sinks using Breach Depressions algorithm
breach_dem_path <- 
  here("R", "02_dem-of-difference", "data", "GE_2022_dem-breach.tiff")

wbt_breach_depressions(
  dem = dem_2022_path,
  output = breach_dem_path,
  compress_rasters = TRUE
)

# The following code chunks will find all streams with catchment 
# area more than 2m². Then calculate Strahler Order, vectorize 
# them and keep only significant streams.

# D8 Flow pointer aka Flow Directions
d8_pointer_path <- 
  here("R", "02_dem-of-difference", "data", "GE_2022_dem-d8pointer.tiff")

wbt_d8_pointer(
  dem = breach_dem_path,
  output = d8_pointer_path,
  compress_rasters = TRUE
)

# Flow Accumulation
d8_fa_path <- 
  here("R", "02_dem-of-difference", "data", "GE_2022_dem-d8fa.tiff")
  
wbt_d8_flow_accumulation(
  input = d8_pointer_path,
  output = d8_fa_path,
  pntr = TRUE,
  out_type = 'catchment area',
  compress_rasters = TRUE
)

# Extract streams
streams_path <- 
  here("R", "02_dem-of-difference", "data", "GE_2022_dem-d8streams2.tiff")

wbt_extract_streams(
  flow_accum = d8_fa_path,
  output = streams_path,
  threshold = 2, # 2m² threshold
  compress_rasters = TRUE
)

# Strahler Order
str_ord_path <- 
  here("R", "02_dem-of-difference", "data", "GE_2022_dem-d8strahler.tiff")

wbt_strahler_stream_order(
  streams = streams_path,
  d8_pntr = d8_pointer_path,
  output = str_ord_path,
  compress_rasters = TRUE
)

# Vectorize streams
vct_streams_path <- 
  here("R", "02_dem-of-difference", "data", "vct_streams", "GE_2022_streams-2m.shp")

wbt_raster_streams_to_vector(
  streams = str_ord_path,
  d8_pntr = d8_pointer_path,
  output = vct_streams_path
)

# 3) Mean erosion rates ---------------------------------------------------
# Load streams vectors
streams <- 
  st_read(vct_streams_path) %>% 
  sf::st_set_crs(32638) 

# Estimate buffer as mean stream width (1m)
streams_buffer <- 
  streams %>% 
  mutate(sss = STRM_VAL) %>% 
  dplyr::filter(sss > 2) %>% 
  sf::st_union() %>% 
  st_buffer(0.5)

buffer_smpl <- 
  rmapshaper::ms_simplify(streams_buffer, keep = 0.3)

# crop to AOI
buffer_aoi <- 
  buffer_smpl %>% 
  vect() %>% 
  crop(mask)

# Save buffer
buffer_aoi %>%
  st_as_sf() %>% 
  st_write(
    here("R", "02_dem-of-difference", "data", "vct_buffer",
         "vct_buffer.shp"),
    delete_dsn = TRUE
  )

# Gully erosion
# Extract all DoD values within a buffer area
dod_vals <- 
  exact_extract(
    raster::raster(dod_rast), 
    st_as_sf(buffer_aoi)
  )

# Summary statistics in m/yr of 
dod_vals %>%
  pluck(1) %>% 
  as_tibble() %>% 
  filter(value < -2*dod_sd | value > 2*dod_sd) %>%
  mutate(rate = value / dod_difftime) %>% 
  rstatix::get_summary_stats(
    rate,
    show = c("mean", "sd", "median", "mad")
  )

# 4) Plots ----------------------------------------------------------------
dod_vals_2sd <- 
  dod_vals %>%
  pluck(1) %>% 
  as_tibble() %>% 
  mutate(bin = ntile(value, 102),
         bin = as_factor(bin)) %>% 
  filter(value < -2*dod_sd | value > 2*dod_sd)
  
dod_streams_plot <-
  dod_vals %>%
  pluck(1) %>% 
  ggplot() +
  geom_histogram(aes(value,
                     fill = ..x..),
                 binwidth = 0.025) +
  geom_histogram(
    data = dod_vals_2sd %>% 
      filter(!between(value, -0.0375, 0.0375)),
    aes(value),
    fill = NA,
    binwidth = 0.025, 
    linewidth = 0.3,
    color = "black"
  ) +
  scale_y_continuous(
    labels = label_number(scale = 1/10000),
    expand =  expansion(mult = 0.01)
  ) +
  scale_x_continuous(
    expand =  expansion(add = .01), 
    limits = c(-1, 1) 
  ) +
  cols4all::scale_fill_continuous_c4a_div(
    "tableau.red_blue_diverging",
    limits = c(-1, 1),
    breaks = seq(-1, 1, by = 0.25),
    labels = seq(-1, 1, by = 0.25),
    mid = 0
  ) +
  labs(
    x = "DoD value [m]",
    y = expression("Frequency [×10"^4*"]"),
    fill = "DoD [m]"
  ) +
  theme(
    legend.direction = "vertical",
    legend.position = c(0.9, 0.7)
  )

# save --------------------------------------------------------------------
ggsave(
  dod_streams_plot, 
  filename = here("R", "02_dem-of-difference", "figures",
                  "dod_hist_plot.tiff"),
  device = ragg::agg_png,
  dpi = 1000,
  width = 8, 
  height = 7, 
  units = "cm"
)
