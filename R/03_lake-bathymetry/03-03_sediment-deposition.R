# libraries ---------------------------------------------------------------
library(tidyverse)
library(terra)
library(here)
library(truncnorm)
library(atslib)
library(ggdist)
library(qs)

# ggplot2 settings --------------------------------------------------------
theme_set(theme_hp(base_font_family = "Verdana"))

# load data ---------------------------------------------------------------
# Lake bathymetry
bath <- 
  rast(here("data", "raster", "lake_dem", "giz_bathymetry.tif"))

# Thickness
thick <- 
  rast(here("R", "03_lake-bathymetry", "data", "giz_thickness.tiff"))

# Watershed border
ws <- 
  vect(here("data", "vector", "watershed", "giz_2022_ws-dem.shp"))

giz_area <- 1.84

# 2) Statistics -----------------------------------------------------------
# Lake volume [m³]
global(bath, sum, na.rm = T) %>% 
  round()

# Mean depth [m]
global(bath, mean, na.rm = T)
global(bath, sd, na.rm = T)

# Sed volume [m³]
(sed_vol <- 
    global(thick, sum, na.rm = T) %>% 
    round() %>%
    pull(sum))

# Mean thickness [m]
global(thick, mean, na.rm = T) %>%
  round(2)

global(thick, sd, na.rm = T) %>%
  round(2)

# 3) Simulated sediment yield ---------------------------------------------
set.seed(123)
possible_densities <- 
  rtruncnorm(10000,
             a =.75,
             b = 2,
             mean = 1.5,
             sd = 0.5)

sed_mass <- 
  mean_qi(possible_densities * sed_vol) %>% 
  mutate(across(contains("y"), ~round(.x)))

# 1989–2020 ---------------------------------------------------------------
# All possible sedimentation rates, t/yr
sed_rates <- 
  possible_densities * sed_vol / 31

# Sed rate t/yr
sed_mass %>% 
  mutate(across(contains("y"), ~(. / (31)))) %>% 
  mutate(across(contains("y"), ~round(.x)))

# Sediment yield t/km²/yr
sed_mass %>% 
  # mutate(across(contains("y"), ~(. / (20 * giz_area)))) %>%
  mutate(across(contains("y"), ~(. / (31 * giz_area)))) %>%
  mutate(across(contains("y"), ~round(.x)))

# 2001–2020 ---------------------------------------------------------------
# All possible sedimentation rates, t/yr
sed_rates_20 <- 
  possible_densities * sed_vol / 20

# Sed rate t/yr
sed_mass %>% 
  mutate(across(contains("y"), ~(. / (20)))) %>% 
  mutate(across(contains("y"), ~round(.x)))

# Sediment yield t/km²/yr
sed_mass %>% 
  mutate(across(contains("y"), ~(. / (20 * giz_area)))) %>%
  mutate(across(contains("y"), ~round(.x)))

# joint assessment --------------------------------------------------------
all_pos20 <- 
  possible_densities * sed_vol / 20

all_pos31 <- 
  possible_densities * sed_vol / 31

# Sediment yield t/km²/yr
(c(all_pos20, all_pos31)/giz_area) %>% 
  mean_qi() %>% 
  mutate(across(contains("y"), ~round(.x)))

# 4) Save -----------------------------------------------------------------
qsave(sed_rates,
      here("R", "03_lake-bathymetry", "data", "sed_rates.qs"))


