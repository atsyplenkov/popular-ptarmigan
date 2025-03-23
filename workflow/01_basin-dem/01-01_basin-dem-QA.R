# libraries ---------------------------------------------------------------
library(tidyverse)    # CRAN v2.0.0
library(terra)        # CRAN v1.7-18
library(hdrcde)       # CRAN v3.4
library(here)         # CRAN v1.0.1
library(sf)           # CRAN v1.0-12
library(atslib)       # [github::atsyplenkov/atslib]
library(geomtextpath) # CRAN v0.1.1
library(patchwork)    # CRAN v1.1.2

# ggplot2 settings --------------------------------------------------------
theme_set(theme_hp(base_font_family = "Verdana"))

# 1) load data ------------------------------------------------------------
# Basin DEM
dem <- 
  rast(here("data", "raster", "basin_dem", "giz_2022_dem_r5-f-lake.tif"))

# DGPS points to validate DEM quality
dgps <- 
  st_read(here("data", "vector", "dgps_points", "dgps_points.shp")) %>% 
  drop_na() %>% 
  st_transform(crs(dem, proj = T)) %>% 
  mutate(DEM_1 = terra::extract(dem, vect(.))[,2])

# 2) Altitude difference --------------------------------------------------
dgps_diff <- 
  dgps %>% 
  st_drop_geometry() %>% 
  as_tibble() %>% 
  mutate(basis_H = H[which.min(H)],
         basis_dem = DEM_1[which.min(H)]) %>% 
  mutate(H_ = H - basis_H,
         dem = DEM_1 - basis_dem) %>% 
  filter(H_ > 25) %>%
  # filter(between(H_, 25, 70)) %>%
  mutate(diff = dem - H_)

# Altitude difference statistics
dgps_hdr <- 
  hdr(den = density(dgps_diff$diff, bw = .2))

dgps_stat <- 
  dgps_diff %>% 
  summarise(med = median(diff),
            mean = mean(diff),
            sigma = sd(diff),
            mad = mad(diff),
            mode = dgps_hdr$mode)

# 3) Plots ----------------------------------------------------------------
# Scatter plot
dgps_diff_scatter <- 
  dgps_diff %>% 
  ggplot(aes(y = H_, x = dem)) +
  Add_1_line() +
  geom_point(fill = "grey60",
             # color = "white",
             shape = 21) +
  Add_R2(add_line = F) +
  scale_x_continuous(breaks = seq(25, 70, 5),
                     expand = c(0.01, 0.5)) +
  scale_y_continuous(breaks = seq(25, 70, 5),
                     expand = c(0.01, 0.5)) +
  labs(x = "DTM, m above local zero",
       y = "DGPS, m above local zero",
       subtitle = "(a)")

# Density plot
dgps_diff_density <-
  dgps_diff %>% 
  ggplot(aes(x = diff,
             y = ..scaled..)) +
  geom_density(alpha = 0.35,
               n = 200,
               fill = "grey60",
               bw = 0.2) +
  geom_textvline(xintercept = dgps_stat$mode,
                 label = paste0("−", round(dgps_stat$mode, 2)),
                 vjust = -.1,
                 linetype = "dotdash") +
  scale_x_continuous(breaks = seq(-2, 2, by = .4),
                     limits = c(-1.5, 1.5),
                     expand = expansion(add = .01)) +
  scale_y_continuous(expand = expansion(add = .01)) +
  labs(x = expression(atop("Elevation error, m", 
                       italic("(H"["DTM"]-"H"["DGPS"]*")"))),
       y = "Density",
       subtitle = "(b)")

# Plot together
dem_error_plot <- 
  dgps_diff_scatter + dgps_diff_density +
  plot_layout(widths = c(1, 1)) 

dem_error_plot

# save --------------------------------------------------------------------
ggsave(
  dem_error_plot, 
  filename = here("R", "01_basin-dem", "figures",
                  "basin_dem_error_plot.tiff"),
  device = ragg::agg_png,
  dpi = 1000,
  width = 16, 
  height = 9, 
  units = "cm"
)
