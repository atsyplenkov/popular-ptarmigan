# libraries ---------------------------------------------------------------
library(tidyverse)    # CRAN v2.0.0
library(terra)        # CRAN v1.7-18
library(sf)           # CRAN v1.0-12
library(here)         # CRAN v1.0.1
library(atslib)       # [github::atsyplenkov/atslib]
library(geomtextpath) # CRAN v0.1.1
library(patchwork)    # CRAN v1.1.2

# ggplot2 settings --------------------------------------------------------
theme_set(theme_hp(base_font_family = "Verdana"))

# load data ---------------------------------------------------------------
# Shoreline
boundary <- 
  read_sf(here("data", "vector", "shoreline", "gizh_bl_ortho.shp")) %>%
  transmute(source = "boundary",
            depth = 0) %>%
  st_transform(32638) %>%
  st_zm()

# Depth survey
measured_depths <- 
  read_sf(here("data", "vector", "depth-points", "gizh_bath-all.shp")) %>% 
  transmute(source = "measured",
            depth = as.numeric(field_3)) %>% 
  st_transform(32638) %>% 
  distinct(.keep_all = T)

# Boundary points
boundary_points <- 
  st_cast(boundary, "POINT")

depths <- 
  rbind(boundary_points, measured_depths) %>%
  cbind(., st_coordinates(.))

# Interpolated data from Reefmaster
reefmaster_map <- 
  rast(here("data", "raster", "lake_dem", "giz_bathymetry.tif"))

# compare modelled and measured depths ------------------------------------
depths_diff <- 
  measured_depths %>% 
  mutate(rm_depth = extract(reefmaster_map,
                            vect(measured_depths))[,2]) %>% 
  st_drop_geometry() %>% 
  as_tibble() %>% 
  mutate(diff = rm_depth - depth)

# Depth difference statistics
bath_hdr <- 
  hdrcde::hdr(den = density(depths_diff$diff,
                            bw = .13,
                            na.rm = T))
(depth_stat <- 
    depths_diff %>% 
    summarise(m = mean(diff, na.rm = T),
              med = median(diff, na.rm = T),
              mode = bath_hdr$mode,
              mad = mad(diff, na.rm = T),
              sigma = sd(diff, na.rm = T),
              min = min(diff, na.rm = T),
              max = max(diff, na.rm = T))
)

# plots -------------------------------------------------------------------
# Scatter plot
depth_diff_scatter <-
  depths_diff %>% 
  ggplot(aes(y = rm_depth, x = depth)) +
  geom_point(fill = "grey60",
             alpha = .4,
             # color = "white",
             shape = 21) +
  Add_1_line() +
  Add_R2(add_line = F) +
  # scale_x_continuous(breaks = seq(25, 70, 5),
  #                    expand = c(0.01, 0.5)) +
  # scale_y_continuous(breaks = seq(25, 70, 5),
  #                    expand = c(0.01, 0.5)) +
  labs(y = "Interpolated depth, m",
       x = "Measured depth, m",
       subtitle = "(a)")

depth_diff_density <-
  depths_diff %>% 
  ggplot(aes(x = diff,
             y = ..scaled..)) +
  geom_density(alpha = 0.35,
               n = 200,
               fill = "grey60",
               bw = 0.13) +
  geom_textvline(xintercept = bath_hdr$mode,
                 label = paste0("−", round(bath_hdr$mode, 2)),
                 vjust = -.1,
                 linetype = "dotdash") +
  scale_x_continuous(breaks = seq(-3, 1, by = .5),
                     limits = c(-2, 1),
                     expand = expansion(add = .01)) +
  scale_y_continuous(expand = expansion(add = .01)) +
  labs(x = "Depth error, m",
       y = "Density",
       subtitle = "(b)")

# plot together
depth_error_plot <- 
  depth_diff_scatter + depth_diff_density +
  plot_layout(widths = c(1, 1))

# save --------------------------------------------------------------------
ggsave(
  depth_error_plot, 
  filename = here("R", "03_lake-bathymetry", "figures",
                  "lake_dem_error_plot.tiff"),
  device = ragg::agg_png,
  dpi = 1000,
  width = 16, 
  height = 9, 
  units = "cm"
)
