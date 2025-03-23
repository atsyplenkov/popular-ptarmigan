# libraries ---------------------------------------------------------------
library(tidyverse)    # CRAN v2.0.0
library(terra)        # CRAN v1.7-18
library(sf)           # CRAN v1.0-12
library(here)         # CRAN v1.0.1
library(atslib)       # [github::atsyplenkov/atslib]
library(patchwork)    # CRAN v1.1.2
library(mgcv)         # CRAN v1.8-42
library(report)       # CRAN v0.5.7
library(performance)  # CRAN v0.10.3

# ggplot2 settings --------------------------------------------------------
theme_set(theme_hp(base_font_family = "Verdana"))

# load data ---------------------------------------------------------------
# Sediment sampling points
seds <- 
  st_read(here("data", "vector", 
               "sediment-samples", "gizh_sediments.shp")) %>% 
  transmute(
    source = "sediments",
    il
  )

# Interpolated data from Reefmaster
reefmaster_map <- 
  rast(here("data", "raster", "lake_dem", "giz_bathymetry.tif"))

# Extract depth at sediment sampling points
seds_sf <- 
  seds %>% 
  mutate(depth = extract(reefmaster_map,
                         vect(.))[,2])

# Spline interpolation ----------------------------------------------------
# GRID
grid <- 
  reefmaster_map %>%
  as.data.frame(xy = T) %>%
  as_tibble() %>% 
  rename(depth = giz_bathymetry) %>% 
  st_as_sf(coords = c("x", "y"),
           crs = 32638)

# Interpolate
# Thin Plate Regression Spline (TPRS)
seds_df <- 
  seds_sf %>% 
  st_drop_geometry() %>% 
  cbind(st_coordinates(seds_sf)) %>% 
  as_tibble()

grid_df <- 
  grid %>% 
  st_drop_geometry() %>% 
  cbind(st_coordinates(grid)) %>% 
  as_tibble()

# Fit model
fit_gam_reml <- mgcv::gam(il ~ s(X, Y, k = 10) + depth,
                          data = seds_df,
                          method = "REML")

# Model params
r <- report(fit_gam_reml)
as.report_table(r, summary = T)
performance::performance(fit_gam_reml)

# Fit to data
grid_df$TPRS <- predict(fit_gam_reml,
                        newdata = grid_df,
                        type = "response")

# Save as SpatRast object
sed_rast <- 
  grid_df %>% 
  dplyr::select(X, Y, TPRS) %>% 
  rast(type="xyz",
       crs = st_crs(boundary)$proj4string)


# Compare modeled and measured seds ---------------------------------
seds_diff <-
  seds_sf %>% 
  mutate(tprs_il = extract(sed_rast,
                           vect(seds_sf))[,2]) %>% 
  st_drop_geometry() %>% 
  as_tibble() %>% 
  mutate(diff = tprs_il - il)

# Depth difference statistics
(sed_stat <- 
    seds_diff %>% 
    summarise(m = mean(diff, na.rm = T),
              sigma = sd(diff, na.rm = T),
              min = min(diff, na.rm = T),
              max = max(diff, na.rm = T))
)

# Scatter plot
seds_diff_scatter <-
  seds_diff %>% 
  ggplot(aes(y = tprs_il,
             x = il)) +
  geom_point(fill = "grey60",
             alpha = .4,
             # color = "white",
             shape = 21) +
  Add_1_line() +
  Add_R2(add_line = F) +
  labs(y = "Interpolated thickness, m",
       x = "Measured thickness, m",
       subtitle = "(a)")


seds_diff_density <-
  seds_diff %>% 
  ggplot(aes(x = diff,
             y = ..scaled..)) +
  geom_density(alpha = 0.35,
               n = 200,
               fill = "grey60",
               bw = 0.03) +
  scale_y_continuous(expand = expansion(add = .01)) +
  scale_x_continuous(breaks = seq(-1, 1, by = .1),
                     limits = c(-0.3, 0.3),
                     expand = expansion(add = .01)) +
  labs(x = "Thickness error, m",
       y = "Density",
       subtitle = "(b)")

# Plot together
seds_error_plot <- 
  seds_diff_scatter + seds_diff_density +
  plot_layout(widths = c(1, 1))

# save --------------------------------------------------------------------
# Sediment thickness raster
writeRaster(
  x = sed_rast,
  filename = here("R", "03_lake-bathymetry", 
                  "data", "giz_thickness.tiff"),
  gdal = c("COPY_SRC_OVERVIEWS=YES", "COMPRESS=LZW", "TILED=YES"),
  overwrite = TRUE
)

# Error plot
ggsave(
  seds_error_plot, 
  filename = here("R", "03_lake-bathymetry", "figures",
                  "lake_sediments_error_plot.tiff"),
  device = ragg::agg_png,
  dpi = 1000,
  width = 16, 
  height = 9, 
  units = "cm"
)


