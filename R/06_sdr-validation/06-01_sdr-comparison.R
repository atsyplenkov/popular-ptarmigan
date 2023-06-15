library(tidyverse)
library(here)
library(qs)
library(atslib)
library(terra)
library(ggdist)

theme_set(theme_hp())

# 1) Load data from previous sections -------------------------------------
# All possible sedimentation rates, t/yr
sed_rates <- qread("R/03_lake-bathymetry/data/sed_rates.qs")
loss_mc <- qread(here("R", "05_erosion-mapping","data", "loss_t-yr.qs"))

# Watershed
ws <- 
  vect(here("data", "vector", "watershed", "giz_2022_ws-dem.shp"))

# Lake shoreline (polygon)
lake <- 
  vect(here("data", "vector", "shoreline_dem", "gizh_bl_dem.shp"))

# Basin DEM
dem <- 
  rast(here("data", "raster", "basin_dem", "giz_2022_dem_r5-f-lake.tif"))

# SDR
sdr_corr <- 
  rast(here("R", "04_connectivity", "data", "sdr_rockfall_corr.tiff")) %>% 
  resample(dem, "cubic") %>% 
  crop(ws, mask = T) %>% 
  terra::mask(lake,
              inverse = T)

# 2) Estimate SDR based on processes map ----------------------------------
# Select 500 random samples
set.seed(123)
sed_r_sample <- sample(seq_len(1000), 500)

# SDR computed as total loss divided by sedimentation
loss_bath_sdr <- 
  loss_mc %>% 
  rownames_to_column(var = "id") %>% 
  as_tibble() %>% 
  rename(tot_loss = sum) %>% 
  # TEMPORARY!!!
  mutate(sed_rate = sed_rates[sed_r_sample]) %>%
  # mutate(sed_rate = sed_rates) %>% 
  expand(tot_loss, sed_rate) %>% 
  mutate(sdr = sed_rate / tot_loss)

# SDR stats
my_sdr_means <- 
  loss_bath_sdr %>% 
  mean_qi(sdr,
          na.rm = T) %>% 
  mutate(across(where(is.numeric), ~round(.x, 3)))

# 3) Compare SDRs ---------------------------------------------------------
# Vizualize
set.seed(123)
sdr_plot_df <- 
  tibble(sdr_k1 = terra::values(sdr_corr, mat = F)) %>% 
  drop_na() %>% 
  slice_sample(n = nrow(loss_bath_sdr)) %>% 
  bind_cols(loss_bath_sdr) %>% 
  dplyr::select(-tot_loss, -sed_rate) %>% 
  gather(method, SDR)

# Plot
sdr_comparison_plot <-
  sdr_plot_df %>% 
  ggplot(aes(y = method,
             x = SDR)) +
  stat_halfeye(aes(fill = stat(cut_cdf_qi(
    cdf, 
    .width = c(.66, .8, .95),
    labels = scales::percent_format()
  )))) +
  scale_fill_brewer(direction = -1,
                    na.translate = FALSE) +
  scale_x_continuous(breaks = seq(0.00, 0.14, by = 0.02),
                     limits = c(0.01, 0.14)) +
  scale_y_discrete(labels = parse(text = c("SDR[Manual]",
                                           "SDR[Cavalli]"))) +
  labs(x = "Sediment Delivery Ratio (SDR)",
       y = "Method",
       fill = "CI") +
  theme(legend.position = c(0.9, 0.8),
        legend.direction = "vertical")

ggsave("R/06_sdr-validation/figures",
       sdr_comparison_plot,
       dpi = 500,
       w = 7,
       h = 5)


# 4) Are the SDR significantly different or not? --------------------------
# Wilcox test
sdr_plot_df %>% 
  wilcox.test(SDR ~ method,
              data = .)

# Student t-test
set.seed(123)
rr <- sample(seq_len(250000), 2000)

sdr_sampled <- 
  sdr_plot_df %>% 
  group_by(method) %>% 
  mutate(id = seq_len(250000)) %>% 
  ungroup() %>% 
  filter(id %in% rr) %>% 
  dplyr::select(-id)

tt <- 
  t.test(SDR ~ method,
         data = sdr_sampled,
         var.equal = FALSE,
         alternative = "less")

p.adjust(tt$p.value, method = "bonferroni")  

# save --------------------------------------------------------------------
ggsave(
  sdr_comparison_plot, 
  filename = here("R", "06_sdr-validation", "figures",
                  "SDR.tiff"),
  device = ragg::agg_png,
  dpi = 1000,
  width = 16, 
  height = 9, 
  units = "cm"
)
