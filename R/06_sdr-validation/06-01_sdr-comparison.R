library(tidyverse)
library(here)
library(qs)
library(atslib)
library(terra)
library(ggdist)

ggplot2::theme_set(atslib::theme_hp())

# 1) Load data from previous sections -------------------------------------
# All possible sedimentation rates, t/yr
sed_rates <- qs::qread("R/03_lake-bathymetry/data/sed_rates.qs")
loss_mc <-
  qs::qread(
    here::here("R", "05_erosion-mapping", "data", "loss_t-yr.qs")
  )

# Watershed
ws <-
  terra::vect(
    here::here("data", "vector", "watershed", "giz_2022_ws-dem.shp")
  )

# Lake shoreline (polygon)
lake <-
  terra::vect(
    here::here("data", "vector", "shoreline_dem", "gizh_bl_dem.shp")
  )

# Basin DEM
dem <-
  terra::rast(
    here::here(
      "data",
      "raster",
      "basin_dem",
      "giz_2022_dem_r5-f-lake.tif"
    )
  )

# SDR
sdr_corr <-
  #FIXME:
  # replace SDR
  terra::rast(
    here::here(
      "R",
      "04_connectivity",
      "data",
      "sdr_rockfall_corr.tiff"
    )
  ) %>%
  terra::resample(dem, "cubic") %>%
  terra::crop(ws, mask = T) %>%
  terra::mask(lake, inverse = T)

# 2) Estimate SDR based on processes map ----------------------------------
# Select 500 random samples
set.seed(123)
sed_r_sample <- sample(seq_len(1000), 500)

# SDR computed as total loss divided by sedimentation
loss_bath_sdr <-
  loss_mc %>%
  tibble::rownames_to_column(var = "id") %>%
  dplyr::as_tibble() %>%
  dplyr::rename(tot_loss = sum) %>%
  # TEMPORARY!!!
  dplyr::mutate(sed_rate = sed_rates[sed_r_sample]) %>%
  # dplyr::mutate(sed_rate = sed_rates) %>%
  tidyr::expand(tot_loss, sed_rate) %>%
  dplyr::mutate(sdr = sed_rate / tot_loss)

# SDR stats
my_sdr_means <-
  loss_bath_sdr %>%
  ggdist::mean_qi(sdr, na.rm = T) %>%
  dplyr::mutate(dplyr::across(
    dplyr::where(is.numeric),
    ~ terra::round(.x, 3)
  ))

# 3) Compare SDRs ---------------------------------------------------------
# Vizualize
set.seed(123)
sdr_plot_df <-
  dplyr::tibble(sdr_k1 = terra::values(sdr_corr, mat = F)) %>%
  tidyr::drop_na() %>%
  dplyr::slice_sample(n = terra::nrow(loss_bath_sdr)) %>%
  dplyr::bind_cols(loss_bath_sdr) %>%
  dplyr::select(-tot_loss, -sed_rate) %>%
  tidyr::gather(method, SDR)

# Plot
sdr_comparison_plot <-
  sdr_plot_df %>%
  ggplot2::ggplot(ggplot2::aes(y = method, x = SDR)) +
  ggdist::stat_halfeye(ggplot2::aes(
    fill = ggplot2::after_stat(ggdist::cut_cdf_qi(
      cdf,
      .width = c(.66, .8, .95),
      labels = scales::percent_format()
    ))
  )) +
  ggplot2::scale_fill_brewer(direction = -1, na.translate = FALSE) +
  ggplot2::scale_x_continuous(
    breaks = seq(0.00, 0.14, by = 0.02),
    limits = c(0.01, 0.14)
  ) +
  ggplot2::scale_y_discrete(
    labels = parse(text = c("SDR[Manual]", "SDR[Cavalli]"))
  ) +
  ggplot2::labs(
    x = "Sediment Delivery Ratio (SDR)",
    y = "Method",
    fill = "CI"
  ) +
  ggplot2::theme(
    legend.position = "inside",
    legend.position.inside = c(0.9, 0.8),
    legend.direction = "vertical"
  )

# 4) Are the SDR significantly different or not? --------------------------
# Wilcox test
stats::wilcox.test(SDR ~ method, data = sdr_plot_df)

# Student t-test
set.seed(123)
rr <- sample(seq_len(250000), 2000)

sdr_sampled <-
  sdr_plot_df %>%
  dplyr::group_by(method) %>%
  dplyr::mutate(id = seq_len(250000)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(id %in% rr) %>%
  dplyr::select(-id)

tt <-
  stats::t.test(
    SDR ~ method,
    data = sdr_sampled,
    var.equal = FALSE,
    alternative = "less"
  )

stats::p.adjust(tt$p.value, method = "bonferroni")

# save --------------------------------------------------------------------
ggplot2::ggsave(
  sdr_comparison_plot,
  filename = here::here(
    "R",
    "06_sdr-validation",
    "figures",
    "SDR.tiff"
  ),
  device = ragg::agg_png,
  dpi = 1000,
  width = 16,
  height = 9,
  units = "cm"
)
