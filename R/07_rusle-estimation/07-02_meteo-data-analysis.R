library(dplyr)
library(purrr)
library(tidyr)
library(readxl)
library(atslib)
library(lubridate)
library(sf)
library(AOI)
library(terra)
library(here)
library(PNWColors)
library(performance)
library(writexl)
library(ggfx)
library(ggrepel)
library(ggpmisc)
library(stringr)

ggplot2::theme_set(atslib::theme_hp())

databoard <- pins::board_folder("data/tidy", versioned = TRUE)

# 1) Load data from meteostations -----------------------------------------

# NCDC
load(here::here("data", "meteo", "ncdc_monthly.Rdata"))

# Terskol meteostation — data from colleagues
terskol_toropov <-
  readxl::read_xlsx(
    here::here("data", "meteo", "toropov_TERSKOL.xlsx"),
    sheet = 2
  ) %>%
  dplyr::rename(date = 1, temp = 2, prec = 3) %>%
  dplyr::mutate(date = lubridate::as_date(date)) %>%
  dplyr::mutate(dplyr::across(
    dplyr::where(is.numeric),
    ~ dplyr::na_if(.x, -999.9)
  ))

# Get Terskol's coords
terskol_sf <-
  AOI::geocode(c("Terskol"), pt = TRUE) %>%
  dplyr::rename(name = request)

# Measured rainfall erosivity (from Larionov et al., 1993)
epo <-
  readxl::read_excel(here::here(
    "data",
    "meteo",
    "ЭПО_прикавказье.xlsx"
  )) %>%
  magrittr::set_colnames(
    #fmt:skip
    c("st_id", "station_name",
      "y", "x", "period", "EPO")
  ) %>%
  sf::st_as_sf(coords = c("x", "y"), crs = 4326) %>%
  dplyr::mutate(EPO = as.numeric(EPO))

# Calculate AOI extent with buffer
ext2 <-
  epo %>%
  sf::st_bbox() %>%
  sf::st_as_sfc() %>%
  sf::st_transform(6931) %>%
  sf::st_buffer(10^4) %>%
  sf::st_transform(4326) %>%
  terra::vect()

# Save
epo %>%
  dplyr::filter(
    !st_id %in%
      c(
        28,
        121,
        # 116,
        113
        # 91
      )
  ) %>%
  dplyr::mutate(
    station_name = stringi::stri_trans_general(
      station_name,
      "russian-latin/bgn"
    ),
    station_name = ifelse(
      stringr::str_detect(station_name, "Pyatigorsk"),
      "Pyatigorsk",
      station_name
    )
  ) %>%
  sf::st_write(
    "data/spatial/gizh_viz.gpkg",
    layer = "meteostations",
    delete_layer = T
  )


# 2) Read TerraClimate dataset --------------------------------------------
# Read NetCDF files and stack them
precipitation_rasters <-
  list.files(
    path = here::here("data", "terraclim"),
    pattern = ".ppt.",
    full.names = T
  ) %>%
  purrr::map(terra::rast)

# Crop rasters and subset to summer months
precipitation_rasters_sub <-
  precipitation_rasters %>%
  purrr::map(~ terra::crop(.x, ext2)) %>%
  purrr::map(~ terra::subset(.x, 5:10))

# Calculate Modified Fouriner Index
precipitation_mfi <-
  purrr::map(
    .x = precipitation_rasters_sub,
    ~ terra::app(.x, fun = function(x) {
      ss <- sum(x)
      sum(x^2 / ss)
    })
  ) %>%
  purrr::set_names(seq(from = 1958, to = 2020, by = 1))

# Filter and keep years to which measured rainfall erosivity is
# available (according to Larionov et al., 1993)
precipitation_mfi_larionov <-
  precipitation_mfi %>%
  purrr::keep(terra::names(precipitation_mfi) %in% c(1961:1983)) %>%
  terra::rast() %>%
  terra::median()

# Create buffer around meteostations
stations_buf <- epo %>%
  sf::st_buffer(dist = 5000)

# Extract MFI values
stations_mfi <-
  stations_buf %>%
  dplyr::mutate(
    mfi = terra::extract(
      precipitation_mfi_larionov,
      terra::vect(stations_buf),
      median
    )[, 2]
  ) %>%
  sf::st_drop_geometry() %>%
  dplyr::filter(
    !st_id %in%
      c(
        28,
        121,
        # 116,
        113
        # 91
      )
  ) %>%
  dplyr::mutate(EPO = EPO * 58.8) #  convert to MJ mm / ha h yr

# formula <- y ~ poly(x, 2, raw = T)
formula <- y ~ x

mfi_mod <-
  stations_mfi %>%
  stats::lm(EPO ~ mfi, data = .)

epo_mfi_plot <-
  stations_mfi %>%
  cbind(
    terra::predict(mfi_mod, stations_mfi, interval = "prediction") %>%
      terra::as.data.frame() %>%
      dplyr::rename_with(~ paste0(., "_pred"))
  ) %>%
  cbind(
    terra::predict(mfi_mod, stations_mfi, interval = "confidence") %>%
      terra::as.data.frame() %>%
      dplyr::rename_with(~ paste0(., "_CI"))
  ) %>%
  # rename stations
  dplyr::mutate(
    station_name = stringi::stri_trans_general(
      station_name,
      "russian-latin/bgn"
    ),
    station_name = ifelse(
      stringr::str_detect(station_name, "Pyatigorsk"),
      "Pyatigorsk",
      station_name
    )
  ) %>%
  ggplot2::ggplot(ggplot2::aes(x = mfi, y = EPO)) +
  # Prediction
  ggplot2::geom_line(ggplot2::aes(x = mfi, y = fit_pred)) +
  # Prediciton interval
  ggplot2::geom_line(
    ggplot2::aes(y = lwr_pred, color = "Prediction Interval"),
    linetype = "dashed"
  ) +
  ggplot2::geom_line(
    ggplot2::aes(y = upr_pred, color = "Prediction Interval"),
    linetype = "dashed"
  ) +
  # CI interval
  ggplot2::geom_line(
    ggplot2::aes(y = lwr_CI, color = "95% Confidence Interval"),
    linetype = "dashed"
  ) +
  ggplot2::geom_line(
    ggplot2::aes(y = upr_CI, color = "95% Confidence Interval"),
    linetype = "dashed"
  ) +
  # Points
  ggplot2::geom_point(alpha = .5) +
  ggrepel::geom_text_repel(
    ggplot2::aes(label = station_name),
    show.legend = F,
    family = "Roboto Condensed",
    size = 2.5
  ) +
  ggpmisc::stat_poly_eq(
    ggplot2::aes(
      label = paste(
        ggplot2::after_stat(eq.label),
        ggplot2::after_stat(rr.label),
        ggplot2::after_stat(p.value.label),
        sep = "*\", \"*"
      )
    ),
    eq.x.rhs = "italic(MFI)",
    eq.with.lhs = "italic(R)~`=`~",
    small.p = T,
    formula = formula
  ) +
  ggplot2::scale_x_continuous(breaks = seq(50, 120, 10)) +
  ggplot2::scale_color_manual(values = c("#0099cc", "#ff3030")) +
  ggplot2::labs(
    x = "Modified Fournier Index (MFI)",
    y = "Rainfall erosivity (R)",
    color = NULL
  ) +
  ggplot2::theme(
    legend.position = "inside",
    legend.position.inside = c(0.20, 0.80)
  )

# 3) Compare measured rainfall and MFI with TerraClimate ------------------
# Get TerraClim values for NCDC meteostations
rain_tc <-
  precipitation_rasters %>%
  purrr::map(~ terra::extract(.x, terra::vect(ncdc_sf), mean)) %>%
  purrr::map(
    ~ dplyr::mutate(
      .x,
      ID = terra::values(terra::vect(ncdc_sf))$id,
      name = terra::values(terra::vect(ncdc_sf))$name,
      .after = ID
    )
  ) %>%
  purrr::set_names(seq(from = 1958, to = 2020, by = 1)) %>%
  dplyr::bind_rows(.id = "year") %>%
  tidyr::gather(month, rain_tc, -year, -ID, -name) %>%
  dplyr::mutate(month = stringr::str_remove(month, "ppt_")) %>%
  dplyr::transmute(
    station = ID,
    name,
    date = lubridate::make_date(year, month, 1),
    rain_tc
  ) %>%
  dplyr::arrange(date) %>%
  tidyr::as_tibble()

# Extract TerraClimate values for Terskol meteostation
rain_tc_terskol <-
  precipitation_rasters %>%
  purrr::map(~ terra::extract(.x, terra::vect(terskol_sf), mean)) %>%
  purrr::map(
    ~ dplyr::mutate(
      .x,
      name = terra::values(terra::vect(terskol_sf))$name,
      .before = ppt_1
    )
  ) %>%
  purrr::set_names(seq(from = 1958, to = 2020, by = 1)) %>%
  dplyr::bind_rows(.id = "year") %>%
  tidyr::gather(month, rain_tc, -year, -ID, -name) %>%
  dplyr::mutate(month = stringr::str_remove(month, "ppt_")) %>%
  dplyr::transmute(
    station = ID,
    name,
    date = lubridate::make_date(year, month, 1),
    rain_tc
  ) %>%
  dplyr::arrange(date) %>%
  tidyr::as_tibble()

# Join Terskol's measured and reanalysis
terskol_toropov_tc <-
  terskol_toropov %>%
  dplyr::group_by(
    year = lubridate::year(date),
    month = lubridate::month(date)
  ) %>%
  dplyr::summarise(
    precip = sum(prec, na.rm = T),
    .groups = "drop"
  ) %>%
  dplyr::transmute(
    date = lubridate::make_date(year, month, 1),
    precip
  ) %>%
  dplyr::left_join(rain_tc_terskol, by = "date") %>%
  dplyr::select(-station)

# Join NCDC, Terskol and TerraClim datasets
ncdc_tc_terskol <-
  caucasus_ncdc %>%
  dplyr::left_join(rain_tc, by = c("station", "date")) %>%
  dplyr::select(-station) %>%
  dplyr::bind_rows(terskol_toropov_tc) %>%
  dplyr::mutate(
    name = stringr::str_remove(name, ", RS| AMSG, RS"),
    name = stringr::str_to_upper(name)
  )

# Compare NCDC and TerraClimate
ncdc_tc_plot <-
  ncdc_tc_terskol %>%
  ggplot2::ggplot(ggplot2::aes(x = rain_tc, y = precip)) +
  ggplot2::geom_point(ggplot2::aes(color = name)) +
  Add_1_line() +
  ggplot2::geom_smooth(
    method = "lm",
    color = "dimgrey",
    linetype = "dashed",
    se = F,
    formula = y ~ x
  ) +
  ggpmisc::stat_poly_eq(
    ggplot2::aes(
      label = paste(
        ggplot2::after_stat(rr.label),
        ggplot2::after_stat(n.label),
        sep = "*\", \"*"
      )
    ),
    formula = y ~ x
  ) +
  ggplot2::scale_color_manual(
    values = PNWColors::pnw_palette("Bay", 6)
  ) +
  ggplot2::facet_wrap(~name, scales = "free") +
  ggplot2::labs(
    color = "Meteostation",
    y = "Measured precipitation, mm",
    x = "Terra Climate precipitation, mm"
  )

ncdc_tc_terskol %>%
  dplyr::group_by(name) %>%
  tidyr::nest() %>%
  dplyr::mutate(
    model = purrr::map(data, ~ stats::lm(precip ~ rain_tc, data = .x))
  ) %>%
  dplyr::mutate(
    prfm = purrr::map(model, ~ performance::model_performance(.x))
  ) %>%
  tidyr::unnest("prfm") %>%
  dplyr::ungroup() %>%
  dplyr::summarise(
    r2_mean = terra::mean(R2),
    r2_min = min(R2),
    r2_max = max(R2),
    rmse_mean = terra::mean(RMSE)
  )

ncdc_tc_table <-
  ncdc_tc_terskol %>%
  dplyr::group_by(name) %>%
  tidyr::nest() %>%
  dplyr::mutate(
    mindate = purrr::map_chr(data, ~ as.character(min(.x$date))),
    maxdate = purrr::map_chr(data, ~ as.character(max(.x$date))),
    n = purrr::map_int(data, ~ terra::nrow(.x))
  ) %>%
  dplyr::mutate(
    model = purrr::map(data, ~ stats::lm(precip ~ rain_tc, data = .x))
  ) %>%
  dplyr::mutate(
    prfm = purrr::map(
      model,
      ~ performance::model_performance(.x, metrics = c("R2", "RMSE"))
    )
  ) %>%
  tidyr::unnest("prfm") %>%
  dplyr::ungroup() %>%
  dplyr::left_join(
    ncdc_sf %>%
      dplyr::mutate(
        name = stringr::str_remove(name, ", RS| AMSG, RS")
      ) %>%
      dplyr::mutate(
        longitude = sf::st_coordinates(.)[, 1],
        latitude = sf::st_coordinates(.)[, 2]
      ) %>%
      sf::st_drop_geometry() %>%
      dplyr::select(name, elevation, longitude, latitude),
    by = "name"
  ) %>%
  dplyr::select(-data, -model) %>%
  dplyr::select(name, elevation:latitude, mindate:n, R2:RMSE)


# 4) Calculate erosivity time-series for Gizhgit --------------------------
ws <-
  sf::st_read(
    "data/vector/watershed/giz_2022_ws-dem.shp",
    quiet = TRUE
  ) %>%
  sf::st_transform(terra::crs(precipitation_mfi[[1]], proj = T))

giz_mfi <-
  precipitation_mfi %>%
  purrr::map(~ terra::extract(.x, terra::vect(ws), median)) %>%
  dplyr::bind_rows(.id = "year") %>%
  dplyr::transmute(year = as.integer(year), mfi = lyr.1) %>%
  tidyr::as_tibble()

giz_r_timeseries <-
  giz_mfi %>%
  # 0.5 confidence level
  dplyr::bind_cols(
    terra::predict(
      mfi_mod,
      .,
      interval = "confidence",
      level = 0.5
    ) %>%
      terra::as.data.frame() %>%
      dplyr::rename_all(~ paste0(., "_50"))
  ) %>%
  # 0.8 confidence level
  dplyr::bind_cols(
    terra::predict(
      mfi_mod,
      .,
      interval = "confidence",
      level = 0.8
    ) %>%
      terra::as.data.frame() %>%
      dplyr::rename_all(~ paste0(., "_80"))
  ) %>%
  # 0.95
  dplyr::bind_cols(
    terra::predict(
      mfi_mod,
      .,
      interval = "confidence",
      level = 0.95
    ) %>%
      terra::as.data.frame() %>%
      dplyr::rename_all(~ paste0(., "_95"))
  ) %>%
  dplyr::select(-fit_50, -fit_80) %>%
  tidyr::gather(conflevel, bounds, -fit_95, -mfi, -year) %>%
  dplyr::mutate(
    type = ifelse(
      stringr::str_detect(conflevel, "lwr"),
      "lwr",
      "upr"
    ),
    conflevel = stringr::str_extract(conflevel, "\\d+")
  ) %>%
  tidyr::pivot_wider(names_from = type, values_from = bounds) %>%
  dplyr::mutate(
    conflevel = as.numeric(conflevel),
    conflevel = conflevel / 100,
    conflevel = forcats::as_factor(conflevel),
    conflevel = forcats::fct_rev(conflevel)
  ) %>%
  ggplot2::ggplot(ggplot2::aes(x = year, y = fit_95)) +
  ggplot2::geom_ribbon(
    ggplot2::aes(
      ymin = lwr,
      ymax = upr,
      group = conflevel,
      alpha = conflevel
    ),
    fill = "#8fa1cb"
  ) +
  ggplot2::geom_line(color = "#7670b2") +
  ggplot2::geom_smooth(
    ggplot2::aes(color = "Trend line"),
    method = "lm",
    linetype = "f8",
    lwd = .7,
    se = F
  ) +
  ggplot2::scale_x_continuous(
    breaks = seq(1958, 2020, by = 10),
    minor_breaks = seq(1958, 2020, by = 1),
    expand = c(0.01, 0)
  ) +
  ggplot2::scale_y_continuous(
    breaks = seq(400, 1400, by = 200),
    expand = c(0.01, 0)
  ) +
  ggplot2::scale_color_manual(values = "#0099cc", name = NULL) +
  ggplot2::labs(
    y = "Rainfall erosivity (R)",
    x = "",
    alpha = "Confidence level"
  )

# Rainfall Erosivity for Gizhgit cathcment during
# various years
giz_r <-
  giz_mfi %>%
  dplyr::bind_cols(
    terra::predict(
      mfi_mod,
      .,
      interval = "confidence",
      level = 0.95
    ) %>%
      terra::as.data.frame()
  )

# Explore temporal trend in Rainfall Erosivity
giz_r_year_mod <- lm(fit ~ year, data = giz_r)
100 * coefficients(giz_r_year_mod)[2] / mean(giz_r$fit)
Kendall::MannKendall(giz_r$fit)

# 5) Calculate mean rainfall erosivity values for periods -----------------
mfi_for_period <-
  function(rasters, periods) {
    period_string <- paste0(min(periods), "-", max(periods))

    rasters %>%
      purrr::keep(terra::names(rasters) %in% periods) %>%
      terra::rast() %>%
      terra::median() %>%
      terra::extract(terra::vect(ws), median) %>%
      dplyr::transmute(period = period_string, mfi = median) %>%
      dplyr::bind_cols(
        terra::predict(mfi_mod, ., interval = "confidence") %>%
          terra::as.data.frame()
      )
  }

periods <- list(
  c(1961L:1983L),
  c(1961L:1989L),
  c(1989L:2000L),
  c(2000L:2020L),
  c(1989L:2020L)
)

giz_erosivity <- purrr::map_dfr(
  periods,
  ~ mfi_for_period(precipitation_mfi, .x)
) %>%
  tidyr::as_tibble()

r_table <- giz_erosivity %>%
  dplyr::rename(R = fit) %>%
  dplyr::rename_all(~ stringr::str_to_sentence(.)) %>%
  dplyr::mutate_if(is.numeric, ~ atslib::smart_round(.))

# Panagos' Rainfall Erosivity --------------------------------------------
glored <-
  ws |>
  sf::st_transform(6931) |>
  sf::st_buffer(500) |>
  sf::st_transform(4326) |> 
  terra::vect() |>
  rusleR::get_glored(warp = TRUE)

library(exactextractr)

# Get GLORED summary stats
exactextractr::exact_extract(glored, ws, fun = "mean")
exactextractr::exact_extract(glored, ws, "quantile", quantiles = c(0.025, 0.975))

# Save --------------------------------------------------------------------
save("giz_erosivity", file = here("data", "tidy", "giz_r.Rdata"))

pins::pin_write(
  board = databoard,
  x = giz_erosivity,
  name = "giz_r",
  type = "csv",
  description = stringr::str_squish(
    "Rainfall Erosivity (in SI values) predicted for the 
    Gitche-Gizhgit Catchment using the empirical relationship
    between Modified Fournier Index and measured 
    Rainfall Intensity during the 1961:1983 period.
    
    Listed mean and lower/upper confidence intervals
    derived from the linear model:
      R = -65.9 + 8.82*MFI"
  )
)

pins::pin_write(
  board = databoard,
  x = ncdc_tc_table,
  name = "ncdc_tc_table",
  type = "csv",
  description = stringr::str_squish(
    "TerraClimate quality assessment of the daily rainfall sums
     for the region"
  )
)

# Save figures and tables ------------------------------------------------
fs::dir_create("figures/supplementary")

ggsave(
  filename = here("figures", "supplementary", "r-mfi_plt.png"),
  epo_mfi_plot,
  dpi = 500,
  w = 7,
  h = 4.5
)

# Figure 4.1
ggsave(
  filename = here("figures", "r_timeseries.png"),
  giz_r_timeseries,
  dpi = 500,
  w = 8,
  h = 5
)

# Figure S1
ggsave(
  filename = here("figures", "supplementary", "ncdc_tc_plt.png"),
  ncdc_tc_plot,
  dpi = 500,
  w = 9,
  h = 8
)
