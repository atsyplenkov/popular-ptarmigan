#FIXME:
# Check that all paths are working

library(tidyverse)
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
library(stringi)

theme_set(theme_hp())

# 1) Load data from meteostations -----------------------------------------

# NCDC
load(here("data", "meteo", "ncdc_monthly.Rdata"))

# Terskol meteostation — data from colleagues
terskol_toropov <- read_xlsx(
  here("data", "meteo", "toropov_TERSKOL.xlsx"),
  sheet = 2
) %>%
  rename(date = 1, temp = 2, prec = 3) %>%
  mutate(date = as_date(date)) %>%
  mutate(across(where(is.numeric), ~ na_if(.x, -999.9)))

# Get Terskol's coords
terskol_sf <- AOI::geocode(c("Terskol"), pt = TRUE) %>%
  rename(name = request)

# Measured rainfall erosivity (from Larionov et al., 1993)
epo <- read_excel(here("data", "meteo", "ЭПО_прикавказье.xlsx")) %>%
  magrittr::set_colnames(
    #fmt:skip
    c("st_id", "station_name",
      "y", "x", "period", "EPO")
  ) %>%
  st_as_sf(coords = c("x", "y"), crs = 4326) %>%
  mutate(EPO = as.numeric(EPO))

# Calculate AOI extent
ext2 <- epo %>%
  st_bbox() %>%
  st_as_sfc() %>%
  st_transform(6931) %>%
  st_buffer(10^4) %>%
  st_transform(4326) %>%
  vect()

# Save
epo %>%
  filter(
    !st_id %in%
      c(
        28,
        121,
        # 116,
        113
        # 91
      )
  ) %>%
  mutate(
    station_name = stri_trans_general(
      station_name,
      "russian-latin/bgn"
    ),
    station_name = ifelse(
      str_detect(station_name, "Pyatigorsk"),
      "Pyatigorsk",
      station_name
    )
  ) %>%
  st_write(
    "data/spatial/gizh_viz.gpkg",
    layer = "meteostations",
    delete_layer = T
  )


# 2) Read TerraClimate dataset --------------------------------------------
# Read NetCDF files and stack them
precipitation_rasters <-
  list.files(
    path = here("data", "terraclim"),
    pattern = ".ppt.",
    full.names = T
  ) %>%
  map(rast)

# Crop rasters and subset to summer months
precipitation_rasters_sub <-
  precipitation_rasters %>%
  map(~ crop(.x, ext2)) %>%
  map(~ subset(.x, 5:10))

# Calculate Modified Fouriner Index
precipitation_mfi <- map(
  .x = precipitation_rasters_sub,
  ~ app(.x, fun = function(x) {
    ss <- sum(x)
    sum(x^2 / ss)
  })
) %>%
  set_names(seq(from = 1958, to = 2020, by = 1))

# Filter and keep years to which measured rainfall erosivity is
# available (according to Larionov et al., 1993)
precipitation_mfi_larionov <-
  precipitation_mfi %>%
  keep(names(precipitation_mfi) %in% c(1961:1983)) %>%
  rast() %>%
  median()

# Create buffer around meteostations
stations_buf <- epo %>%
  st_buffer(dist = 5000)

# Extract MFI values
stations_mfi <-
  stations_buf %>%
  mutate(
    mfi = extract(
      precipitation_mfi_larionov,
      vect(stations_buf),
      median
    )[, 2]
  ) %>%
  st_drop_geometry() %>%
  filter(
    !st_id %in%
      c(
        28,
        121,
        # 116,
        113
        # 91
      )
  ) %>%
  mutate(EPO = EPO * 58.8) #  convert to MJ mm / ha h yr

# formula <- y ~ poly(x, 2, raw = T)
formula <- y ~ x

mfi_mod <-
  stations_mfi %>%
  lm(EPO ~ mfi, data = .)

epo_mfi_plot <-
  stations_mfi %>%
  cbind(
    predict(mfi_mod, stations_mfi, interval = "prediction") %>%
      as.data.frame() %>%
      rename_with(~ paste0(., "_pred"))
  ) %>%
  cbind(
    predict(mfi_mod, stations_mfi, interval = "confidence") %>%
      as.data.frame() %>%
      rename_with(~ paste0(., "_CI"))
  ) %>%
  # rename stations
  mutate(
    station_name = stri_trans_general(
      station_name,
      "russian-latin/bgn"
    ),
    station_name = ifelse(
      str_detect(station_name, "Pyatigorsk"),
      "Pyatigorsk",
      station_name
    )
  ) %>%
  ggplot(aes(x = mfi, y = EPO)) +
  # Prediction
  geom_line(aes(x = mfi, y = fit_pred)) +
  # Prediciton interval
  geom_line(
    aes(y = lwr_pred, color = "Prediction Interval"),
    linetype = "dashed"
  ) +
  geom_line(
    aes(y = upr_pred, color = "Prediction Interval"),
    linetype = "dashed"
  ) +
  # CI interval
  geom_line(
    aes(y = lwr_CI, color = "95% Confidence Interval"),
    linetype = "dashed"
  ) +
  geom_line(
    aes(y = upr_CI, color = "95% Confidence Interval"),
    linetype = "dashed"
  ) +
  # Points
  geom_point(alpha = .5) +
  geom_text_repel(
    aes(label = station_name),
    show.legend = F,
    family = "Roboto Condensed",
    size = 2.5
  ) +
  stat_poly_eq(
    aes(
      label = paste(
        after_stat(eq.label),
        after_stat(rr.label),
        after_stat(p.value.label),
        sep = "*\", \"*"
      )
    ),
    eq.x.rhs = "italic(MFI)",
    eq.with.lhs = "italic(R)~`=`~",
    small.p = T,
    formula = formula
  ) +
  scale_x_continuous(breaks = seq(50, 120, 10)) +
  scale_color_manual(values = c("#0099cc", "#ff3030")) +
  labs(
    x = "Modified Fournier Index (MFI)",
    y = "Rainfall erosivity (R)",
    color = NULL
  ) +
  theme(legend.position = c(0.20, 0.80))

# 3) Compare measured rainfall and MFI with TerraClimate ------------------
# Get TerraClim values for NCDC meteostations
rain_tc <-
  precipitation_rasters %>%
  map(~ extract(.x, vect(ncdc_sf), mean)) %>%
  map(
    ~ mutate(
      .x,
      ID = values(vect(ncdc_sf))$id,
      name = values(vect(ncdc_sf))$name,
      .after = ID
    )
  ) %>%
  set_names(seq(from = 1958, to = 2020, by = 1)) %>%
  bind_rows(.id = "year") %>%
  gather(month, rain_tc, -year, -ID, -name) %>%
  mutate(month = str_remove(month, "ppt_")) %>%
  transmute(
    station = ID,
    name,
    date = make_date(year, month, 1),
    rain_tc
  ) %>%
  arrange(date) %>%
  as_tibble()

# Extract TerraClimate values for Terskol meteostation
rain_tc_terskol <-
  precipitation_rasters %>%
  map(~ extract(.x, vect(terskol_sf), mean)) %>%
  map(
    ~ mutate(
      .x,
      name = values(vect(terskol_sf))$name,
      .before = ppt_1
    )
  ) %>%
  set_names(seq(from = 1958, to = 2020, by = 1)) %>%
  bind_rows(.id = "year") %>%
  gather(month, rain_tc, -year, -ID, -name) %>%
  mutate(month = str_remove(month, "ppt_")) %>%
  transmute(
    station = ID,
    name,
    date = make_date(year, month, 1),
    rain_tc
  ) %>%
  arrange(date) %>%
  as_tibble()

# Join Terskol's measured and reanalysis
terskol_toropov_tc <-
  terskol_toropov %>%
  group_by(year = year(date), month = month(date)) %>%
  summarise(precip = sum(prec, na.rm = T), .groups = "drop") %>%
  transmute(date = make_date(year, month, 1), precip) %>%
  left_join(rain_tc_terskol, by = "date") %>%
  dplyr::select(-station)

# Join NCDC, Terskol and TerraClim datasets
ncdc_tc_terskol <-
  caucasus_ncdc %>%
  left_join(rain_tc, by = c("station", "date")) %>%
  dplyr::select(-station) %>%
  bind_rows(terskol_toropov_tc) %>%
  mutate(
    name = str_remove(name, ", RS| AMSG, RS"),
    name = str_to_upper(name)
  )

# Compare NCDC and TerraClimate
ncdc_tc_plot <-
  ncdc_tc_terskol %>%
  ggplot(aes(x = rain_tc, y = precip)) +
  geom_point(aes(color = name)) +
  Add_1_line() +
  geom_smooth(
    method = "lm",
    color = "dimgrey",
    linetype = "dashed",
    se = F,
    formula = y ~ x
  ) +
  stat_poly_eq(
    aes(
      label = paste(
        after_stat(rr.label),
        after_stat(n.label),
        sep = "*\", \"*"
      )
    ),
    formula = y ~ x
  ) +
  scale_color_manual(values = pnw_palette("Bay", 6)) +
  facet_wrap(~name, scales = "free") +
  labs(
    color = "Meteostation",
    y = "Measured precipitation, mm",
    x = "Terra Climate precipitation, mm"
  )

ncdc_tc_terskol %>%
  group_by(name) %>%
  nest() %>%
  mutate(model = map(data, ~ lm(precip ~ rain_tc, data = .x))) %>%
  mutate(prfm = map(model, ~ performance::model_performance(.x))) %>%
  unnest("prfm") %>%
  ungroup() %>%
  summarise(
    r2_mean = mean(R2),
    r2_min = min(R2),
    r2_max = max(R2),
    rmse_mean = mean(RMSE)
  )

ncdc_tc_table <-
  ncdc_tc_terskol %>%
  group_by(name) %>%
  nest() %>%
  mutate(
    mindate = map_chr(data, ~ as.character(min(.x$date))),
    maxdate = map_chr(data, ~ as.character(max(.x$date))),
    n = map_int(data, ~ nrow(.x))
  ) %>%
  mutate(model = map(data, ~ lm(precip ~ rain_tc, data = .x))) %>%
  mutate(
    prfm = map(
      model,
      ~ performance::model_performance(.x, metrics = c("R2", "RMSE"))
    )
  ) %>%
  unnest("prfm") %>%
  ungroup() %>%
  left_join(
    ncdc_sf %>%
      mutate(name = str_remove(name, ", RS| AMSG, RS")) %>%
      mutate(
        longitude = sf::st_coordinates(.)[, 1],
        latitude = sf::st_coordinates(.)[, 2]
      ) %>%
      st_drop_geometry() %>%
      dplyr::select(name, elevation, longitude, latitude),
    by = "name"
  ) %>%
  dplyr::select(-data, -model) %>%
  dplyr::select(name, elevation:latitude, mindate:n, R2:RMSE)

# 4) Calculate erosivity time-series for Gizhgit --------------------------
ws <- st_read("data/spatial/gizh_ws-topo.shp") %>%
  st_transform(crs(precipitation_mfi[[1]], proj = T))

giz_mfi <-
  precipitation_mfi %>%
  map(~ extract(.x, vect(ws), median)) %>%
  bind_rows(.id = "year") %>%
  transmute(year = as.integer(year), mfi = lyr.1) %>%
  as_tibble()

giz_r_timeseries <- giz_mfi %>%
  # 0.5 confidence level
  bind_cols(
    predict(mfi_mod, ., interval = "confidence", level = 0.5) %>%
      as.data.frame() %>%
      rename_all(~ paste0(., "_50"))
  ) %>%
  # 0.8 confidence level
  bind_cols(
    predict(mfi_mod, ., interval = "confidence", level = 0.8) %>%
      as.data.frame() %>%
      rename_all(~ paste0(., "_80"))
  ) %>%
  # 0.95
  bind_cols(
    predict(mfi_mod, ., interval = "confidence", level = 0.95) %>%
      as.data.frame() %>%
      rename_all(~ paste0(., "_95"))
  ) %>%
  dplyr::select(-fit_50, -fit_80) %>%
  gather(conflevel, bounds, -fit_95, -mfi, -year) %>%
  mutate(
    type = ifelse(str_detect(conflevel, "lwr"), "lwr", "upr"),
    conflevel = str_extract(conflevel, "\\d+")
  ) %>%
  pivot_wider(names_from = type, values_from = bounds) %>%
  mutate(
    conflevel = as.numeric(conflevel),
    conflevel = conflevel / 100,
    conflevel = as_factor(conflevel),
    conflevel = fct_rev(conflevel)
  ) %>%
  ggplot(aes(x = year, y = fit_95)) +
  geom_ribbon(
    aes(ymin = lwr, ymax = upr, group = conflevel, alpha = conflevel),
    fill = "#8fa1cb"
  ) +
  geom_line(color = "#7670b2") +
  geom_smooth(
    aes(color = "Trend line"),
    method = "lm",
    linetype = "f8",
    size = .7,
    se = F
  ) +
  scale_x_continuous(
    breaks = seq(1958, 2020, by = 10),
    minor_breaks = seq(1958, 2020, by = 1),
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(
    breaks = seq(400, 1400, by = 200),
    expand = c(0.01, 0)
  ) +
  scale_color_manual(values = "#0099cc", name = NULL) +
  labs(
    y = "Rainfall erosivity (R)",
    x = "",
    alpha = "Confidence level"
  )

giz_r <- giz_mfi %>%
  # 0.5 confidence level
  bind_cols(
    predict(mfi_mod, ., interval = "confidence", level = 0.95) %>%
      as.data.frame()
  )

giz_r_year_mod <- lm(fit ~ year, data = giz_r)

100 * coefficients(giz_r_year_mod)[2] / mean(giz_r$fit)

Kendall::MannKendall(giz_r$fit)

# 5) Calculate mean rainfall erosivity values for periods -----------------

mfi_for_period <- function(rasters, periods) {
  period_string <- paste0(min(periods), "-", max(periods))

  rasters %>%
    keep(names(rasters) %in% periods) %>%
    rast() %>%
    median() %>%
    extract(vect(ws), median) %>%
    transmute(period = period_string, mfi = median) %>%
    bind_cols(
      predict(mfi_mod, ., interval = "confidence") %>%
        as.data.frame()
    )
}

periods <- list(
  c(1961L:1983L),
  c(1961L:1989L),
  c(1989L:2000L),
  c(2000L:2020L),
  c(1989L:2020L)
)

giz_erosivity <- map_dfr(
  periods,
  ~ mfi_for_period(precipitation_mfi, .x)
) %>%
  as_tibble()

r_table <- giz_erosivity %>%
  rename(R = fit) %>%
  rename_all(~ str_to_sentence(.)) %>%
  mutate_if(is.numeric, ~ smart_round(.))

# Save --------------------------------------------------------------------
save("giz_erosivity", file = here("data", "tidy", "giz_r.Rdata"))


# Figure 2
ggsave(
  filename = here("figures", "fig02_r-mfi_plt.png"),
  epo_mfi_plot,
  dpi = 500,
  w = 7,
  h = 4.5
)

# Figure 4.1
ggsave(
  filename = here("figures", "fig04-1_r_timeseries.png"),
  giz_r_timeseries,
  dpi = 500,
  w = 8,
  h = 5
)

# Figure S1
ggsave(
  filename = here("figures", "figs1_ncdc_tc_plt.png"),
  ncdc_tc_plot,
  dpi = 500,
  w = 9,
  h = 8
)

# Table S1
write_xlsx(ncdc_tc_table, path = here("tables", "table_s1.xlsx"))

# Table 4.1
write_xlsx(r_table, path = here("tables", "table_4-1.xlsx"))
