library(sf)
library(dplyr)

ws_lulc <-
  sf::st_read("data/vector/lulc-ortho/giz_2022_lulc-ortho_25may.shp")


ws_lulc |>
  dplyr::mutate(
    lulc = dplyr::case_when(
      id == 1 ~ "shrub",
      id == 2 ~ "shrub",
      id == 3 ~ "forest",
      id == 4 ~ "natural grassland"
    )
  ) %>%
  mutate(area = as.numeric(st_area(x = .))) |>
  st_drop_geometry() |>
  filter(!is.na(lulc)) |>
  group_by(lulc) |>
  reframe(area = sum(area)) |>
  mutate(area_pct = 100 * area / sum(area))
