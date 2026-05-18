library(sf); library(here); library(dplyr); library(tigris); library(lubridate)
sf_use_s2(FALSE)
here::i_am("analysis/03_emapr_biomass_exploration.qmd")

mtbs <- sf::st_read(here("data","raw","mtbs","mtbs_perimeter_data","mtbs_perims_DD.shp"), quiet=TRUE)
ca   <- tigris::states(cb=TRUE, year=2022, resolution="5m") |> dplyr::filter(STUSPS=="CA")

mtbs_ca <- mtbs |>
  dplyr::filter(lubridate::year(ig_date) %in% c(2000, 2001, 2003)) |>
  sf::st_transform(sf::st_crs(ca)) |>
  sf::st_filter(ca)

mtbs_ca |>
  dplyr::mutate(year = lubridate::year(ig_date)) |>
  dplyr::group_by(year) |>
  dplyr::summarise(
    n_fires      = n(),
    total_acres  = round(sum(burnbndac, na.rm=TRUE)),
    median_acres = round(median(burnbndac, na.rm=TRUE)),
    max_acres    = round(max(burnbndac, na.rm=TRUE))
  ) |>
  dplyr::arrange(year) |>
  as.data.frame() |>
  print()
