source("script/loadpackages.R")

# region = read_rds('data/inputs/shp/GAUL.rds') %>%
# st_simplify(dTolerance = 0.1)

get_growing_season <- function(akey) {
  dates <- map_dfr(c("rf"), function(atype) {
    planting_day <- raster(str_c("data/inputs/crop_calendar/", akey, "_", atype, "_ggcmi_crop_calendar_phase3_v1.01.nc4"), varname = "planting_day")

    growing_season_length <- raster(str_c("data/inputs/crop_calendar/", akey, "_", atype, "_ggcmi_crop_calendar_phase3_v1.01.nc4"), varname = "growing_season_length")

    if (akey %in% names(wheat_sup)) {
      planting_day <-
        planting_day %>%
        mask(wheat_sup[[akey]])

      growing_season_length <-
        growing_season_length %>%
        mask(wheat_sup[[akey]])
    }

    if (akey %in% c("ri1", "ri2")) {
      rice_mask <- raster(str_c("data/inputs/crop_calendar/", akey, "_", atype, "_ggcmi_crop_calendar_phase3_v1.01.nc4"), varname = "fraction_of_harvested_area")

      rice_mask[rice_mask <= 0.2] <- NA

      planting_day <-
        planting_day %>%
        mask(rice_mask)

      growing_season_length <-
        growing_season_length %>%
        mask(rice_mask)
    }

    planting_day <-
      planting_day %>%
      as.data.frame(xy = T) %>%
      as_tibble() %>%
      drop_na() %>%
      `names<-`(c("x", "y", "Planting"))

    growing_season_length <-
      growing_season_length %>%
      as.data.frame(xy = T) %>%
      as_tibble() %>%
      drop_na() %>%
      `names<-`(c("x", "y", "Growing_season_length"))

    dates <- list(planting_day, growing_season_length) %>%
      reduce(left_join, by = join_by(x, y)) %>%
      mutate(
        Planting = as.Date(Planting, origin = "0000-01-01"),
        Haversting = Planting %m+% days(Growing_season_length),
        Planting = if_else(day(Planting) <= 15, floor_date(Planting, unit = "month"), ceiling_date(Planting, unit = "month")),
        Haversting = if_else(day(Haversting) >= 15, ceiling_date(Haversting, unit = "month"), floor_date(Haversting, unit = "month"))
      ) %>%
      drop_na() %>%
      select(x, y, contains("Planting"), contains("Haversting"))

    return(dates)
  })

  dates <- dates %>%
    mutate(
      YOM = pmap(list(Planting, Haversting), function(ap, ah) {
        seq(ap, ah, by = "months") %>% `names<-`(1:length(.))
      }),
      offset = map_dbl(YOM, ~ max(year(.x)))
    ) %>%
    select(x, y, YOM, offset) %>%
    unnest_longer(YOM) %>%
    mutate(YOM_id = as.numeric(YOM_id))

  # how below works
  # expand_grid(head = seq(0, 10, by = 2) * 0.1, tail = seq(0, 10, by = 2) * 0.1) %>%
  #   filter(head < tail) %>%
  #   mutate(data = map2(head, tail, function(ahead, atail){
  #
  #     tibble(a = 1:10) %>%
  #       filter(between(a, quantile(a, ahead), quantile(a, atail)))
  #
  #   }))

  template <- expand_grid(head = seq(0, 10, by = 2) * 0.1, tail = seq(0, 10, by = 2) * 0.1) %>%
    filter(head < tail) %>%
    mutate(prop = str_c(head * 10, tail * 10, sep = "_")) # * 10 must be done for the later raster processing

  # 1 mins 6 cores
  dates <- template %>%
    mutate(data = future_map2(head, tail, function(ahead, atail) {
      dates %>%
        filter(between(YOM_id, quantile(YOM_id, ahead), quantile(YOM_id, atail)), .by = c(x, y))
    }, .progress = F), .keep = "unused") %>%
    unnest(data)

  dates <- dates %>%
    select(-YOM_id) %>%
    funique()

  growing_season <-
    tibble(
      Index_year = 2003:2019,
      data = future_map(Index_year, ~ mutate(dates, YOM = YOM %m+% years(. - offset)))
    ) %>%
    unnest(data) %>%
    mutate(year = year(YOM), month = month(YOM), .keep = "unused") %>%
    select(-offset)

  return(growing_season)
}

wheat_sup <- rast("data/inputs/crop_calendar/wheat_supplement/winter_and_spring_wheat_areas_phase3.nc4")["mask"] %>%
  `names<-`(c("wwh", "swh")) %>%
  stack()

wheat_sup[wheat_sup == 0] <- NA

plan(multisession, workers = 6)

a <- tibble(
  Item =
    c("Barley", "Bean", "Cassava", "Cotton", "Groundnut", "Maize", "Millet", "Pea", "Potato", "Rye", "Rapeseed", "Sorghum", "Soybean", "Sugarbeet", "Sugarcane", "Sunflower", "Wheat", "Wheat", "Rice1", "Rice2"),
  key = c("bar", "bea", "cas", "cot", "nut", "mai", "mil", "pea", "pot", "rye", "rap", "sor", "soy", "sgb", "sgc", "sun", "swh", "wwh", "ri1", "ri2")
) %>%
  mutate(growing_season = map(key, get_growing_season, .progress = T))

plan(sequential)

# for wheat and rice because of two seasons
a <- bind_rows(
  a %>%
    select(-key) %>%
    filter(!Item %in% c("Wheat")),
  a %>%
    select(-key) %>%
    filter(Item %in% c("Wheat")) %>%
    unnest() %>%
    funique() %>%
    nest(growing_season = c(Index_year, prop, x, y, year, month))
) %>%
  mutate(growing_season = map(growing_season, function(ags) {
    ags %>%
      nest(.by = prop, .key = "fdata")
  }))

qsave(a, "data/outputs/gs/growing_season.qs", nthreads = 5)
