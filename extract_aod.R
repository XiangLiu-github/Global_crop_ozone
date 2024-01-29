source("script/loadpackages.R")
source("script/loadfunctions.R")

temp <- rast("data/inputs/crop_calendar/mai_rf_ggcmi_crop_calendar_phase3_v1.01.nc4")[[1]]

data <- map_dfr(2003:2019, function(ayear) {
  datap <- tibble(
    files = list.files("../data_archive/MERRA2_daily_AOD/", full.names = T, pattern = str_c("^", ayear, ".*")),
    files_ = list.files("../data_archive/MERRA2_daily_AOD/", pattern = str_c("^", ayear, ".*"))
  ) %>%
    mutate(
      date = str_remove(files_, ".tif") %>% ymd(),
      data = map(files, function(afile) {
        rast(afile) %>%
          resample(temp) %>%
          mask(temp) %>%
          `names<-`("AOD") %>%
          as.data.frame(xy = T)
      }, .progress = F),
      year = year(date),
      month = month(date),
      .keep = "unused"
    ) %>%
    unnest()

  datap %>%
    fgroup_by(x, y, year, month) %>%
    fsummarise(AOD = fmean(AOD), AOD2 = fmean(AOD^2), AOD3 = fmean(AOD^3))
}, .progress = T)

qsave(data, "data/outputs/aod/aod.qs")

library(qs)
library(raster)
library(tidyverse)
library(exactextractr)

data <- qread("data/outputs/aod/aod.qs")

gs <- qread("data/outputs/gs/growing_season.qs", nthreads = 5)

library(dtplyr)
gs_aod <- gs %>%
  mutate(aod = map(growing_season, function(ags) {
    gc()

    ags %>%
      mutate(fdata = map(fdata, function(afdata) {
        left_join(afdata, data, by = join_by(x, y, year, month)) %>%
          lazy_dt() %>%
          group_by(Index_year, x, y) %>%
          filter(sum(is.na(AOD)) == 0) %>%
          summarise(across(starts_with("AOD"), mean),
            .groups = "drop"
          ) %>%
          as_tibble()
      })) %>%
      unnest(fdata) %>%
      pivot_wider(
        names_from = c(Index_year, prop),
        names_sep = ".",
        values_from = starts_with("AOD")
      )

    # left_join(ags, data, by = join_by(x, y, year, month)) %>%
    #   lazy_dt() %>%
    #   group_by(Index_year, prop, x, y) %>%
    #   filter(sum(is.na(AOD)) == 0) %>%
    #   summarise(across(starts_with('AOD'), mean),
    #             .groups = 'drop') %>%
    #   pivot_wider(names_from = c(Index_year, prop),
    #               names_sep = ".",
    #               values_from = starts_with('AOD')) %>%
    #   as_tibble()
  }, .progress = T)) %>%
  select(-growing_season)

gs_aod <- bind_rows(
  gs_aod %>%
    filter(!str_detect(Item, "Rice")),
  gs_aod %>%
    filter(str_detect(Item, "Rice")) %>%
    mutate(Item = "Rice") %>%
    unnest() %>%
    lazy_dt() %>%
    summarise(across(everything(), mean), .by = c(Item, x, y)) %>%
    as_tibble() %>%
    nest(.by = Item, .key = "aod")
) %>%
  mutate(aod = map(aod, ~ rasterFromXYZ(.x, crs = "+proj=longlat +datum=WGS84 +no_defs")))

qsave(gs_aod, "data/outputs/aod/gs_aod.qs", nthreads = 5)

gs_aod <- qread("data/outputs/aod/gs_aod.qs", nthreads = 5)

region <- read_rds("data/inputs/shp/GAUL.rds")

weights <- qread("data/inputs/cropland/weights.qs")

extracted <- gs_aod %>%
  inner_join(weights) %>%
  mutate(aod = pmap(list(aod, weights), function(abrick, aweight) {
    exact_extract(abrick, region, weights = aweight, fun = "weighted_mean", default_weight = 0, stack_apply = T, progress = F) %>%
      as_tibble() %>%
      bind_cols(region %>% sf::st_drop_geometry()) %>%
      pivot_longer(-c(ADM_CODE, ADM_NAME, FAO_CODE),
        names_to = c("AOD", "year", "prop"), names_sep = "\\.",
        names_transform = list(year = as.integer),
        names_prefix = "weighted_mean.", values_drop_na = T
      ) %>%
      pivot_wider(names_from = AOD)
  }, .progress = T)) %>%
  select(-weights)

# plan(sequential)

saveRDS(extracted, "data/outputs/aod/extracted_aod.rds")

# interact irrigation -----------------------------------------------------

order <- get_season() %>%
  filter(var == "AOD") %>%
  select(Item, prop)

gs_aod <- qread("data/outputs/aod/gs_aod.qs", nthreads = 5)

gs_aod <- inner_join(order, gs_aod) %>%
  mutate(aod = map2(prop, aod, function(aprop, aaod) {
    aaod %>%
      subset(names(aaod)[str_detect(names(aaod), aprop)])
  })) %>%
  select(-prop)

# plan(multisession, workers = 3)

region <- read_rds("data/inputs/shp/GAUL.rds")

weights <- qread("data/inputs/cropland/weights.qs")

irr <- qread("data/inputs/GMIA/irr.qs")

extracted <- gs_aod %>%
  inner_join(weights) %>%
  inner_join(irr) %>%
  mutate(aod = pmap(list(aod, weights, irr_frc), function(abrick, aweight, airr_frc) {
    abrick <- (abrick * airr_frc) %>%
      `names<-`(names(abrick))

    exact_extract(abrick, region, weights = aweight, fun = "weighted_mean", default_weight = 0, stack_apply = T, progress = F) %>%
      as_tibble() %>%
      bind_cols(region %>% sf::st_drop_geometry()) %>%
      pivot_longer(-c(ADM_CODE, ADM_NAME, FAO_CODE),
        names_to = c("aod", "year", "prop"), names_sep = "\\.",
        names_transform = list(year = as.integer),
        names_prefix = "weighted_mean.", values_drop_na = T
      ) %>%
      pivot_wider(names_from = aod)
  }, .progress = T)) %>%
  select(-c(weights, irr_frc))

saveRDS(extracted, "data/outputs/aod/extracted_aod_itr_irr.rds")
