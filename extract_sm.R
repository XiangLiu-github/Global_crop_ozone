source("script/loadpackages.R")
source("script/loadfunctions.R")

# https://essd.copernicus.org/articles/14/4473/2022/essd-14-4473-2022.html

temp <- rast("data/inputs/crop_calendar/mai_rf_ggcmi_crop_calendar_phase3_v1.01.nc4")[[1]]

data <- map_dfr(2003:2019, function(ayear) {
  datap <- zip::zip_list(str_c("../data_archive/SGD-SM2.0/", ayear, ".zip")) %>%
    filter(compressed_size != 0) %>%
    select(filename) %>%
    mutate(
      date = str_remove(filename, as.character(ayear)),
      date = str_remove(date, "/SGD_SM_") %>% ymd(),
      year = year(date), month = month(date), day = day(date),
      filename = str_c("/vsizip/../data_archive/SGD-SM2.0/", ayear, ".zip/", filename),
      data = map(filename, function(afile) {
        rast(afile)["reconstructed_sm"] %>%
          `ext<-`(ext(temp)) %>%
          `names<-`("sm") %>%
          resample(temp) %>%
          as.data.frame(xy = T)
      }, .progress = F)
    ) %>%
    select(-date, -filename) %>%
    unnest()

  datap %>%
    fgroup_by(x, y, year, month) %>%
    fsummarise(sm = fmean(sm), sm2 = fmean(sm^2), sm3 = fmean(sm^3))
}, .progress = T)

qsave(data, "data/outputs/sm/sm.qs")

library(qs)
library(raster)
library(tidyverse)
library(exactextractr)

data <- qread("data/outputs/sm/sm.qs")

gs <- qread("data/outputs/gs/growing_season.qs", nthreads = 5)

library(dtplyr)
gs_sm <- gs %>%
  mutate(sm = map(growing_season, function(ags) {
    gc()

    ags %>%
      mutate(fdata = map(fdata, function(afdata) {
        left_join(afdata, data, by = join_by(x, y, year, month)) %>%
          lazy_dt() %>%
          group_by(Index_year, x, y) %>%
          filter(sum(is.na(sm)) == 0) %>%
          summarise(across(starts_with("sm"), mean),
            .groups = "drop"
          ) %>%
          as_tibble()
      })) %>%
      unnest(fdata) %>%
      pivot_wider(
        names_from = c(Index_year, prop),
        names_sep = ".",
        values_from = starts_with("sm")
      )

    # left_join(ags, data, by = join_by(x, y, year, month)) %>%
    #   lazy_dt() %>%
    #   group_by(Index_year, prop, x, y) %>%
    #   filter(sum(is.na(sm)) == 0) %>%
    #   summarise(across(starts_with('sm'), mean),
    #             .groups = 'drop') %>%
    #   pivot_wider(names_from = c(Index_year, prop),
    #               names_sep = ".",
    #               values_from = starts_with('sm')) %>%
    #   as_tibble()
  }, .progress = T)) %>%
  select(-growing_season)

gc()

gs_sm <- bind_rows(
  gs_sm %>%
    filter(!str_detect(Item, "Rice")),
  gs_sm %>%
    filter(str_detect(Item, "Rice")) %>%
    mutate(Item = "Rice") %>%
    unnest() %>%
    lazy_dt() %>%
    summarise(across(everything(), mean), .by = c(Item, x, y)) %>%
    as_tibble() %>%
    nest(.by = Item, .key = "sm")
) %>%
  mutate(sm = map(sm, ~ rasterFromXYZ(.x, crs = "+proj=longlat +datum=WGS84 +no_defs")))

qsave(gs_sm, "data/outputs/sm/gs_sm.qs", nthreads = 5)

gs_sm <- qread("data/outputs/sm/gs_sm.qs", nthreads = 5)

# plan(multisession, workers = 3)

region <- read_rds("data/inputs/shp/GAUL.rds")

weights <- qread("data/inputs/cropland/weights.qs")

extracted <- gs_sm %>%
  inner_join(weights) %>%
  mutate(sm = pmap(list(sm, weights), function(abrick, aweight) {
    exact_extract(abrick, region, weights = aweight, fun = "weighted_mean", default_weight = 0, stack_apply = T, progress = F) %>%
      as_tibble() %>%
      bind_cols(region %>% sf::st_drop_geometry()) %>%
      pivot_longer(-c(ADM_CODE, ADM_NAME, FAO_CODE),
        names_to = c("sm", "year", "prop"), names_sep = "\\.",
        names_transform = list(year = as.integer),
        names_prefix = "weighted_mean.", values_drop_na = T
      ) %>%
      pivot_wider(names_from = sm)
  }, .progress = T)) %>%
  select(-weights)

saveRDS(extracted, "data/outputs/sm/extracted_sm.rds")

# interact irrigation -----------------------------------------------------

order <- get_season() %>%
  filter(var == "sm") %>%
  select(Item, prop)

gs_sm <- qread("data/outputs/sm/gs_sm.qs", nthreads = 5)

gs_sm <- inner_join(order, gs_sm) %>%
  mutate(sm = map2(prop, sm, function(aprop, asm) {
    asm %>%
      subset(names(asm)[str_detect(names(asm), aprop)])
  })) %>%
  select(-prop)

# plan(multisession, workers = 3)

region <- read_rds("data/inputs/shp/GAUL.rds")

weights <- qread("data/inputs/cropland/weights.qs")

irr <- qread("data/inputs/GMIA/irr.qs")

extracted <- gs_sm %>%
  inner_join(weights) %>%
  inner_join(irr) %>%
  mutate(sm = pmap(list(sm, weights, irr_frc), function(abrick, aweight, airr_frc) {
    abrick <- (abrick * airr_frc) %>%
      `names<-`(names(abrick))

    exact_extract(abrick, region, weights = aweight, fun = "weighted_mean", default_weight = 0, stack_apply = T, progress = F) %>%
      as_tibble() %>%
      bind_cols(region %>% sf::st_drop_geometry()) %>%
      pivot_longer(-c(ADM_CODE, ADM_NAME, FAO_CODE),
        names_to = c("sm", "year", "prop"), names_sep = "\\.",
        names_transform = list(year = as.integer),
        names_prefix = "weighted_mean.", values_drop_na = T
      ) %>%
      pivot_wider(names_from = sm)
  }, .progress = T)) %>%
  select(-c(weights, irr_frc))

saveRDS(extracted, "data/outputs/sm/extracted_sm_itr_irr.rds")
