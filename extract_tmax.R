source("script/loadpackages.R")
source("script/loadfunctions.R")

############################################
temp <- rast("data/inputs/crop_calendar/mai_rf_ggcmi_crop_calendar_phase3_v1.01.nc4")[[1]]

data <- map_dfr(2003:2019, function(ayear) {
  rst <- rast(str_c("../data_archive/CPC_daily_climate/tmax.", ayear, ".nc"))

  map_dfr(1:nlyr(rst), function(index) {
    rst[[index]] %>%
      rotate() %>%
      resample(temp) %>%
      `names<-`("tmax") %>%
      as.data.frame(xy = T) %>%
      mutate(date = time(rst[[index]]))
  }) %>%
    mutate(year = year(date), month = month(date), .keep = "unused") %>%
    fgroup_by(x, y, year, month) %>%
    fsummarise(tmax = fmean(tmax), tmax2 = fmean(tmax^2), tmax3 = fmean(tmax^3))
}, .progress = T)

qsave(data, "data/outputs/tmax/tmax.qs")

library(qs)
library(raster)
library(tidyverse)
library(exactextractr)

data <- qread("data/outputs/tmax/tmax.qs")

gs <- qread("data/outputs/gs/growing_season.qs", nthreads = 5)

library(dtplyr)
gs_tmax <- gs %>%
  mutate(tmax = map(growing_season, function(ags) {
    gc()

    ags %>%
      mutate(fdata = map(fdata, function(afdata) {
        left_join(afdata, data, by = join_by(x, y, year, month)) %>%
          lazy_dt() %>%
          group_by(Index_year, x, y) %>%
          filter(sum(is.na(tmax)) == 0) %>%
          summarise(across(starts_with("tmax"), mean),
            .groups = "drop"
          ) %>%
          as_tibble()
      })) %>%
      unnest(fdata) %>%
      pivot_wider(
        names_from = c(Index_year, prop),
        names_sep = ".",
        values_from = starts_with("tmax")
      )

    # left_join(ags, data, by = join_by(x, y, year, month)) %>%
    #   lazy_dt() %>%
    #   group_by(Index_year, prop, x, y) %>%
    #   filter(sum(is.na(tmax)) == 0) %>%
    #   summarise(across(starts_with('tmax'), mean),
    #             .groups = 'drop') %>%
    #   pivot_wider(names_from = c(Index_year, prop),
    #               names_sep = ".",
    #               values_from = starts_with('tmax')) %>%
    #   as_tibble()
  }, .progress = T)) %>%
  select(-growing_season)

gc()

gs_tmax <- bind_rows(
  gs_tmax %>%
    filter(!str_detect(Item, "Rice")),
  gs_tmax %>%
    filter(str_detect(Item, "Rice")) %>%
    mutate(Item = "Rice") %>%
    unnest() %>%
    lazy_dt() %>%
    summarise(across(everything(), mean), .by = c(Item, x, y)) %>%
    as_tibble() %>%
    nest(.by = Item, .key = "tmax")
) %>%
  mutate(tmax = map(tmax, ~ rasterFromXYZ(.x, crs = "+proj=longlat +datum=WGS84 +no_defs")))

qsave(gs_tmax, "data/outputs/tmax/gs_tmax.qs", nthreads = 5)

gs_tmax <- qread("data/outputs/tmax/gs_tmax.qs", nthreads = 5)

region <- read_rds("data/inputs/shp/GAUL.rds")

weights <- qread("data/inputs/cropland/weights.qs")

# plan(multisession, workers = 3)

extracted <- gs_tmax %>%
  inner_join(weights) %>%
  mutate(tmax = pmap(list(tmax, weights), function(abrick, aweight) {
    exact_extract(abrick, region, weights = aweight, fun = "weighted_mean", default_weight = 0, stack_apply = T, progress = F) %>%
      as_tibble() %>%
      bind_cols(region %>% sf::st_drop_geometry()) %>%
      pivot_longer(-c(ADM_CODE, ADM_NAME, FAO_CODE),
        names_to = c("tmax", "year", "prop"), names_sep = "\\.",
        names_transform = list(year = as.integer),
        names_prefix = "weighted_mean.", values_drop_na = T
      ) %>%
      pivot_wider(names_from = tmax)
  }, .progress = T)) %>%
  select(-weights)

# plan(sequential)

saveRDS(extracted, "data/outputs/tmax/extracted_tmax.rds")


# counterfactual 2degree --------------------------------------------------

data <- map_dfr(2003:2019, function(ayear) {
  rst <- rast(str_c("../data_archive/CPC_daily_climate/tmax.", ayear, ".nc"))

  map_dfr(1:nlyr(rst), function(index) {
    rst[[index]] %>%
      rotate() %>%
      resample(temp) %>%
      `names<-`("tmax") %>%
      as.data.frame(xy = T) %>%
      mutate(date = time(rst[[index]]))
  }) %>%
    fmutate(year = year(date), month = month(date), day = day(date), .keep = "unused", tmax = tmax + 2) %>%
    fgroup_by(x, y, year, month) %>%
    fsummarise(tmax = fmean(tmax), tmax2 = fmean(tmax^2), tmax3 = fmean(tmax^3))
}, .progress = T)

qsave(data, "data/outputs/tmax/tmax_count2.qs")

data <- qread("data/outputs/tmax/tmax_count2.qs")

order <- get_season() %>%
  filter(var == "tmax") %>%
  select(Item, prop)

order <- bind_rows(
  order %>% filter(Item != "Rice"),
  order %>%
    filter(Item == "Rice") %>%
    mutate(Item = map(Item, ~ str_c("Rice", 1:2))) %>%
    unnest()
)

gs <- qread("data/outputs/gs/growing_season.qs", nthreads = 5) %>%
  inner_join(order) %>%
  mutate(growing_season = map2(growing_season, prop, function(ags, aprop) {
    gc()

    ags <- ags %>%
      filter(prop == aprop) %>%
      pull()

    ags[[1]]
  }, .progress = T))

library(dtplyr)
gs_tmax <- gs %>%
  mutate(tmax = map2(growing_season, prop, function(ags, aprop) {
    left_join(ags, data, by = join_by(x, y, year, month)) %>%
      lazy_dt() %>%
      group_by(Index_year, x, y) %>%
      filter(sum(is.na(tmax)) == 0) %>%
      summarise(across(starts_with("tmax"), mean),
        .groups = "drop"
      ) %>%
      pivot_wider(
        names_from = Index_year,
        names_sep = ".",
        values_from = starts_with("tmax")
      ) %>%
      as_tibble()
  }, .progress = T)) %>%
  select(-growing_season)

gs_tmax <- bind_rows(
  gs_tmax %>%
    filter(!str_detect(Item, "Rice")),
  gs_tmax %>%
    filter(str_detect(Item, "Rice")) %>%
    mutate(Item = "Rice") %>%
    unnest() %>%
    lazy_dt() %>%
    summarise(across(everything(), mean), .by = c(Item, prop, x, y)) %>%
    as_tibble() %>%
    nest(.by = c(Item, prop), .key = "tmax")
) %>%
  mutate(tmax = map(tmax, ~ rasterFromXYZ(.x, crs = "+proj=longlat +datum=WGS84 +no_defs")))

qsave(gs_tmax, "data/outputs/tmax/gs_tmax_count2.qs", nthreads = 5)

gs_tmax <- qread("data/outputs/tmax/gs_tmax_count2.qs", nthreads = 5)

region <- read_rds("data/inputs/shp/GAUL.rds")

weights <- qread("data/inputs/cropland/weights.qs")

# plan(multisession, workers = 3)

extracted <- gs_tmax %>%
  inner_join(weights) %>%
  mutate(tmax = pmap(list(tmax, weights), function(abrick, aweight) {
    exact_extract(abrick, region, weights = aweight, fun = "weighted_mean", default_weight = 0, stack_apply = T, progress = F) %>%
      as_tibble() %>%
      bind_cols(region %>% st_drop_geometry()) %>%
      pivot_longer(-c(ADM_CODE, ADM_NAME, FAO_CODE),
        names_to = c("tmax", "year"), names_sep = "\\.",
        names_transform = list(year = as.integer),
        names_prefix = "weighted_mean.", values_drop_na = T
      ) %>%
      pivot_wider(names_from = tmax)
  }, .progress = T)) %>%
  select(-weights)

# plan(sequential)

saveRDS(extracted, "data/outputs/tmax/extracted_tmax_count2.rds")


# counterfactual 2degree whole --------------------------------------------------

order <- get_season() %>%
  filter(var == "tmax") %>%
  select(Item, prop) %>%
  mutate(prop = "0_10")

order <- bind_rows(
  order %>% filter(Item != "Rice"),
  order %>%
    filter(Item == "Rice") %>%
    mutate(Item = map(Item, ~ str_c("Rice", 1:2))) %>%
    unnest()
)

data <- qread("data/outputs/tmax/tmax_count2.qs")

gs <- qread("data/outputs/gs/growing_season.qs", nthreads = 5) %>%
  inner_join(order) %>%
  mutate(growing_season = map2(growing_season, prop, function(ags, aprop) {
    gc()

    ags <- ags %>%
      filter(prop == aprop) %>%
      pull()

    ags[[1]]
  }, .progress = T))

library(dtplyr)
gs_tmax <- gs %>%
  mutate(tmax = map2(growing_season, prop, function(ags, aprop) {
    left_join(ags, data, by = join_by(x, y, year, month)) %>%
      lazy_dt() %>%
      group_by(Index_year, x, y) %>%
      filter(sum(is.na(tmax)) == 0) %>%
      summarise(across(starts_with("tmax"), mean),
        .groups = "drop"
      ) %>%
      pivot_wider(
        names_from = Index_year,
        names_sep = ".",
        values_from = starts_with("tmax")
      ) %>%
      as_tibble()
  }, .progress = T)) %>%
  select(-growing_season)

gs_tmax <- bind_rows(
  gs_tmax %>%
    filter(!str_detect(Item, "Rice")),
  gs_tmax %>%
    filter(str_detect(Item, "Rice")) %>%
    mutate(Item = "Rice") %>%
    unnest() %>%
    lazy_dt() %>%
    summarise(across(everything(), mean), .by = c(Item, prop, x, y)) %>%
    as_tibble() %>%
    nest(.by = c(Item, prop), .key = "tmax")
) %>%
  mutate(tmax = map(tmax, ~ rasterFromXYZ(.x, crs = "+proj=longlat +datum=WGS84 +no_defs")))

qsave(gs_tmax, "data/outputs/tmax/gs_tmax_count2_whole.qs", nthreads = 5)

gs_tmax <- qread("data/outputs/tmax/gs_tmax_count2_whole.qs", nthreads = 5)

region <- read_rds("data/inputs/shp/GAUL.rds")

weights <- qread("data/inputs/cropland/weights.qs")

# plan(multisession, workers = 3)

extracted <- gs_tmax %>%
  inner_join(weights) %>%
  mutate(tmax = pmap(list(tmax, weights), function(abrick, aweight) {
    exact_extract(abrick, region, weights = aweight, fun = "weighted_mean", default_weight = 0, stack_apply = T, progress = F) %>%
      as_tibble() %>%
      bind_cols(region %>% st_drop_geometry()) %>%
      pivot_longer(-c(ADM_CODE, ADM_NAME, FAO_CODE),
        names_to = c("tmax", "year"), names_sep = "\\.",
        names_transform = list(year = as.integer),
        names_prefix = "weighted_mean.", values_drop_na = T
      ) %>%
      pivot_wider(names_from = tmax)
  }, .progress = T)) %>%
  select(-weights)

# plan(sequential)

saveRDS(extracted, "data/outputs/tmax/extracted_tmax_count2_whole.rds")


# interact irrigation -----------------------------------------------------

order <- get_season() %>%
  filter(var == "tmax") %>%
  select(Item, prop)

gs_tmax <- qread("data/outputs/tmax/gs_tmax.qs", nthreads = 5)

gs_tmax <- inner_join(order, gs_tmax) %>%
  mutate(tmax = map2(prop, tmax, function(aprop, atmax) {
    atmax %>%
      subset(names(atmax)[str_detect(names(atmax), aprop)])
  })) %>%
  select(-prop)

region <- read_rds("data/inputs/shp/GAUL.rds")

weights <- qread("data/inputs/cropland/weights.qs")

irr <- qread("data/inputs/GMIA/irr.qs")

# plan(multisession, workers = 3)

extracted <- gs_tmax %>%
  inner_join(weights) %>%
  inner_join(irr) %>%
  mutate(tmax = pmap(list(tmax, weights, irr_frc), function(abrick, aweight, airr_frc) {
    abrick <- (abrick * airr_frc) %>%
      `names<-`(names(abrick))

    exact_extract(abrick, region, weights = aweight, fun = "weighted_mean", default_weight = 0, stack_apply = T, progress = F) %>%
      as_tibble() %>%
      bind_cols(region %>% sf::st_drop_geometry()) %>%
      pivot_longer(-c(ADM_CODE, ADM_NAME, FAO_CODE),
        names_to = c("tmax", "year", "prop"), names_sep = "\\.",
        names_transform = list(year = as.integer),
        names_prefix = "weighted_mean.", values_drop_na = T
      ) %>%
      pivot_wider(names_from = tmax)
  }, .progress = T)) %>%
  select(-c(weights, irr_frc))

# plan(sequential)

saveRDS(extracted, "data/outputs/tmax/extracted_tmax_itr_irr.rds")


# future scenario ----------------------------------------------------

data <- tibble(
  model = c("GFDL-ESM4_r1i1p1f1", "MPI-ESM1-2-HR_r1i1p1f1", "MRI-ESM2-0_r1i1p1f1", "NorESM2-LM_r1i1p1f1", "NorESM2-MM_r1i1p1f1", "UKESM1-0-LL_r1i1p1f2"),
  data = map(model, function(amodel) {
    arast <- rast(str_c("D:/Data/CMIP6/remapped/tasmax_", amodel, "_ssp370.nc")) %>%
      `names<-`(time(.))

    map_dfr(c(2044:2050, 2094:2100), function(ayear) {
      arast[[year(time(arast)) == ayear]] %>%
        as.data.frame(xy = T) %>%
        pivot_longer(-c(x, y),
          names_to = c("year", "month", "day"), names_sep = "-", names_transform = list(year = as.integer, month = as.integer, day = as.integer),
          values_to = c("tmax")
        ) %>%
        fgroup_by(x, y, year, month) %>%
        fsummarise(tmax = fmean(tmax), tmax2 = fmean(tmax^2), tmax3 = fmean(tmax^3))
    }) %>%
      mutate(
        indx = fifelse(year <= 2050, "2045-2049", "2095-2099"),
        year = fifelse(year <= 2050, year - 40, year - 90)
      )
  }, .progress = T)
)

data <- data %>%
  unnest() %>%
  nest(.by = c(model, indx))

qsave(data, "data/outputs/tmax/tmax_future.qs")

data <- qread("data/outputs/tmax/tmax_future.qs")

order <- get_season() %>%
  filter(var == "tmax") %>%
  select(Item, prop)

order <- bind_rows(
  order %>% filter(Item != "Rice"),
  order %>%
    filter(Item == "Rice") %>%
    mutate(Item = map(Item, ~ str_c("Rice", 1:2))) %>%
    unnest()
)

gs <- qread("data/outputs/gs/growing_season.qs") %>%
  inner_join(order) %>%
  mutate(growing_season = map2(growing_season, prop, function(ags, aprop) {
    gc()

    ags <- ags %>%
      filter(prop == aprop) %>%
      pull()

    ags[[1]]
  }, .progress = T))

region <- read_rds("data/inputs/shp/GAUL.rds")

weights <- qread("data/inputs/cropland/weights.qs")

library(dtplyr)
data <- data %>%
  mutate(data = map(data, function(adata) {
    gs_tmax <- gs %>%
      mutate(tmax = map2(growing_season, prop, function(ags, aprop) {
        left_join(ags, adata, by = join_by(x, y, year, month)) %>%
          lazy_dt() %>%
          group_by(Index_year, x, y) %>%
          filter(sum(is.na(tmax)) == 0) %>%
          summarise(across(starts_with("tmax"), mean),
            .groups = "drop"
          ) %>%
          pivot_wider(
            names_from = Index_year,
            names_sep = ".",
            values_from = starts_with("tmax")
          ) %>%
          as_tibble()
      }, .progress = F)) %>%
      select(-growing_season)

    gs_tmax <- bind_rows(
      gs_tmax %>%
        filter(!str_detect(Item, "Rice")),
      gs_tmax %>%
        filter(str_detect(Item, "Rice")) %>%
        mutate(Item = "Rice") %>%
        unnest() %>%
        lazy_dt() %>%
        summarise(across(everything(), mean), .by = c(Item, prop, x, y)) %>%
        as_tibble() %>%
        nest(.by = c(Item, prop), .key = "tmax")
    ) %>%
      mutate(tmax = map(tmax, ~ rasterFromXYZ(.x, crs = "+proj=longlat +datum=WGS84 +no_defs")))

    # plan(multisession, workers = 3)

    extracted <- gs_tmax %>%
      inner_join(weights, by = join_by(Item)) %>%
      mutate(tmax = pmap(list(tmax, weights), function(abrick, aweight) {
        exact_extract(abrick, region, weights = aweight, fun = "weighted_mean", default_weight = 0, stack_apply = T, progress = F) %>%
          as_tibble() %>%
          bind_cols(region %>% st_drop_geometry()) %>%
          pivot_longer(-c(ADM_CODE, ADM_NAME, FAO_CODE),
            names_to = c("tmax", "year"), names_sep = "\\.",
            names_transform = list(year = as.integer),
            names_prefix = "weighted_mean.", values_drop_na = T
          ) %>%
          pivot_wider(names_from = tmax)
      }, .progress = F)) %>%
      select(-weights)

    return(extracted)
  }, .progress = T))

# qsave(data, 'data/outputs/tmax/gs_tmax_future.qs')
saveRDS(data, "data/outputs/tmax/extracted_tmax_future.rds")
