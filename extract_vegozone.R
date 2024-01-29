source("script/loadpackages.R")
source("script/loadfunctions.R")

###########################################################

ozone_situ = read_rds('data/inputs/US_ozone.rds') # ppb, and ppb h

model_aot40 <- bam(AOT40 ~ s(rollingo3), data = ozone_situ, family = tw(), discrete = T, nthreads = 9)

saveRDS(model_aot40, "data/outputs/model_aot40.rds")

model_aot40 <- read_rds("data/outputs/model_aot40.rds")

a <- brick("../pre_NS2021/data/results/emsb_o3_grouped_smoothed.grd") %>%
  as.data.frame(xy = T) %>%
  drop_na() %>%
  as_tibble() %>%
  pivot_longer(-c(x, y),
    names_to = c("year", "month"),
    names_transform = list(year = as.integer, month = as.double),
    names_prefix = "X",
    names_sep = "_",
    values_to = "rollingo3"
  ) %>%
  mutate(rollingo3 = rollingo3 * 1e3) %>% 
  mutate(aot40 = exp(predict(model_aot40, .))) %>%
  select(-rollingo3)

qsave(a, "data/outputs/ozone/aot40.qs")

library(qs)
library(raster)
library(tidyverse)
library(exactextractr)

data <- qread("data/outputs/ozone/aot40.qs")

gs <- qread("data/outputs/gs/growing_season.qs", nthreads = 5)

library(dtplyr)
gs_aot40 <- gs %>%
  mutate(aot40 = map(growing_season, function(ags) {
    gc()

    ags %>%
      mutate(fdata = map(fdata, function(afdata) {
        left_join(afdata, data, by = join_by(x, y, year, month)) %>%
          lazy_dt() %>%
          group_by(Index_year, x, y) %>%
          filter(sum(is.na(aot40)) == 0) %>%
          summarise(across(starts_with("aot40"), mean),
            .groups = "drop"
          ) %>%
          as_tibble()
      })) %>%
      unnest(fdata) %>%
      pivot_wider(
        names_from = c(Index_year, prop),
        names_sep = ".",
        values_from = starts_with("aot40")
      )

    # left_join(ags, data, by = join_by(x, y, year, month)) %>%
    #   lazy_dt() %>%
    #   group_by(Index_year, prop, x, y) %>%
    #   filter(sum(is.na(aot40)) == 0) %>%
    #   summarise(aot40 = sum(aot40),
    #             .groups = 'drop') %>%
    #   pivot_wider(names_from = c(Index_year, prop),
    #               names_sep = ".",
    #               values_from = c(aot40)) %>%
    #   as_tibble()
  }, .progress = T)) %>%
  select(-growing_season)

gc()

gs_aot40 <- bind_rows(
  gs_aot40 %>%
    filter(!str_detect(Item, "Rice")),
  gs_aot40 %>%
    filter(str_detect(Item, "Rice")) %>%
    mutate(Item = "Rice") %>%
    unnest() %>%
    lazy_dt() %>%
    summarise(across(everything(), mean), .by = c(Item, x, y)) %>%
    as_tibble() %>%
    nest(.by = Item, .key = "aot40")
) %>%
  mutate(aot40 = map(aot40, ~ rasterFromXYZ(.x, crs = "+proj=longlat +datum=WGS84 +no_defs")))

map_dbl(gs_aot40$aot40, inMemory)

qsave(gs_aot40, "data/outputs/ozone/gs_aot40.qs", nthreads = 5)

gs_aot40 <- qread("data/outputs/ozone/gs_aot40.qs", nthreads = 5)

# plan(multisession, workers = 3)

region <- read_rds("data/inputs/shp/GAUL.rds")

weights <- qread("data/inputs/cropland/weights.qs")

extracted <- gs_aot40 %>%
  inner_join(weights) %>%
  mutate(aot40 = pmap(list(aot40, weights), function(abrick, aweight) {
    exact_extract(abrick, region, weights = aweight, fun = "weighted_mean", default_weight = 0, stack_apply = T, progress = F) %>%
      as_tibble() %>%
      bind_cols(region %>% sf::st_drop_geometry()) %>%
      pivot_longer(-c(ADM_CODE, ADM_NAME, FAO_CODE),
        names_to = c("year", "prop"), names_transform = list(year = as.integer), names_prefix = "weighted_mean.X", names_sep = "\\.",
        values_to = "aot40", values_drop_na = T
      )
  }, .progress = T)) %>%
  select(-weights)

saveRDS(extracted, "data/outputs/ozone/extracted_aot40.rds")

# interact irrigation -----------------------------------------------------

order <- get_season() %>%
  filter(var == "aot40") %>%
  select(Item, prop)

gs_aot40 <- qread("data/outputs/ozone/gs_aot40.qs", nthreads = 5)

gs_aot40 <- inner_join(order, gs_aot40) %>%
  mutate(aot40 = map2(prop, aot40, function(aprop, aaot40) {
    aaot40 %>%
      subset(names(aaot40)[str_detect(names(aaot40), aprop)])
  })) %>%
  select(-prop)

# plan(multisession, workers = 3)

region <- read_rds("data/inputs/shp/GAUL.rds")

weights <- qread("data/inputs/cropland/weights.qs")

irr <- qread("data/inputs/GMIA/irr.qs")

extracted <- gs_aot40 %>%
  inner_join(weights) %>%
  inner_join(irr) %>%
  mutate(aot40 = pmap(list(aot40, weights, irr_frc), function(abrick, aweight, airr_frc) {
    abrick <- (abrick * airr_frc) %>%
      `names<-`(names(abrick))

    exact_extract(abrick, region, weights = aweight, fun = "weighted_mean", default_weight = 0, stack_apply = T, progress = F) %>%
      as_tibble() %>%
      bind_cols(region %>% sf::st_drop_geometry()) %>%
      pivot_longer(-c(ADM_CODE, ADM_NAME, FAO_CODE),
        names_to = c("year", "prop"), names_transform = list(year = as.integer), names_prefix = "weighted_mean.X", names_sep = "\\.",
        values_to = "aot40", values_drop_na = T
      )
  }, .progress = T)) %>%
  select(-c(weights, irr_frc))

saveRDS(extracted, "data/outputs/ozone/extracted_aot40_itr_irr.rds")

# future ----------------------------------------------------

ozone_situ = read_rds('data/inputs/US_ozone.rds') # ppb, and ppb h

model_aot40 <- read_rds("data/outputs/model_aot40.rds")

rate <- feols(log(rollingo3) ~ log(o3) | year + month + Longitude^Latitude, data = ozone_situ) %>%
  coef() %>%
  .[["log(o3)"]]

obs <- rast("../pre_NS2021/data/results/Satellite_MDA8.nc") %>% # ppb
  `time<-`(seq(ymd("2003-01-01"), ymd("2019-12-01"), "months")) %>%
  `names<-`(time(.))

obs <- obs[[year(time(obs)) %in% 2004:2010]]

data <- tibble(
  model = c("GFDL-ESM4_r1i1p1f1", "MPI-ESM1-2-HR_r1i1p1f1", "MRI-ESM2-0_r1i1p1f1", "NorESM2-LM_r1i1p1f1", "NorESM2-MM_r1i1p1f1", "UKESM1-0-LL_r1i1p1f2"),
  data = map(model, function(amodel) {
    fut_unadj <- rast(str_c("D:/data/CMIP6/remapped/o3_", amodel, "_ssp370.nc"), "o3")
    hist_unadj <- rast(str_c("D:/data/CMIP6/remapped/o3_", amodel, "_historical.nc"), "o3")
    hist_unadj <- c(hist_unadj, hist_unadj)

    arast <- (fut_unadj / hist_unadj) %>%
      crop(obs)

    fut <- ((arast ^ rate) * c(obs, obs)) %>%
      `names<-`(time(arast))

    fut %>%
      as.data.frame(xy = T) %>%
      pivot_longer(-c(x, y),
        names_to = c("year", "month", NA), names_transform = list(year = as.integer, month = as.integer), names_sep = "-",
        values_to = "rollingo3"
      ) %>%
      mutate(aot40 = exp(predict(model_aot40, ., cluster = cl))) %>%
      select(-rollingo3) %>%
      mutate(
        indx = fifelse(year <= 2050, "2045-2049", "2095-2099"),
        year = fifelse(year <= 2050, year - 40, year - 90)
      )
  }, .progress = T)
)

data <- data %>%
  unnest() %>%
  nest(.by = c(model, indx))

qsave(data, "data/outputs/ozone/aot40_future.qs")

data <- qread("data/outputs/ozone/aot40_future.qs")

order <- get_season() %>%
  filter(var == "aot40") %>%
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
    gs_aot40 <- gs %>%
      mutate(aot40 = map2(growing_season, prop, function(ags, aprop) {
        left_join(ags, adata, by = join_by(x, y, year, month)) %>%
          lazy_dt() %>%
          group_by(Index_year, x, y) %>%
          filter(sum(is.na(aot40)) == 0) %>%
          summarise(
            aot40 = sum(aot40),
            .groups = "drop"
          ) %>%
          pivot_wider(
            names_from = Index_year,
            names_sep = ".",
            values_from = c(aot40)
          ) %>%
          as_tibble()
      }, .progress = F)) %>%
      select(-growing_season)

    gs_aot40 <- bind_rows(
      gs_aot40 %>%
        filter(!str_detect(Item, "Rice")),
      gs_aot40 %>%
        filter(str_detect(Item, "Rice")) %>%
        mutate(Item = "Rice") %>%
        unnest() %>%
        lazy_dt() %>%
        summarise(across(everything(), mean), .by = c(Item, prop, x, y)) %>%
        as_tibble() %>%
        nest(.by = c(Item, prop), .key = "aot40")
    ) %>%
      mutate(aot40 = map(aot40, ~ rasterFromXYZ(.x, crs = "+proj=longlat +datum=WGS84 +no_defs")))

    # plan(multisession, workers = 3)

    extracted <- gs_aot40 %>%
      inner_join(weights, by = join_by(Item)) %>%
      mutate(aot40 = pmap(list(aot40, weights), function(abrick, aweight) {
        exact_extract(abrick, region, weights = aweight, fun = "weighted_mean", default_weight = 0, stack_apply = T, progress = F) %>%
          as_tibble() %>%
          bind_cols(region %>% sf::st_drop_geometry()) %>%
          pivot_longer(-c(ADM_CODE, ADM_NAME, FAO_CODE),
            names_to = "year", names_transform = list(year = as.integer), names_prefix = "weighted_mean.X",
            values_to = "aot40", values_drop_na = T
          )
      }, .progress = F)) %>%
      select(-weights)

    return(extracted)
  }, .progress = T))

saveRDS(data, "data/outputs/ozone/extracted_aot40_future.rds")
