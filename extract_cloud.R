source("script/loadpackages.R")
source("script/loadfunctions.R")

# a = crossing(year = 2000:2017,
#              month = str_pad(1:12, 2, pad = '0')) %>%
#   mutate(url = str_c('https://www.ncei.noaa.gov/thredds/ncss/cdr/isccp_hgm_agg/files/isccp-basic/hgm/ISCCP-Basic.HGM.v01r00.GLOBAL.', year, '.', month, '.99.9999.GPC.10KM.CS00.EA1.00.nc?var=tau,cldamt'),
#          destfile = str_c('ISCCP-Basic.HGM.v01r00.GLOBAL.', year, '.', month, '.99.9999.GPC.10KM.CS00.EA1.00.nc'))
#
# pro_walk2(a$url, a$destfile, curl::curl_fetch_disk)

a <- nc_open("D:/Data/Desktop/CER_SSF1deg-Day_Aqua-MODIS_Edition4A_400405.200207.hdf")
# 4 dimensions:
#   day_of_month  Size:31
#   long_name: Day of the Month
#   units: N/A
#   latitude  Size:180
#   long_name: Latitude
#   units: degrees_north
#   longitude  Size:360
#   long_name: Longitude
#   units: degrees_east
#   cloud_layer  Size:5
#   long_name: Index of Cloud Layers Stratified by Pressure
#   units: N/A
#   comment: 1 = High (50-300 mb), 2 = UpperMid (300-500 mb),3 = LowerMid (500-700 mb), 4 = Low (700 mb-Surface),5 = Total (50 mb-Surface)

# check variables
map_dfc(a$var, "longname") %>%
  pivot_longer(everything()) %>%
  view()

ncdf4::ncvar_get(a, "cld_od") %>%
  aperm(c(2, 1, 3, 4)) %>%
  .[, , , 5] %>%
  rast() %>%
  focal(w = 9, fun = mean, na.policy = "only", na.rm = T) %>%
  plot(20)


# https://asdc.larc.nasa.gov/data/CERES/SSF1deg-Day/Aqua-MODIS_Edition4A/
# set http_proxy=http://172.28.208.1:1080
# set https_proxy=http://172.28.208.1:1080
#
# wget.exe --header "Authorization: Bearer eyJ0eXAiOiJKV1QiLCJvcmlnaW4iOiJFYXJ0aGRhdGEgTG9naW4iLCJzaWciOiJlZGxqd3RwdWJrZXlfb3BzIiwiYWxnIjoiUlMyNTYifQ.eyJ0eXBlIjoiVXNlciIsInVpZCI6InhpYW5nbGl1IiwiZXhwIjoxNjgyMTM5ODcyLCJpYXQiOjE2NzY5NTU4NzIsImlzcyI6IkVhcnRoZGF0YSBMb2dpbiJ9.f4m-WWIMA-AVc5SI1K7Pl12XGOe9SRUX52y5PEeUXWuUqOcWa12LUGJ-LjzBXywr0dNVNK2CWHLRwDQdOObQgZoCdrS3ae-lHRzGNktdfx3y-rLxOGnE2o47miaxezHZpK9zH64AV55BL9wSVPMnvNBa9IuMDxuVQ97QnptTx_7eNGhQlpEU_hKkiEmlYN4SSHWeAOS6TduxuRnyOvIQ_28I8_VuHh83ke4e7k2C3FSNXypMva5DUAEmoCaNHzXxpcgumsUTOVmA2BPU-wtjVJLuicxEcb-cJCqj-ZMOhiQ5jmcCiMubozlJkUQTKrvlUGmF-zRLTGYvynYQgQGjjQ" --recursive --no-parent --reject "index.html*" --execute robots=off --proxy=on https://asdc.larc.nasa.gov/data/CERES/SSF1deg-Day/Aqua-MODIS_Edition4A/2003/02/CER_SSF1deg-Day_Aqua-MODIS_Edition4A_400405.200302.hdf
# ~ 200Gb ~ 13h

temp <- raster("data/inputs/crop_calendar/mai_rf_ggcmi_crop_calendar_phase3_v1.01.nc4")

a <- stack(
  stack("data/inputs/cloud/CERES_SYN1deg-Day_Terra-Aqua-MODIS_Ed4.1_Subset_20030101-20130112.nc", varname = "cldtau_total_daily"),
  stack("data/inputs/cloud/CERES_SYN1deg-Day_Terra-Aqua-MODIS_Ed4.1_Subset_20130113-20191231.nc", varname = "cldtau_total_daily")
)

b <- stack(
  stack("data/inputs/cloud/CERES_SYN1deg-Day_Terra-Aqua-MODIS_Ed4.1_Subset_20030101-20130112.nc", varname = "cldarea_total_daily"),
  stack("data/inputs/cloud/CERES_SYN1deg-Day_Terra-Aqua-MODIS_Ed4.1_Subset_20130113-20191231.nc", varname = "cldarea_total_daily")
)

plan(multisession, workers = 5)

data <- tibble(
  name = names(a),
  date = str_remove(name, "X") %>% ymd(),
  year = year(date), month = month(date), day = day(date)
) %>%
  nest(.by = year) %>%
  mutate(data = map(data, function(adata) {
    adata %>%
      mutate(data = future_map(name, function(aname) {
        rotate(a[[aname]] * b[[aname]] / 100) %>%
          resample(temp, method = "ngb") %>%
          mask(temp) %>%
          `names<-`("cloud") %>%
          as.data.frame(xy = T, na.rm = T)
      }, .progress = F)) %>%
      unnest() %>%
      fgroup_by(x, y, month) %>%
      fsummarise(cloud = fmean(cloud), cloud2 = fmean(cloud^2), cloud3 = fmean(cloud^3))
  }, .progress = T)) %>%
  unnest()

plan(sequential)

qsave(data, "data/outputs/cloud/cloud.qs")

library(qs)
library(raster)
library(tidyverse)
library(exactextractr)

data <- qread("data/outputs/cloud/cloud.qs")

gs <- qread("data/outputs/gs/growing_season.qs", nthreads = 5)

library(dtplyr)
gs_cloud <- gs %>%
  mutate(cloud = map(growing_season, function(ags) {
    gc()

    ags %>%
      mutate(fdata = map(fdata, function(afdata) {
        left_join(afdata, data, by = join_by(x, y, year, month)) %>%
          lazy_dt() %>%
          group_by(Index_year, x, y) %>%
          filter(sum(is.na(cloud)) == 0) %>%
          summarise(across(starts_with("cloud"), mean),
            .groups = "drop"
          ) %>%
          as_tibble()
      })) %>%
      unnest(fdata) %>%
      pivot_wider(
        names_from = c(Index_year, prop),
        names_sep = ".",
        values_from = starts_with("cloud")
      )

    # left_join(ags, data, by = join_by(x, y, year, month)) %>%
    #   lazy_dt() %>%
    #   group_by(Index_year, prop, x, y) %>%
    #   filter(sum(is.na(cloud)) == 0) %>%
    #   summarise(across(starts_with('cloud'), mean),
    #             .groups = 'drop') %>%
    #   pivot_wider(names_from = c(Index_year, prop),
    #               names_sep = ".",
    #               values_from = starts_with('cloud')) %>%
    #   as_tibble()
  }, .progress = T)) %>%
  select(-growing_season)

gs_cloud <- bind_rows(
  gs_cloud %>%
    filter(!str_detect(Item, "Rice")),
  gs_cloud %>%
    filter(str_detect(Item, "Rice")) %>%
    mutate(Item = "Rice") %>%
    unnest() %>%
    lazy_dt() %>%
    summarise(across(everything(), mean), .by = c(Item, x, y)) %>%
    as_tibble() %>%
    nest(.by = Item, .key = "cloud")
) %>%
  mutate(cloud = map(cloud, ~ rasterFromXYZ(.x, crs = "+proj=longlat +datum=WGS84 +no_defs")))

qsave(gs_cloud, "data/outputs/cloud/gs_cloud.qs", nthreads = 5)

gs_cloud <- qread("data/outputs/cloud/gs_cloud.qs", nthreads = 5)

# plan(multisession, workers = 3)

region <- read_rds("data/inputs/shp/GAUL.rds")

weights <- qread("data/inputs/cropland/weights.qs")

extracted <- gs_cloud %>%
  inner_join(weights) %>%
  mutate(cloud = pmap(list(cloud, weights), function(abrick, aweight) {
    exact_extract(abrick, region, weights = aweight, fun = "weighted_mean", default_weight = 0, stack_apply = T, progress = F) %>%
      as_tibble() %>%
      bind_cols(region %>% sf::st_drop_geometry()) %>%
      pivot_longer(-c(ADM_CODE, ADM_NAME, FAO_CODE),
        names_to = c("cloud", "year", "prop"), names_sep = "\\.",
        names_transform = list(year = as.integer),
        names_prefix = "weighted_mean.", values_drop_na = T
      ) %>%
      pivot_wider(names_from = cloud)
  }, .progress = T)) %>%
  select(-weights)

saveRDS(extracted, "data/outputs/cloud/extracted_cloud.rds")


# interact irrigation -----------------------------------------------------

order <- get_season() %>%
  filter(var == "cloud") %>%
  select(Item, prop)

gs_cloud <- qread("data/outputs/cloud/gs_cloud.qs", nthreads = 5)

gs_cloud <- inner_join(order, gs_cloud) %>%
  mutate(cloud = map2(prop, cloud, function(aprop, acloud) {
    acloud %>%
      subset(names(acloud)[str_detect(names(acloud), aprop)])
  })) %>%
  select(-prop)

# plan(multisession, workers = 3)

region <- read_rds("data/inputs/shp/GAUL.rds")

weights <- qread("data/inputs/cropland/weights.qs")

irr <- qread("data/inputs/GMIA/irr.qs")

extracted <- gs_cloud %>%
  inner_join(weights) %>%
  inner_join(irr) %>%
  mutate(cloud = pmap(list(cloud, weights, irr_frc), function(abrick, aweight, airr_frc) {
    abrick <- (abrick * airr_frc) %>%
      `names<-`(names(abrick))

    exact_extract(abrick, region, weights = aweight, fun = "weighted_mean", default_weight = 0, stack_apply = T, progress = F) %>%
      as_tibble() %>%
      bind_cols(region %>% sf::st_drop_geometry()) %>%
      pivot_longer(-c(ADM_CODE, ADM_NAME, FAO_CODE),
        names_to = c("cloud", "year", "prop"), names_sep = "\\.",
        names_transform = list(year = as.integer),
        names_prefix = "weighted_mean.", values_drop_na = T
      ) %>%
      pivot_wider(names_from = cloud)
  }, .progress = T)) %>%
  select(-c(weights, irr_frc))

saveRDS(extracted, "data/outputs/cloud/extracted_cloud_itr_irr.rds")
