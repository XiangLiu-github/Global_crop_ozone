source("script/loadpackages.R")
# source("script/loadfunctions.R")

# shp <- read_rds("data/inputs/shp/GAUL.rds") %>%
#   st_simplify(dTolerance = 0.1)
# 
# derice <- rast("data/inputs/crop_calendar/ri2_rf_ggcmi_crop_calendar_phase3_v1.01.nc4")["fraction_of_harvested_area"] %>%
#   as.data.frame(xy = T) %>%
#   filter(fraction_of_harvested_area >= 0.25) %>%
#   st_as_sf(coords = c("x", "y"), crs = 4326, remove = F) %>%
#   st_intersection(shp) %>%
#   st_drop_geometry() %>%
#   distinct(ADM_CODE, ADM_NAME, FAO_CODE) %>%
#   mutate(Item = "Rice")

yields <- read_rds("data/outputs/yields.rds") # %>% anti_join(derice)
irg <- read_rds("data/outputs/irg.rds")
cloud <- read_rds("data/outputs/cloud/extracted_cloud.rds") %>% unnest()
aot40 <- read_rds("data/outputs/ozone/extracted_aot40.rds") %>% unnest()
# tmean = read_rds('data/outputs/tmean/extracted_tmean.rds') %>% unnest()
# tmean_counter = read_rds('data/outputs/tmean/extracted_tmean_count2.rds') %>% unnest() %>%
# rename_with(~ str_c('counter_', .x), starts_with('tmean'))
tmax <- read_rds("data/outputs/tmax/extracted_tmax.rds") %>% unnest()
# tmax_counter = read_rds('data/outputs/tmax/extracted_tmax_count2.rds') %>% unnest() %>%
#   rename_with(~ str_c('counter_', .x), starts_with('tmax'))
aod <- read_rds("data/outputs/aod/extracted_aod.rds") %>% unnest()
# aod_counter = read_rds('data/outputs/aod/extracted_aod_count0.1.rds') %>% unnest() %>%
#   rename_with(~ str_c('counter_', .x), starts_with('AOD'))
sm <- read_rds("data/outputs/sm/extracted_sm.rds") %>% unnest()

a <- list(yields, irg, aot40, tmax, aod, sm, cloud) %>%
  reduce(left_join) %>%
  drop_na(Yield, prop, irg_fraction, aot40, tmax, AOD, sm, cloud) %>%
  mutate(
    PP = fifelse(Item %in% c("Maize", "Sugarcane", "Sorghum"), "C4", "C3"),
    continent = countrycode::countrycode(as.integer(ADM_CODE), "gaul", "continent", custom_match = c("151" = "Africa", "147296" = "Asia", "267" = "Asia", "147295" = "Asia", "2" = "Asia", "15" = "Asia", "40760" = "Africa", "136" = "Asia", "40762" = "Africa", "52" = "Asia", "91" = "Asia", "61013" = "Africa", "40781" = "Asia"))
  ) %>%
  filter(n() >= 10, .by = c(FAO_CODE, ADM_NAME, ADM_CODE, prop, Item))

a %>%
  group_by(Item, prop) %>%
  tally()

saveRDS(a, "data/tidied.rds")
