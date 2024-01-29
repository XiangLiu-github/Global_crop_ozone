source("script/loadPackages.R")

library(dtplyr)
plan(multisession, workers = 2)
US <- tibble(
  year = 2000:2019,
  data = future_map(year, function(ayear) {
    a = fread(str_c("D:\\Data\\Site_air_quality\\US\\hourly_44201_", ayear, ".csv"), 
              select = c('Latitude', 'Longitude', 'Date Local', 'Time Local', 'Sample Measurement'))
    
    a = a %>% 
      filter(`Sample Measurement` > 0) %>% 
      mutate(date = ymd_hm(stringi::stri_join(`Date Local`, ' ', `Time Local`)),
             site = stringi::stri_join(Latitude, '_', Longitude),
             o3 = `Sample Measurement` * 1e3, .keep = 'unused') %>% 
      openair::rollingMean(align = 'left') %>% 
      ungroup() %>% 
      mutate(year = year(date),
             month = month(date), 
             day = day(date), 
             hour = hour(date), .keep = 'unused') %>% 
      as.data.table()
    
    # MDA8
    MDA8 = a %>% 
      lazy_dt() %>% 
      drop_na() %>% 
      # daily
      group_by(year, month, day, site) %>% 
      filter(n() > 18) %>% 
      summarise(o3_max = max(rolling8o3), .groups = 'drop') %>% 
      # monthly
      group_by(year, month, site) %>% 
      filter(n() > days_in_month(ymd(stringi::stri_join(year, month, day, sep = '-'))) * 0.9) %>% 
      summarise(rollingo3 = mean(o3_max), .groups = 'drop') %>% 
      as_tibble()
    
    Mean = a %>% 
      lazy_dt() %>% 
      # daily
      group_by(year, month, day, site) %>% 
      filter(n() > 18) %>% 
      summarise(o3 = mean(o3), .groups = 'drop') %>% 
      # monthly
      group_by(year, month, site) %>% 
      filter(n() > days_in_month(ymd(stringi::stri_join(year, month, day, sep = '-'))) * 0.9) %>% 
      summarise(o3 = mean(o3), .groups = 'drop') %>% 
      as_tibble()
    
    Veg = a %>% 
      lazy_dt() %>% 
      filter(hour %between% c(8, 19)) %>% 
      mutate(AOT40 = fifelse(o3 >= 40, o3 - 40, 0),
             W126 = o3 / (1 + 4403 * exp(- 126 * o3 / 1000))) %>% 
      #daily
      group_by(site, year, month, day) %>% 
      summarise(AOT40 = sum(AOT40), W126 = sum(W126), data_capture = length(o3)/12, .groups = 'drop') %>%
      filter(data_capture >= 0.75) %>% 
      mutate(AOT40 = AOT40 / data_capture, W126 = W126 / data_capture) %>% 
      # mutate(across(c(AOT40, W126), ~ . / data_capture)) %>% 
      select(-data_capture) %>% 
      as.data.table() %>% 
      #monthly
      fgroup_by(site, year, month) %>% 
      fsummarise(AOT40 = sum(AOT40), W126 = sum(W126), count = length(AOT40)) %>% 
      fmutate(data_capture = count / days_in_month(ym(str_c(year, month, '-')))) %>% 
      filter(data_capture >= 0.75) %>% 
      mutate(across(c(AOT40, W126), ~ . / data_capture)) %>% 
      fselect(-data_capture, -count)
    
    list(MDA8, Mean, Veg) %>% 
      reduce(inner_join, by = join_by(year, month, site))
  }, .progress = T)
)
plan(sequential)

US = US %>% 
  select(-year) %>% 
  unnest() %>% 
  separate(site, c('Latitude', 'Longitude'), '_', convert = T)

# ppb and ppb h
US = US %>% 
  group_by(Latitude, Longitude) %>% 
  filter(n() > 5) %>% 
  ungroup()

saveRDS(US, "data/inputs/US_ozone.rds")



