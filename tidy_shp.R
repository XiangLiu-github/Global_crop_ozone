source("script/loadPackages.R")

# all one china
table <- read_xls("data/inputs/shp/G2015_InternationalCountryCodesAttributes.xls")

# all tidy china
table <- table %>%
  select(ADM_NAME, ADM_CODE, FAO_CODE) %>%
  separate_rows(ADM_CODE, sep = "\\+") %>%
  mutate(
    across(everything(), as.character),
    ADM_NAME = case_when(
      ADM_CODE == 147295 ~ "China, mainland",
      ADM_CODE == 147296 ~ "China, Taiwan Province of",
      TRUE ~ ADM_NAME
    ),
    FAO_CODE = case_when(
      ADM_CODE == 147295 ~ "41",
      ADM_CODE == 147296 ~ "214",
      TRUE ~ FAO_CODE
    )
  )

# ADM_CODE two china
gaul <- st_read("https://data.apps.fao.org/map/gsrv/gsrv1/gaul/wfs?service=WFS&version=1.1.0&request=GetFeature&typeName=gaul:g2015_2014_0&srsName=EPSG%3A4326&maxFeatures=1000", type = 6)

gaul <- gaul %>%
  mutate(across(-wkb_geometry, as.character)) %>%
  inner_join(table, c("adm0_code" = "ADM_CODE")) %>%
  select(adm0_code, ADM_NAME, FAO_CODE) %>%
  rename(ADM_CODE = adm0_code)

saveRDS(gaul, "data/inputs/shp/GAUL.rds")

gaul_sp <- gaul %>%
  st_simplify(dTolerance = 0.1)

gaul_sp %>%
  filter(adm0_code %in% a$ADM_CODE) %>%
  select(adm0_name, adm0_code) %>%
  plot()

plot(gaul_sp$wkb_geometry, add = T)



shp %>%
  st_drop_geometry() %>%
  view()


# ArcGIS Hub --------------------------------------------------------------

# https://hub.arcgis.com/datasets/esri::world-countries-generalized

shp <- st_read("data/inputs/shp/World_Countries_(Generalized).geojson") %>%
  filter(COUNTRY != "Antarctica") %>%
  summarise(geometry = st_union(geometry), .by = c(COUNTRYAFF, AFF_ISO)) %>%
  rename(COUNTRY = COUNTRYAFF, ISO2 = AFF_ISO) %>%
  mutate(
    COUNTRY = fifelse(COUNTRY == "Turkiye", "Turkey", COUNTRY),
    ISO3 = countrycode::countrycode(ISO2, "iso2c", "iso3c")
  ) %>%
  relocate(COUNTRY, ISO2, ISO3)

saveRDS(shp, "data/inputs/shp/tidied.rds")

plot(shp %>% select(-COUNTRY))

shp %>%
  st_drop_geometry() %>%
  select(COUNTRY, ISO2) %>%
  mutate(
    ISO2_ = countrycode::countrycode(COUNTRY, "country.name", "iso2c"),
    COUNTRY_ = countrycode::countrycode(ISO2, "iso2c", "country.name")
  ) %>%
  # filter(ISO != ISO_) %>%
  filter(COUNTRY != COUNTRY_) %>%
  view()

fao <- read_csv("data/inputs/FAO_stats/Production_Crops_Livestock_E_All_Data.csv", locale = locale(encoding = "ISO-8859-1")) %>%
  select(`Area Code`, `Area Code (M49)`, Area, Item, Element, num_range("Y", 1961:2021), num_range("Y", 1961:2021, "F")) %>%
  rename(
    FAO_CODE = `Area Code`,
    FAO_NAME = Area,
    M49 = `Area Code (M49)`
  ) %>%
  mutate(M49 = str_remove(M49, "'"))

faoinfo <- read_csv("data/inputs/shp/FAO_country_info.csv") %>%
  select(`Country Code`, Country, `M49 Code`, `ISO2 Code`, `ISO3 Code`) %>%
  rename(
    FAO_CODE = `Country Code`,
    FAO_NAME = Country,
    M49 = `M49 Code`,
    ISO2 = `ISO2 Code`,
    ISO3 = `ISO3 Code`
  )


anti_join(
  fao %>% distinct(FAO_CODE, FAO_NAME, M49),
  faoinfo
) # works

anti_join(
  faoinfo,
  fao %>% distinct(FAO_CODE, FAO_NAME, M49)
) %>% view()
