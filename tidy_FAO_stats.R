source("script/loadpackages.R")
source("script/loadfunctions.R")

tidy_FAO <- function(aFAO_file) {
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

  # remove Sudan (former) duplicated
  table <- table %>% filter(FAO_CODE != 206)

  # FAO code tow china
  stats <- read_csv(aFAO_file) %>%
    mutate(across(everything(), as.character))

  # remove duplicates
  temp <- inner_join(stats, table, c("Area Code" = "FAO_CODE")) %>%
    distinct(`Area Code`, Area, ADM_NAME, ADM_CODE)

  dup <- temp$`Area Code`[duplicated(temp$`Area Code`)] %>% unique()

  temp <- bind_rows(
    temp %>%
      filter(
        `Area Code` %in% dup,
        ADM_NAME %in% c("Australia", "France", "Portugal", "Réunion", "U.K. of Great Britain and Northern", "Venezuela")
      ),
    temp %>%
      filter(!`Area Code` %in% dup)
  )

  stats <- inner_join(stats, temp) %>%
    select(`Area Code`, Item, Element, num_range("Y", 1961:2020), num_range("Y", 1961:2021, "F"), ADM_NAME, ADM_CODE) %>%
    rename(FAO_CODE = `Area Code`)

  return(stats)
}

# yield -------------------------------------------------------------------

stats <- tidy_FAO("data/inputs/FAO_stats/Production_Crops_Livestock_E_All_Data.csv")

# # A tibble: 3 × 3
# Element        Unit       n
# <chr>          <chr>  <int>
# 1 Area harvested ha      2006
# 2 Production     tonnes  2007
# 3 Yield          hg/ha   1913

datas <- stats %>%
  filter(Item %in% c(
    "Barley", "Cassava, fresh", "Groundnuts, excluding shelled", "Beans, dry",
    "Maize (corn)", "Millet", "Potatoes", "Rape or colza seed", "Rye", "Peas, dry",
    "Seed cotton, unginned", "Sorghum", "Soya beans", "Sugar beet",
    "Sugar cane", "Sunflower seed", "Wheat", "Rice"
  ) &
    Element %in% c("Yield", "Production", "Area harvested")) %>%
  pivot_longer(starts_with("Y"), names_prefix = "Y") %>%
  separate(name, c("year", "temp"), "F") %>%
  mutate(
    year = as.numeric(year),
    temp = fifelse(is.na(temp), "raw", "flag")
  ) %>%
  pivot_wider(names_from = temp) %>%
  drop_na(raw) %>%
  rename(value = raw)

datas <- datas %>%
  filter(flag %in% c("A")) %>%
  # filter(flag %in% c('A', 'E')) %>%
  select(-flag) %>%
  mutate(value = as.numeric(value)) %>%
  pivot_wider(names_from = Element, values_from = value) %>%
  rename(Area_harvested = `Area harvested`) %>%
  filter(year >= 2003) %>%
  mutate(Item = case_when(
    Item == "Maize (corn)" ~ "Maize",
    Item == "Cassava, fresh" ~ "Cassava",
    Item == "Beans, dry" ~ "Bean",
    Item == "Peas, dry" ~ "Pea",
    Item == "Seed cotton, unginned" ~ "Cotton",
    Item == "Groundnuts, excluding shelled" ~ "Groundnut",
    Item == "Potatoes" ~ "Potato",
    Item == "Soya beans" ~ "Soybean",
    Item == "Sugar beet" ~ "Sugarbeet",
    Item == "Sugar cane" ~ "Sugarcane",
    Item == "Sunflower seed" ~ "Sunflower",
    Item == "Rape or colza seed" ~ "Rapeseed",
    TRUE ~ Item
  ))

datas %>%
  group_by(Item) %>%
  tally()

saveRDS(datas, "data/outputs/yields.rds")

# irrigation --------------------------------------------------------------

stats <- tidy_FAO("data/inputs/FAO_stats/Inputs_LandUse_E_All_Data_NOFLAG.csv")

stats <- stats %>%
  filter(Item %in% c("Agricultural land", "Land area equipped for irrigation")) %>%
  pivot_longer(starts_with("Y"), names_to = "year", names_prefix = "Y", names_transform = list(year = as.integer), values_transform = list(value = as.double)) %>%
  filter(year >= 2003) %>%
  pivot_wider(names_from = Item, values_from = value) %>%
  replace_NA() %>%
  set_names(~ str_replace_all(.x, " ", "_")) %>%
  mutate(irg_fraction = Land_area_equipped_for_irrigation / Agricultural_land, .keep = "unused")

saveRDS(stats, "data/outputs/irg.rds")


# food balance region-level -----------------------------------------------

indi <- read_csv("data/inputs/FAO_stats/FoodBalanceSheets_E_ItemCodes.csv") %>%
  slice(2:99)

fbdata <- read_csv("data/inputs/FAO_stats/FoodBalanceSheets_E_All_Area_Groups_NOFLAG.csv") %>%
  select(Area, Item, `Item Code`, Element, num_range("Y", 2010:2020)) %>%
  filter(Area %in% c(
    "Eastern Africa", "Middle Africa", "Northern Africa", "Southern Africa",
    "Western Africa", "Northern America", "Central America", "Caribbean",
    "South America", "Central Asia", "Eastern Asia", "Southern Asia",
    "South-Eastern Asia", "Western Asia", "Eastern Europe", "Northern Europe",
    "Southern Europe", "Western Europe", "Oceania"
  )) %>%
  pivot_longer(num_range("Y", 2010:2020),
    names_prefix = "Y", names_to = "Year", names_transform = list(Year = as.integer)
  )

data <- list(
  fbdata %>%
    filter(Item == "Grand Total", Element == "Food supply (kcal/capita/day)"),
  fbdata %>%
    filter(Item == "Population"),
  fbdata %>%
    filter(`Item Code` %in% indi$`Item Code`, Element == "Food supply (kcal/capita/day)")
) %>%
  map(~ select(.x, -`Item Code`)) %>%
  map(~ select(pivot_wider(.x, names_from = Item, values_fn = min, values_fill = 0), -c(Element))) %>%
  reduce(inner_join) %>%
  mutate(Food_supply = `Grand Total`, .keep = "unused") %>%
  janitor::clean_names() %>%
  replace_NA()

# check whether calculated sum equals reported sum
tibble(
  a = select(data, -c(area, year, food_supply, population)) %>% rowSums(),
  b = data$food_supply,
  area = data$area
) %>%
  ggplot(aes(x = a, y = b, color = area)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  stat_poly_eq(use_label(c("rr.label", "eq.label")))

# check people at hunger risk

FS <- fread("data/inputs/FAO_stats/Food_Security_Data_E_All_Data_NOFLAG.csv") %>%
  as_tibble() %>%
  filter(Area %in% c(
    "Eastern Africa", "Middle Africa", "Northern Africa", "Southern Africa",
    "Western Africa", "Northern America", "Central America", "Caribbean",
    "South America", "Central Asia", "Eastern Asia", "Southern Asia",
    "South-Eastern Asia", "Western Asia", "Eastern Europe", "Northern Europe",
    "Southern Europe", "Western Europe", "Oceania"
  ))

FS %>%
  filter(Item == "Number of severely food insecure people (million) (annual value)") %>%
  select(Area, Element, num_range("Y", 2014:2020)) %>%
  mutate(across(starts_with("Y"), as.double))

# CV only available at country-level
FS %>%
  filter(Item == "Coefficient of variation of habitual caloric consumption distribution (real number)") %>%
  view()

FS %>%
  filter(Item == "Minimum dietary energy requirement  (kcal/cap/day)") %>%
  view()

# food balance country-level ----------------------------------------------

countries <- get_FAO_groups()

indi <- read_csv("data/inputs/FAO_stats/FoodBalanceSheets_E_ItemCodes.csv") %>%
  slice(2:99)

fbdata1 <- read_csv("data/inputs/FAO_stats/FoodBalanceSheets_E_All_Data_NOFLAG.csv") %>%
  select(Area, Item, `Item Code`, Element, num_range("Y", 2010:2020)) %>%
  filter(Area %in% countries$ADM_NAME) %>%
  pivot_longer(num_range("Y", 2010:2020),
    names_prefix = "Y", names_to = "Year", names_transform = list(Year = as.integer)
  )

fbdata2 <- read_csv("data/inputs/FAO_stats/FoodBalanceSheetsHistoric_E_All_Data_NOFLAG.csv") %>%
  select(Area, Item, `Item Code`, Element, num_range("Y", 1999:2009)) %>%
  filter(Area %in% countries$ADM_NAME) %>%
  pivot_longer(num_range("Y", 1990:2009),
    names_prefix = "Y", names_to = "Year", names_transform = list(Year = as.integer)
  )

fbdata <- bind_rows(fbdata1, fbdata2)

data <- list(
  fbdata %>%
    filter(Item == "Grand Total", Element == "Food supply (kcal/capita/day)"),
  fbdata %>%
    filter(Element == "Total Population - Both sexes"),
  fbdata %>%
    filter(`Item Code` %in% indi$`Item Code`, Element == "Food supply (kcal/capita/day)")
) %>%
  map(~ select(.x, -`Item Code`)) %>%
  map(~ select(pivot_wider(.x, names_from = Item, values_fn = min, values_fill = 0), -c(Element))) %>%
  reduce(inner_join) %>%
  mutate(Food_supply = `Grand Total`, .keep = "unused") %>%
  replace_NA() %>%
  filter(Food_supply != 0)

# check whether calculated sum equals reported sum
# historical data is not balanced
tibble(
  a = select(data, -c(Area, Year, Food_supply, Population)) %>% rowSums(),
  b = data$Food_supply,
  area = data$Area,
  year = factor(data$Year)
) %>%
  ggplot(aes(x = a, y = b, color = year)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  stat_poly_eq(use_label(c("rr.label", "eq.label")))

# check people at hunger risk

FS <- fread("data/inputs/FAO_stats/Food_Security_Data_E_All_Data_NOFLAG.csv") %>%
  as_tibble() %>%
  filter(Area %in% countries$ADM_NAME)

# country level only has three-year average
reported <- FS %>%
  filter(Item == "Prevalence of undernourishment (percent) (3-year average)") %>%
  select(Area, Element, num_range("Y", str_c(2000:2019, 2002:2021))) %>%
  mutate(across(starts_with("Y"), as.double)) %>%
  pivot_longer(starts_with("Y"),
    names_prefix = "Y", names_to = "year", names_transform = list(year = as.integer),
    values_transform = list(value = as.double)
  ) %>%
  mutate(
    start = as.integer(str_sub(year, end = 4)),
    end = as.integer(str_sub(year, start = 5)),
    year = (start + end) / 2
  )

# CV only available at country-level
para <- FS %>%
  filter(Item %in% c("Coefficient of variation of habitual caloric consumption distribution (real number)", "Minimum dietary energy requirement  (kcal/cap/day)")) %>%
  select(Area, Item, num_range("Y", 2000:2021)) %>%
  pivot_longer(starts_with("Y"),
    names_prefix = "Y", names_to = "year", names_transform = list(year = as.integer),
    values_transform = list(value = as.double)
  ) %>%
  pivot_wider(names_from = Item) %>%
  rename(
    CV = `Coefficient of variation of habitual caloric consumption distribution (real number)`,
    rl = `Minimum dietary energy requirement  (kcal/cap/day)`
  )

caled <- inner_join(
  data %>%
    select(Area, Year, Food_supply),
  para, join_by(Area, Year == year)
) %>%
  mutate(
    sigma = sqrt(log(CV^2 + 1)),
    mu = log(Food_supply) - sigma^2 / 2,
    pn = pnorm((log(rl) - mu) / sigma) * 100
  ) %>%
  mutate(pn_3year = slide_dbl(pn, mean, .before = 1, .after = 1, .complete = T), .by = Area)

# for somewhat reason, the calculation is smaller than FAO reported
inner_join(caled, reported, join_by(Area, Year == year)) %>%
  view() %>%
  ggplot(aes(x = pn_3year, y = value)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  stat_poly_eq(use_label(c("rr.label", "eq.label"))) +
  scale_y_log10() +
  scale_x_log10()


cv <- 0.23
x <- 2964
rl <- 1845

sigma <- sqrt(log(cv^2 + 1))
mu <- log(x) - sigma^2 / 2
pnorm((log(rl) - mu) / sigma) * 100
