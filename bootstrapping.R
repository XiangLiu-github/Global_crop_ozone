source("script/loadpackages.R")
source("script/loadfunctions.R")

season <- get_season() %>%
  select(Item, var, prop)
order <- get_order()

data <- get_data() %>%
  pivot_longer(c(aot40, starts_with("tmax"), starts_with("AOD"), starts_with("sm"), starts_with("cloud"))) %>%
  mutate(var = str_remove(name, "[23]")) %>%
  inner_join(season, by = join_by(Item, prop, var)) %>%
  select(-c(prop, var)) %>%
  pivot_wider() %>%
  drop_na() %>%
  nest(.by = Item, .key = "fdata")

tmax_rmax <- data %>%
  unnest() %>%
  summarise(across(c(tmax), mean), .by = c(Item, FAO_CODE, ADM_CODE, ADM_NAME)) %>%
  slice_min(tmax, prop = 1 - 0.1, by = Item) %>%
  select(-tmax) %>%
  inner_join(data %>% unnest()) %>%
  nest(.by = Item, .key = "fdata")

aot40_rmax <- data %>%
  unnest() %>%
  summarise(across(c(aot40), mean), .by = c(Item, FAO_CODE, ADM_CODE, ADM_NAME)) %>%
  slice_min(aot40, prop = 1 - 0.1, by = Item) %>%
  select(-aot40) %>%
  inner_join(data %>% unnest()) %>%
  nest(.by = Item, .key = "fdata")

tmax_rmin <- data %>%
  unnest() %>%
  summarise(across(c(tmax), mean), .by = c(Item, FAO_CODE, ADM_CODE, ADM_NAME)) %>%
  slice_max(tmax, prop = 1 - 0.1, by = Item) %>%
  select(-tmax) %>%
  inner_join(data %>% unnest()) %>%
  nest(.by = Item, .key = "fdata")

aot40_rmin <- data %>%
  unnest() %>%
  summarise(across(c(aot40), mean), .by = c(Item, FAO_CODE, ADM_CODE, ADM_NAME)) %>%
  slice_max(aot40, prop = 1 - 0.1, by = Item) %>%
  select(-aot40) %>%
  inner_join(data %>% unnest()) %>%
  nest(.by = Item, .key = "fdata")

n <- 1e3

plan(multisession, workers = 8)

res <- data %>%
  mutate(coefs = future_map(fdata, function(afdata) {
    set.seed(2023)

    bootstraps <- group_bootstraps(afdata, ADM_NAME, times = n)

    map_dfr(bootstraps$splits, function(asplit) {
      aadata <- analysis(asplit)

      feols(log(Yield) ~ aot40 + tmax + tmax2 + tmax3 + AOD + AOD2 + AOD3 + sm + sm2 + sm3 + cloud + cloud2 + cloud3 | ADM_CODE[year + year^2], aadata, lean = T) %>%
        coef()
    }, .id = "id")
  }, .progress = T, .options = furrr_options(seed = T)), .keep = "unused")

qsave(res, "data/outputs/boot_coefs.qs")

res <- data %>%
  mutate(coefs = future_map(fdata, function(afdata) {
    set.seed(2023)

    bootstraps <- group_bootstraps(afdata, ADM_NAME, times = n)

    map_dfr(bootstraps$splits, function(asplit) {
      aadata <- analysis(asplit)

      feols(log(Yield) ~ tmax + tmax2 + tmax3 + AOD + AOD2 + AOD3 + sm + sm2 + sm3 + cloud + cloud2 + cloud3 | ADM_CODE[year + year^2], aadata, lean = T) %>%
        coef()
    }, .id = "id")
  }, .progress = T, .options = furrr_options(seed = T)), .keep = "unused")

qsave(res, "data/outputs/boot_coefs_noaot40.qs")

res <- data %>%
  mutate(coefs = future_map(fdata, function(afdata) {
    set.seed(2023)

    bootstraps <- group_bootstraps(afdata, ADM_NAME, times = n)

    map_dfr(bootstraps$splits, function(asplit) {
      aadata <- analysis(asplit)

      feols(log(Yield) ~ aot40 + AOD + AOD2 + AOD3 + sm + sm2 + sm3 + cloud + cloud2 + cloud3 | ADM_CODE[year + year^2], aadata, lean = T) %>%
        coef()
    }, .id = "id")
  }, .progress = T, .options = furrr_options(seed = T)), .keep = "unused")

qsave(res, "data/outputs/boot_coefs_notmax.qs")

res <- data %>%
  mutate(coefs = future_map(fdata, function(afdata) {
    set.seed(2023)

    bootstraps <- group_bootstraps(afdata, ADM_NAME, times = n)

    map_dfr(bootstraps$splits, function(asplit) {
      aadata <- analysis(asplit)

      feols(log(Yield) ~ aot40 + tmax + tmax2 + tmax3 + AOD + AOD2 + AOD3 + sm + sm2 + sm3 + cloud + cloud2 + cloud3 | ADM_CODE[year], aadata, lean = T) %>%
        coef()
    }, .id = "id")
  }, .progress = T, .options = furrr_options(seed = T)), .keep = "unused")

qsave(res, "data/outputs/boot_coefs_linearFE.qs")

res <- data %>%
  mutate(coefs = future_map(fdata, function(afdata) {
    set.seed(2023)

    bootstraps <- group_bootstraps(afdata, ADM_NAME, times = n)

    map_dfr(bootstraps$splits, function(asplit) {
      aadata <- analysis(asplit)

      feols(log(Yield) ~ aot40 + tmax + tmax2 + tmax3 + AOD + AOD2 + AOD3 + sm + sm2 + sm3 + cloud + cloud2 + cloud3 | ADM_CODE[year], aadata, lean = T, weights = ~ log(Area_harvested), notes = F) %>%
        coef()
    }, .id = "id")
  }, .progress = T, .options = furrr_options(seed = T)), .keep = "unused")

qsave(res, "data/outputs/boot_coefs_area_weighted.qs")

tmax_rmax <- tmax_rmax %>%
  mutate(coefs = future_map(fdata, function(afdata) {
    set.seed(2023)

    bootstraps <- group_bootstraps(afdata, ADM_NAME, times = n)

    map_dfr(bootstraps$splits, function(asplit) {
      aadata <- analysis(asplit)

      feols(log(Yield) ~ aot40 + tmax + tmax2 + tmax3 + AOD + AOD2 + AOD3 + sm + sm2 + sm3 + cloud + cloud2 + cloud3 | ADM_CODE[year + year^2], aadata, lean = T) %>%
        coef()
    }, .id = "id")
  }, .progress = T, .options = furrr_options(seed = T)), .keep = "unused")

qsave(tmax_rmax, "data/outputs/boot_coefs_tmax_0.1max_removed.qs")

aot40_rmax <- aot40_rmax %>%
  mutate(coefs = future_map(fdata, function(afdata) {
    set.seed(2023)

    bootstraps <- group_bootstraps(afdata, ADM_NAME, times = n)

    map_dfr(bootstraps$splits, function(asplit) {
      aadata <- analysis(asplit)

      feols(log(Yield) ~ aot40 + tmax + tmax2 + tmax3 + AOD + AOD2 + AOD3 + sm + sm2 + sm3 + cloud + cloud2 + cloud3 | ADM_CODE[year + year^2], aadata, lean = T) %>%
        coef()
    }, .id = "id")
  }, .progress = T, .options = furrr_options(seed = T)), .keep = "unused")

qsave(aot40_rmax, "data/outputs/boot_coefs_aot40_0.1max_removed.qs")

tmax_rmin <- tmax_rmin %>%
  mutate(coefs = future_map(fdata, function(afdata) {
    set.seed(2023)

    bootstraps <- group_bootstraps(afdata, ADM_NAME, times = n)

    map_dfr(bootstraps$splits, function(asplit) {
      aadata <- analysis(asplit)

      feols(log(Yield) ~ aot40 + tmax + tmax2 + tmax3 + AOD + AOD2 + AOD3 + sm + sm2 + sm3 + cloud + cloud2 + cloud3 | ADM_CODE[year + year^2], aadata, lean = T) %>%
        coef()
    }, .id = "id")
  }, .progress = T, .options = furrr_options(seed = T)), .keep = "unused")

qsave(tmax_rmin, "data/outputs/boot_coefs_tmax_0.1min_removed.qs")

aot40_rmin <- aot40_rmin %>%
  mutate(coefs = future_map(fdata, function(afdata) {
    set.seed(2023)

    bootstraps <- group_bootstraps(afdata, ADM_NAME, times = n)

    map_dfr(bootstraps$splits, function(asplit) {
      aadata <- analysis(asplit)

      feols(log(Yield) ~ aot40 + tmax + tmax2 + tmax3 + AOD + AOD2 + AOD3 + sm + sm2 + sm3 + cloud + cloud2 + cloud3 | ADM_CODE[year + year^2], aadata, lean = T) %>%
        coef()
    }, .id = "id")
  }, .progress = T, .options = furrr_options(seed = T)), .keep = "unused")

qsave(aot40_rmin, "data/outputs/boot_coefs_aot40_0.1min_removed.qs")

plan(sequential)
