source("script/loadpackages.R")
source("script/loadfunctions.R")

season <- get_season() %>%
  select(Item, var, prop) %>%
  mutate(prop = "0_10")

data <- get_data() %>%
  pivot_longer(c(aot40, starts_with("tmax"), starts_with("AOD"), starts_with("sm"), starts_with("cloud"))) %>%
  mutate(var = str_remove(name, "[23]")) %>%
  inner_join(season, by = join_by(Item, prop, var)) %>%
  select(-c(prop, var)) %>%
  pivot_wider() %>%
  drop_na() %>%
  nest(.by = Item, .key = "fdata")

n <- 1e3

plan(multisession, workers = 8)

data <- data %>%
  mutate(coefs = future_map(fdata, function(afdata) {
    set.seed(2023)

    bootstraps <- group_bootstraps(afdata, ADM_NAME, times = n)

    map_dfr(bootstraps$splits, function(asplit) {
      aadata <- analysis(asplit)

      feols(log(Yield) ~ aot40 + tmax + tmax2 + tmax3 + AOD + AOD2 + AOD3 + sm + sm2 + sm3 + cloud + cloud2 + cloud3 | ADM_CODE[year + year^2], aadata, lean = T) %>%
        coef()
    }, .id = "id")
  }, .progress = T, .options = furrr_options(seed = T)), .keep = "unused")

qsave(data, "data/outputs/boot_coefs_whole.qs")

plan(sequential)
