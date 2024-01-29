source("script/loadpackages.R")
source("script/loadfunctions.R")

season <- get_season() %>%
  select(Item, var, prop) %>%
  mutate(prop = "0_10")

data <- get_data() %>%
  mutate(logYield = log(Yield))

data <- data %>%
  pivot_longer(c(aot40, starts_with("tmax"), starts_with("AOD"), starts_with("sm"), starts_with("cloud"))) %>%
  mutate(var = str_remove(name, "[23]")) %>%
  inner_join(season, by = join_by(Item, prop, var)) %>%
  select(-c(prop, var)) %>%
  pivot_wider() %>%
  drop_na() %>%
  nest(.by = Item, .key = "fdata")

n <- 1e3

plan(multisession, workers = 8)

cvdata <- data %>%
  mutate(
    results = future_map(fdata, function(afdata) {
      afdata <- feols(c(logYield, aot40, tmax, tmax2, tmax3, AOD, AOD2, AOD3, sm, sm2, sm3, cloud, cloud2, cloud3) ~ 1 | ADM_CODE[year + year^2], data = afdata) %>%
        map_dfc("residuals") %>%
        purrr::set_names(~ str_remove(.x, "lhs: ")) %>%
        mutate(ADM_CODE = afdata$ADM_CODE)

      set.seed(2023)

      bootstraps <- group_bootstraps(afdata, ADM_CODE, times = n)
      # bootstraps = mc_cv(afdata, times = n)

      map_dfr(bootstraps$splits, function(asplit) {
        analysisd <- analysis(asplit)
        assessmentd <- assessment(asplit)

        model <- feols(logYield ~ aot40 + tmax + tmax2 + tmax3 + AOD + AOD2 + AOD3 + sm + sm2 + sm3 + cloud + cloud2 + cloud3, data = analysisd)

        tibble(
          Training = cor(analysisd$logYield, predict(model, analysisd)),
          Testing = cor(assessmentd$logYield, predict(model, assessmentd))
        )
      }, .id = "id")
    }, .progress = T, .options = furrr_options(seed = T)), .keep = "unused"
  )

plan(sequential)

qsave(cvdata, "data/outputs/oob_whole.qs")

cvdata <- qread("data/outputs/oob_whole.qs")

summed <- cvdata %>%
  unnest() %>%
  group_by(Item) %>%
  median_qi(Testing, .width = uncertainty[2]) %>%
  group_by(.lower < 0) %>%
  arrange(desc(Testing), .by_group = T)

order <- summed %>%
  pull(Item)

alpha <- summed %>%
  filter(`.lower < 0` == F) %>%
  pull(Item)

plot_tbl <- tibble(
  x = rep(alpha[length(alpha)], 2),
  y = -0.25,
  label = c("Reliable crop models", "Unreliable crop models"),
  hjust = c(1, -0.3)
)

p1 <- cvdata %>%
  unnest() %>%
  pivot_longer(-c(Item, id)) %>%
  inner_join(color_tbl) %>%
  mutate(
    alpha = fifelse(Item %in% alpha, "1", "0"),
    Item = factor(Item, order)
  ) %>%
  ggplot(aes(x = Item, y = value)) +
  geom_hline(yintercept = 0, linetype = "longdash", color = "grey") +
  geom_vline(aes(xintercept = stage("any", after_stat = length(alpha) + 0.5)), color = "#F94A29", linewidth = 1) +
  stat_eye(aes(fill = factor(name, c("Training", "Testing")), alpha = alpha), position = position_dodge(0.8), .width = uncertainty) +
  geom_text(aes(x = Item, y = 0.6, label = Item, color = color), data = . %>% slice_head(n = 1, by = Item), size = 4.5, family = "Roboto Condensed") +
  geom_text(aes(x = x, y = y, label = label, hjust = hjust), data = plot_tbl, inherit.aes = F, size = 5, family = "Roboto Condensed") +
  scale_y_continuous(name = "Correlation coefficient") +
  scale_alpha_manual(values = c(0.4, 1), guide = "none") +
  scale_fill_manual(values = c("#008500", "#FF8500")) +
  scale_color_identity() +
  theme_half_open(font_size = 16, font_family = "Roboto Condensed") +
  background_grid(major = "x") +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.1, 0.2),
    axis.title.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

coefs <- qread("data/outputs/boot_coefs_whole.qs")

p2 <- inner_join(coefs, data) %>%
  mutate(aot40 = map2(coefs, fdata, function(acoef, afdata) {
    acoef <- acoef %>%
      select(starts_with("aot40"))

    X <- afdata %>%
      select(starts_with("aot40"))

    counter_X <- X %>%
      mutate(aot40 = aot40 - 1e3)

    rel_X <- counter_X - X

    tcrossprod(as.matrix(rel_X), as.matrix(acoef)) %>%
      expm1() %>%
      colWeightedMeans(w = afdata$Area_harvested)
  }), .keep = "unused") %>%
  unnest() %>%
  inner_join(color_tbl) %>%
  inner_join(data %>% unnest() %>% summarise(aot40_mean = mean(aot40) / 1e3, .by = Item)) %>%
  mutate(
    alpha = fifelse(Item %in% alpha, "1", "0"),
    Item = factor(Item, order),
    # aot40 = aot40 * aot40_mean
  ) %>%
  ggplot(aes(x = Item, y = aot40, alpha = alpha))

counter_tmax <- read_rds("data/outputs/tmax/extracted_tmax_count2_whole.rds") %>%
  select(-prop)

p3 <- list(coefs, data, counter_tmax) %>%
  reduce(inner_join) %>%
  mutate(tmax = pmap(list(coefs, fdata, tmax), function(acoef, afdata, atdata) {
    atdata <- atdata %>%
      rename_with(~ str_c("counter_", .x), starts_with("tmax"))

    afdata <- inner_join(afdata, atdata, by = join_by(FAO_CODE, ADM_NAME, ADM_CODE, year))

    acoef <- acoef %>%
      select(starts_with("tmax"))

    X <- afdata %>%
      select(starts_with("tmax"))

    counter_X <- afdata %>%
      select(starts_with("counter_tmax"))

    rel_X <- counter_X - X

    tcrossprod(as.matrix(rel_X), as.matrix(acoef)) %>%
      expm1() %>%
      colWeightedMeans(w = afdata$Area_harvested)
  }), .keep = "unused") %>%
  unnest() %>%
  inner_join(color_tbl) %>%
  mutate(
    alpha = fifelse(Item %in% alpha, "1", "0"),
    Item = factor(Item, order)
  ) %>%
  ggplot(aes(x = Item, y = tmax, alpha = alpha))

make_plot <- function(ap) {
  ap +
    geom_hline(yintercept = 0, linetype = "longdash", color = "grey") +
    geom_vline(aes(xintercept = stage("any", after_stat = length(alpha) + 0.5)), color = "#F94A29", linewidth = 1) +
    stat_summary(aes(fill = color), color = "black", geom = "bar", fun = median, width = 0.5) +
    stat_pointinterval(.width = uncertainty) +
    scale_y_percent(name = "Change in yield") +
    scale_color_identity() +
    scale_fill_identity() +
    scale_alpha_manual(values = c(0.4, 1), guide = "none") +
    theme_half_open(font_size = 16, font_family = "Roboto Condensed") +
    background_grid(major = "x") +
    theme(
      legend.title = element_blank(),
      legend.position = c(0.1, 0.3),
      axis.title.x = element_blank(),
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      strip.background = element_blank(),
      strip.text.x = element_blank()
    )
}

p2 <- make_plot(p2)
p3 <- make_plot(p3)

patch <- wrap_plots(p1, p2, p3, heights = c(1, 0.5, 0.5)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 30))

ggsave("figures/oob_whole.pdf", patch, width = 4, height = 2.5, scale = 3.5)
