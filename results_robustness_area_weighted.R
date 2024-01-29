source("script/loadpackages.R")
source("script/loadfunctions.R")

season <- get_season() %>%
  select(Item, var, prop)

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

coefs <- qread("data/outputs/boot_coefs.qs")
coefs_w <- qread("data/outputs/boot_coefs_area_weighted.qs")

p1 <- list(coefs, coefs_w) %>%
  map_dfr(function(acoefs) {
    inner_join(acoefs, data) %>%
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
      unnest()
  }, .id = "model") %>%
  group_by(model, Item) %>%
  median_qi(aot40, .width = uncertainty) %>%
  ggplot(aes(x = factor(Item, get_order()), y = aot40, ymax = .upper, ymin = .lower, color = factor(model, labels = c("Perfered specification", "Regression weighted by logarithm of harvested area"))))


counter_tmax <- read_rds("data/outputs/tmax/extracted_tmax_count2.rds") %>%
  select(-prop)

p2 <- list(coefs, coefs_w) %>%
  map_dfr(function(acoefs) {
    list(acoefs, data, counter_tmax) %>%
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
      unnest()
  }, .id = "model") %>%
  group_by(model, Item) %>%
  median_qi(tmax, .width = uncertainty) %>%
  ggplot(aes(x = factor(Item, get_order()), y = tmax, ymax = .upper, ymin = .lower, color = factor(model, labels = c("Perfered specification", "Regression weighted by logarithm of harvested area"))))

patch <- wrap_plots(p1, p2, ncol = 1, guides = "collect") &
  plot_annotation(tag_levels = "a") &
  geom_hline(yintercept = 0, linetype = "longdash", color = "grey30") &
  geom_pointinterval(position = position_dodge(0.4)) &
  scale_y_percent("Change in yield") &
  scale_color_manual(values = c("#EF4343", "#1399B2")) &
  theme_half_open(font_size = 16, font_family = "Roboto Condensed") &
  background_grid(major = "x") &
  theme(
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12),
    legend.position = "bottom",
    legend.justification = 0.5,
    plot.tag = element_text(size = 20)
  )

ggsave("figures/robustness_areaweighted.pdf", patch, width = 4, height = 2.5, scale = 3.5)
