source("script/loadpackages.R")
source("script/loadfunctions.R")

season <- get_season() %>%
  select(Item, var, prop)
order <- get_order()

# pooled ------------------------------------------------------------

data <- get_data() %>%
  mutate(logYield = log(Yield))

data <- data %>%
  pivot_longer(c(aot40, starts_with("tmax"), starts_with("AOD"), starts_with("sm"), starts_with("cloud"))) %>%
  mutate(var = str_remove(name, "[23]")) %>%
  inner_join(season, by = join_by(Item, prop, var)) %>%
  select(-c(prop, var)) %>%
  pivot_wider() %>%
  drop_na()

saved <- qread("data/outputs/oob.qs") %>%
  unnest() %>%
  group_by(Item) %>%
  median_qi(Testing, .width = uncertainty[2]) %>%
  group_by(.lower < 0) %>%
  arrange(desc(Testing), .by_group = T) %>%
  filter(`.lower < 0` == F) %>%
  pull(Item)

data <- data %>%
  filter(Item %in% saved)

models <- list(
  fml_full = logYield ~ aot40 + tmax + tmax2 + tmax3 + AOD + AOD2 + AOD3 + sm + sm2 + sm3 + cloud + cloud2 + cloud3 | ADM_CODE^Item[year + year^2],
  fml_no_ozone = logYield ~ tmax + tmax2 + tmax3 + AOD + AOD2 + AOD3 + sm + sm2 + sm3 + cloud + cloud2 + cloud3 | ADM_CODE^Item[year + year^2],
  fml_no_tmax = logYield ~ aot40 + AOD + AOD2 + AOD3 + sm + sm2 + sm3 + cloud + cloud2 + cloud3 | ADM_CODE^Item[year + year^2]
) %>%
  map(feols, data = data, vcov = ~ ADM_CODE)

# How to determine whether the estimated coeficiences are significantly different between two models? For example, a model is `lm(yield ~ tmax + I(tmax ^ 2) + precipitation + ozone, data = data)` and another is `lm(yield ~ tmax + I(tmax ^ 2) + precipitation, data = data)`, I want to know whether adding the ozone variable will significantly alter the coefficients of temperature (tmax) on yield. Provide the R code.

ozone_seq <- seq(min(data$aot40), 12e3, by = 100)
tmax_seq <- seq(min(data$tmax), 45, by = 0.1)

pooled <- bind_rows(
  map_dfr(models, relpred,
    newdata = tibble(aot40 = ozone_seq),
    baseline = tibble(aot40 = rep(median(ozone_seq), length(ozone_seq))),
    x = ozone_seq, .id = "model"
  ) %>%
    filter(model != "fml_no_ozone") %>%
    mutate(var = "AOT40 (ppb h)"),
  map_dfr(models, relpred,
    newdata = tibble(tmax = tmax_seq, tmax2 = tmax^2, tmax3 = tmax^3),
    baseline = tibble(tmax = rep(median(tmax_seq), length(tmax_seq)), tmax2 = tmax^2, tmax3 = tmax^3),
    x = tmax_seq, .id = "model"
  ) %>%
    filter(model != "fml_no_tmax") %>%
    mutate(var = "Temperature (°C)")
) %>%
  mutate(across(c(fit, lwr, upr), expm1),
    model = factor(model, labels = c("Full model", "Full model without ozone", "Full model without temperature"))
  )

labels <- pooled %>%
  filter(x %in% c(ozone_seq[length(ozone_seq)], tmax_seq[length(tmax_seq)])) %>%
  select(model, x, fit, var) %>%
  pivot_wider(names_from = model, values_from = fit) %>%
  pivot_longer(-c(x, var, `Full model`), values_drop_na = T)

arrow <- arrow(
  length = unit(0.015, "npc"),
  ends = "last",
  type = "open"
)

p1 <- pooled %>%
  ggplot(aes(x = x, y = fit)) +
  # facet_wrap(~ TeX(var, output = "character"), labeller = label_parsed, scales = 'free', strip.position = 'bottom') + # slow
  facet_wrap(~var, scales = "free", strip.position = "bottom") +
  geom_hline(yintercept = 0, linetype = "longdash", color = "grey") +
  geom_ribbon(aes(ymax = upr, ymin = lwr, fill = model), alpha = 0.2) +
  geom_line(aes(color = model), show.legend = F, linewidth = 1.5) +
  geom_xsidehistogram(aes(x = value), data = data %>% select(tmax, aot40) %>% rename("Temperature (°C)" = tmax, "AOT40 (ppb h)" = aot40) %>% pivot_longer(everything(), names_to = "var"), fill = "grey", color = "black", inherit.aes = F, bins = 50) +
  geom_segment(aes(x = x * 1.02, xend = x * 1.02, y = value, yend = `Full model`), arrow = arrow, inherit.aes = F, data = labels) +
  geom_text_repel(aes(x = x * 1.18, y = (value + `Full model`) / 2, label = label), data = labels %>% mutate(label = str_c(round_pettry(-(value - `Full model`) * 100), "%")), inherit.aes = F, min.segment.length = Inf, box.padding = 0, family = "Roboto Condensed", size = 4) +
  geom_text(aes(x = x, y = y, label = str_c("P value = ", round_pettry(pv))),
    inherit.aes = F, family = "Roboto Condensed", size = 4,
    data = pooled %>%
      mutate(y = quantile(fit, 0.9), .by = c(model, var)) %>%
      filter(x %in% c(ozone_seq[length(ozone_seq)], tmax_seq[length(tmax_seq)]), model == "Full model") %>% select(x, y, var) %>%
      mutate(pv = c(cal_significance(models$fml_full, models$fml_no_tmax, type = "ozone"), cal_significance(models$fml_full, models$fml_no_ozone, type = "tmax")))
  ) +
  scale_y_percent("Change in yield") +
  scale_x_continuous(expand = expansion(mult = 0.01)) +
  scale_xsidey_continuous(labels = NULL, breaks = NULL) +
  scale_fill_manual(values = c("gray30", "#C15555", "#89AEE5")) +
  scale_color_manual(values = c("gray30", "#C15555", "#89AEE5")) +
  ggside(x.pos = "bottom", scales = "free") +
  theme_half_open(font_size = 16, font_family = "Roboto Condensed") +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.justification = 0.5,
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.title.x = element_blank(),
    ggside.panel.scale = 0.2
  )

ggsave("figures/response_comparison_pooled.pdf", p1, width = 2, height = 1, scale = 5)



# single crops ------------------------------------------------------------

data <- get_data() %>%
  mutate(logYield = log(Yield))

data <- data %>%
  pivot_longer(c(aot40, starts_with("tmax"), starts_with("AOD"), starts_with("sm"), starts_with("cloud"))) %>%
  mutate(var = str_remove(name, "[23]")) %>%
  inner_join(season, by = join_by(Item, prop, var)) %>%
  select(-c(prop, var)) %>%
  pivot_wider() %>%
  drop_na()

data <- data %>%
  nest(.by = Item) %>%
  arrange(factor(Item, order))

walk2(data$Item, data$data, function(aItem, adata) {
  models <- list(
    fml_full = logYield ~ aot40 + tmax + tmax2 + tmax3 + AOD + AOD2 + AOD3 + sm + sm2 + sm3 + cloud + cloud2 + cloud3 | ADM_CODE[year + year^2],
    fml_no_ozone = logYield ~ tmax + tmax2 + tmax3 + AOD + AOD2 + AOD3 + sm + sm2 + sm3 + cloud + cloud2 + cloud3 | ADM_CODE[year + year^2],
    fml_no_tmax = logYield ~ aot40 + AOD + AOD2 + AOD3 + sm + sm2 + sm3 + cloud + cloud2 + cloud3 | ADM_CODE[year + year^2]
  ) %>%
    map(feols, data = adata)

  # How to determine whether the estimated coeficiences are significantly different between two models? For example, a model is `lm(yield ~ tmax + I(tmax ^ 2) + precipitation + ozone, data = data)` and another is `lm(yield ~ tmax + I(tmax ^ 2) + precipitation, data = data)`, I want to know whether adding the ozone variable will significantly alter the coefficients of temperature (tmax) on yield. Provide the R code.

  ozone_seq <- seq(min(adata$aot40), max(adata$aot40), by = 100)
  tmax_seq <- seq(min(adata$tmax), max(adata$tmax), by = 0.1)

  pooled <- bind_rows(
    map_dfr(models, relpred,
      newdata = tibble(aot40 = ozone_seq),
      baseline = tibble(aot40 = rep(median(ozone_seq), length(ozone_seq))),
      x = ozone_seq, .id = "model"
    ) %>%
      filter(model != "fml_no_ozone") %>%
      mutate(var = "AOT40 (ppb h)"),
    map_dfr(models, relpred,
      newdata = tibble(tmax = tmax_seq, tmax2 = tmax^2, tmax3 = tmax^3),
      baseline = tibble(tmax = rep(median(tmax_seq), length(tmax_seq)), tmax2 = tmax^2, tmax3 = tmax^3),
      x = tmax_seq, .id = "model"
    ) %>%
      filter(model != "fml_no_tmax") %>%
      mutate(var = "Daily maximum temperature ($\\degree C$)")
  ) %>%
    mutate(across(c(fit, lwr, upr), expm1),
      model = factor(model, labels = c("Full model", "Full model without ozone", "Full model without temperature"))
    )

  labels <- pooled %>%
    filter(x %in% c(ozone_seq[length(ozone_seq)], tmax_seq[length(tmax_seq)])) %>%
    select(model, x, fit, var) %>%
    pivot_wider(names_from = model, values_from = fit) %>%
    pivot_longer(-c(x, var, `Full model`), values_drop_na = T)

  arrow <- arrow(
    length = unit(0.015, "npc"),
    ends = "last",
    type = "open"
  )

  p1 <- pooled %>%
    ggplot(aes(x = x, y = fit)) +
    # facet_wrap(~ TeX(var, output = "character"), labeller = label_parsed, scales = 'free', strip.position = 'bottom') + # slow
    facet_wrap(~ var, scales = "free", strip.position = "bottom") +
    geom_hline(yintercept = 0, linetype = "longdash", color = "grey") +
    geom_ribbon(aes(ymax = upr, ymin = lwr, fill = model), alpha = 0.2) +
    geom_line(aes(color = model), show.legend = F, linewidth = 1.5) +
    geom_xsidehistogram(aes(x = value), data = adata %>% select(tmax, aot40) %>% rename("Daily maximum temperature ($\\degree C$)" = tmax, "AOT40 (ppb h)" = aot40) %>% pivot_longer(everything(), names_to = "var"), fill = "grey", color = "black", inherit.aes = F, bins = 50) +
    geom_segment(aes(x = x * 1.02, xend = x * 1.02, y = value, yend = `Full model`), arrow = arrow, inherit.aes = F, data = labels) +
    # geom_text_repel(aes(x = x * 1.18, y = (value + `Full model`) / 2, label = label), data = labels %>% mutate(label = str_c(round_pettry(-(value - `Full model`) * 100), "%")), inherit.aes = F, min.segment.length = Inf, box.padding = 0, family = "Roboto Condensed", size = 4) +
    geom_text(aes(x = x, y = y, label = str_c("P value = ", round_pettry(pv))),
      inherit.aes = F, family = "Roboto Condensed", size = 4,
      data = pooled %>%
        mutate(y = quantile(fit, 0.9), .by = c(model, var)) %>%
        filter(x %in% c(ozone_seq[length(ozone_seq)], tmax_seq[length(tmax_seq)]), model == "Full model") %>% select(x, y, var) %>%
        mutate(pv = c(cal_significance(models$fml_full, models$fml_no_tmax, type = "ozone"), cal_significance(models$fml_full, models$fml_no_ozone, type = "tmax")))
    ) +
    scale_y_percent("Change in yield") +
    scale_x_continuous(expand = expansion(mult = 0.01)) +
    scale_xsidey_continuous(labels = NULL, breaks = NULL) +
    scale_fill_manual(values = c("gray30", "#C15555", "#89AEE5")) +
    scale_color_manual(values = c("gray30", "#C15555", "#89AEE5")) +
    ggside(x.pos = "bottom", scales = "free") +
    labs(title = aItem) +
    theme_half_open(font_size = 16, font_family = "Roboto Condensed") +
    theme(
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.justification = 0.5,
      strip.background = element_blank(),
      strip.placement = "outside",
      axis.title.x = element_blank(),
      ggside.panel.scale = 0.2,
      plot.title = element_text(
        size = 20, hjust = 0.5,
        margin = margin(0.5, 0, 0.5, 0, "lines")
      )
    )

  print(p1)
})
