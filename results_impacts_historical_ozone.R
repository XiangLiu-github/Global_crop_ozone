source("script/loadpackages.R")
source("script/loadfunctions.R")

coefs <- qread("data/outputs/boot_coefs.qs")

region <- get_shp()

shp <- region %>%
  st_simplify(dTolerance = 0.1) %>%
  filter(!ADM_NAME %in% c("Antarctica"))

extracted_weights <- qread("data/inputs/cropland/weights.qs") %>%
  mutate(weights = map(weights, function(aweight) {
    region %>%
      st_drop_geometry() %>%
      mutate(harvest_area = exact_extract(aweight, region, "sum", progress = F))
  }, .progress = T))

season <- get_season() %>%
  select(Item, var, prop) %>%
  filter(var == "aot40") %>%
  select(-var)

# ozone ---------------------------------------------------------------------

aot40_hist <- read_rds("data/outputs/ozone/extracted_aot40.rds") %>%
  unnest() %>%
  inner_join(season)

# counterfactual is zero
aot40 <- aot40_hist %>%
  nest(.by = Item)

aot40 <- list(coefs, aot40, extracted_weights) %>%
  reduce(inner_join) %>%
  mutate(data = pmap(list(coefs, data, weights), function(acoef, adata, aweight) {
    acoef <- acoef %>% select(starts_with("aot40"))

    rel_X <- adata %>%
      select(starts_with("aot40"))

    result <- tcrossprod(as.matrix(rel_X), as.matrix(acoef)) %>%
      expm1() %>%
      as_tibble() %>%
      bind_cols(adata %>% select(-starts_with("aot40"))) %>%
      pivot_longer(num_range("V", 1:1e3)) %>%
      inner_join(aweight, by = join_by(ADM_CODE, ADM_NAME, FAO_CODE))

    formap <- result %>%
      fgroup_by(ADM_CODE, ADM_NAME, FAO_CODE, ISO, NAME, name, harvest_area) %>%
      fsummarise(value = mean(value)) %>%
      group_by(ADM_CODE, ADM_NAME, FAO_CODE, ISO, NAME, harvest_area) %>%
      median_qi(value, .width = uncertainty)

    # weighted mean by region is needed
    fortrend <- result %>%
      fgroup_by(year, name) %>%
      fsummarise(value = weighted.mean(value, harvest_area)) %>%
      group_by(year) %>%
      median_qi(value, .width = uncertainty)

    lst(formap, fortrend)
  }, .progress = T), .keep = "unused")

qsave(aot40, "data/outputs/impact_historical_ozone.qs")

aot40 <- qread("data/outputs/impact_historical_ozone.qs")

p1 <- aot40 %>%
  unnest_wider(data) %>%
  select(-fortrend) %>%
  unnest() %>%
  mutate(sig = ifelse((.lower * .upper >= 0), "grey30", "grey90")) %>%
  nest(.by = c(Item)) %>%
  arrange(factor(Item, get_order())) %>%
  pmap(function(Item, data) {
    data <- data %>%
      filter(.width == uncertainty[2]) %>%
      inner_join(shp, ., by = join_by(ADM_CODE, ADM_NAME, FAO_CODE, ISO, NAME))

    plot_data <- overlap_CN(data)

    limits <- quantile(data %>% pull(value), c(0.05, 0.95)) %>%
      ggpmisc::symmetric_limits()

    p1 <- plot_data %>%
      ggplot(aes(fill = value, color = sig)) +
      geom_sf(color = NA, fill = "grey", size = 0.2, data = shp, inherit.aes = F) +
      geom_sf(size = 0.1) +
      labs(title = Item) +
      scale_fill_gradientn(
        name = "Change in yield",
        limits = limits,
        label = label_percent(),
        oob = squish,
        colors = rev(ozone_trend_pal),
        breaks = breaks_extended(n = 5),
        guide = guide_coloursteps(title.position = "top")
      ) +
      scale_color_identity() +
      theme_void(base_family = "Roboto Condensed", base_size = 16) +
      theme(
        plot.title = element_text(
          size = 25, hjust = 0.5,
          margin = margin(0.5, 0, 0.5, 0, "lines")
        ),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.width = unit(5, "lines"),
        legend.key.height = unit(0.5, "lines"),
        legend.title.align = 0.5
      ) +
      coord_sf(expand = F, crs = "+proj=robin")

    p2 <- data %>%
      ggplot(aes(y = value)) +
      geom_histogram(aes(x = after_stat(..count..) / max(after_stat(..count..)), fill = ..y..), bins = 20, color = "black", show.legend = F) +
      # geom_vline(xintercept = median(data$cor), linetype = 'longdash', linewidth = 1) +
      # geom_text_repel(aes(y = 1.4, label = round_pettry(cor)), data = . %>% st_drop_geometry() %>% summarise(cor = round(median(cor), 2)), nudge_x = 0.2, family = "Roboto Condensed", size = 5, min.segment.length = Inf) +
      scale_fill_gradientn(
        limits = limits,
        oob = squish,
        colors = rev(ozone_trend_pal)
      ) +
      theme_half_open(font_size = 16, font_family = "Roboto Condensed") +
      scale_x_continuous(expand = expansion()) +
      scale_y_continuous(labels = label_percent()) +
      theme(
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank()
      )

    p1 + inset_element(p2, 0, 0, 0.25, 0.6)
  }) %>%
  wrap_plots(ncol = 3)

ggsave("figures/impacts_map_aot40_historical.pdf", p1, width = 1, height = 1.3, scale = 20)

p2 <- aot40 %>%
  unnest_wider(data) %>%
  select(-formap) %>%
  inner_join(color_tbl) %>%
  unnest() %>%
  ggplot(aes(x = year, y = value, ymax = .upper, ymin = .lower, fill = color, color = color)) +
  facet_wrap(~ factor(Item, get_order()), ncol = 3, scales = "free_y") +
  geom_hline(yintercept = 0, linetype = "longdash", color = "grey") +
  geom_ribbon(data = . %>% filter(.width == uncertainty[2]), alpha = 0.2, lineend = "round") +
  geom_ribbon(data = . %>% filter(.width == uncertainty[1]), alpha = 0.6, lineend = "round") +
  ggborderline::geom_borderline(data = . %>% filter(.width == uncertainty[1]), linewidth = 1.5, lineend = "round", bordercolour = "black") +
  stat_poly_eq(use_label(c("eq", "P")), data = . %>% filter(.width == uncertainty[1]) %>% mutate(value = value * 100), color = "black", label.y = 0.9, family = "Roboto Condensed") +
  scale_color_identity() +
  scale_fill_identity() +
  scale_y_percent(name = "Change in Yield", expand = expansion(mult = c(0.1, 0.3))) +
  scale_x_continuous("Year", breaks = c(2003, 2011, 2019)) +
  theme_half_open(font_size = 16, font_family = "Roboto Condensed") +
  background_grid(major = "x", minor = "x") +
  theme(strip.background = element_blank())

ggsave("figures/impacts_trend_aot40_historical.pdf", p2, width = 1, height = 1.3, scale = 12)
