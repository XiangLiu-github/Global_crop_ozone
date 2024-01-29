source("script/loadpackages.R")
source("script/loadfunctions.R")

shp <- read_rds("data/inputs/shp/GAUL.rds") %>%
  st_simplify(dTolerance = 0.1) %>%
  filter(!ADM_NAME %in% c("Antarctica"))

aot40 <- read_rds("data/outputs/ozone/extracted_aot40.rds") %>% unnest()
tmax <- read_rds("data/outputs/tmax/extracted_tmax.rds") %>% unnest()

season <- get_season() %>%
  select(Item, var, prop)
order <- get_order()

data <- inner_join(
  season %>% filter(var == "aot40") %>% inner_join(aot40) %>% select(-c(var, prop)),
  season %>% filter(var == "tmax") %>% inner_join(tmax) %>% select(-c(var, prop))
)

feols(tmax ~ aot40 | ADM_NAME[year + year^2], data, split = ~Item) %>%
  map_dfr(r2, type = "wr2", .id = "crop") %>%
  separate(crop, c(NA, NA, "crop"), ":") %>%
  mutate(
    crop = str_trim(crop),
    crop = factor(crop, get_order())
  ) %>%
  arrange(crop)

p <- data %>%
  mutate(across(c(tmax, aot40), ~ detrend(.x, year)), .by = c(Item, ADM_CODE, ADM_NAME, FAO_CODE)) %>%
  nest(.by = c(Item, ADM_CODE, ADM_NAME, FAO_CODE)) %>%
  mutate(data = map(data, ~ tidy(cor.test(.x$tmax, .x$aot40)))) %>%
  unnest() %>%
  mutate(sig = ifelse(p.value <= 0.05, "grey30", "grey90")) %>%
  inner_join(shp, .) %>%
  nest(.by = Item) %>%
  arrange(factor(Item, order)) %>%
  pmap(function(Item, data) {
    plot_data <- overlap_CN(data)

    p1 <- plot_data %>%
      ggplot(aes(fill = estimate, color = sig)) +
      geom_sf(color = NA, fill = "grey", size = 0.2, data = shp, inherit.aes = F) +
      geom_sf(size = 0.1) +
      labs(title = Item) +
      scale_fill_gradientn(
        name = "Correlation coefficient",
        limits = c(-0.8, 0.8),
        oob = squish,
        colors = NatParksPalettes::natparks.pals("Arches", n = 11),
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
      ggplot(aes(y = estimate)) +
      geom_histogram(aes(x = after_stat(..count..) / max(after_stat(..count..)), fill = ..y..), bins = 20, color = "black", show.legend = F) +
      # geom_vline(xintercept = median(data$cor), linetype = 'longdash', linewidth = 1) +
      # geom_text_repel(aes(y = 1.4, label = round_pettry(cor)), data = . %>% st_drop_geometry() %>% summarise(cor = round(median(cor), 2)), nudge_x = 0.2, family = "Roboto Condensed", size = 5, min.segment.length = Inf) +
      scale_fill_gradientn(
        name = "Correlation coefficient",
        limits = c(-0.8, 0.8),
        oob = squish,
        colors = NatParksPalettes::natparks.pals("Arches", n = 10)
      ) +
      theme_half_open(font_size = 16, font_family = "Roboto Condensed") +
      scale_x_continuous(expand = expansion()) +
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

ggsave("figures/correlation_map_ozone_tmax.pdf", p, width = 1, height = 1.3, scale = 20)

p3 <- data %>%
  filter(ADM_NAME %in% c("United States of America", "France", "China, mainland"), Item %in% c('Maize', "Wheat", 'Soybean', 'Rapeseed')) %>%
  select(ADM_NAME, Item, year, aot40, tmax) %>%
  pivot_longer(-c(ADM_NAME, Item, year)) %>%
  mutate(value = fscale(detrend(value, year)), .by = c(ADM_NAME, Item, name)) %>%
  ggplot(aes(x = year, y = value, color = factor(name, labels = c("Ozone", "Temperature")))) +
  facet_grid(vars(factor(ADM_NAME, c("United States of America", "France", "China, mainland"), c("USA", "France", "China"))), vars(factor(Item, order))) +
  geom_hline(yintercept = 0, linetype = "longdash", color = "grey") +
  geom_point_s(size = 2) +
  geom_line() +
  geom_text_repel(aes(x = 2019, y = 3, label = label),
    data = . %>%
      pivot_wider() %>%
      summarise(
        p.value = cor.test(aot40, tmax)$p.value,
        cor = cor.test(aot40, tmax)$estimate,
        .by = c(ADM_NAME, Item)
      ) %>%
      mutate(
        p.value = case_when(
          p.value <= 0.05 ~ "**",
          p.value <= 0.1 ~ "*",
          TRUE ~ ""
        ),
        label = str_c(round_pettry(cor), p.value)
      ),
    inherit.aes = F, family = "Roboto Condensed", size = 5
  ) +
  scale_x_continuous("Year", breaks = c(2003, 2011, 2019)) +
  scale_y_continuous("Normalized & detrended anomaly") +
  scale_color_manual(values = c("#008500", "#FF8500")) +
  theme_half_open(font_size = 16, font_family = "Roboto Condensed") +
  background_grid(major = "x", minor = "x") +
  theme(
    strip.background = element_blank(),
    legend.position = "bottom",
    legend.justification = 0.5,
    legend.title = element_blank()
  )

ggsave("figures/correlation_line_ozone_tmax.pdf", p3, width = 4.5, height = 2, scale = 3)
