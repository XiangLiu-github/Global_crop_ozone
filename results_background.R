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
  drop_na()

shp <- read_rds("data/inputs/shp/GAUL.rds") %>%
  st_simplify(dTolerance = 0.1) %>%
  filter(!ADM_NAME %in% c("Antarctica"))

datam <- data %>%
  mutate(Yield = log(Yield), aot40 = aot40 / 1e3, AOD = AOD * 10) %>%
  summarise(across(c(Yield, aot40, AOD, tmax, sm, cloud, Area_harvested), fmean),
    .by = c(FAO_CODE, Item, ADM_NAME, ADM_CODE)
  )

data_nested <- datam %>%
  inner_join(shp, ., by = join_by(ADM_CODE, ADM_NAME, FAO_CODE)) %>%
  nest(.by = Item) %>%
  arrange(factor(Item, order))

pwalk(
  tibble(
    avar = c("Yield", "aot40", "AOD", "tmax", "sm", "cloud"),
    aname = list("log(Yield)", "AOT40 (ppm h)", TeX("AOD $\\times$ 10"), TeX("Daily maximum temperature ($\\degree C$)"), "Soil moisture (%)", "COD"),
    apal = list(yield_pal, ozone_pal, aerosol_pal, temperature_pal, moisture_pal, cloud_pal)
  ),
  function(avar, aname, apal) {
    patch <- data_nested %>%
      pmap(function(Item, data) {
        limits <- quantile(data %>% pull(avar), c(0.1, 0.9))

        plot_data <- overlap_CN(data)

        p1 <- plot_data %>%
          ggplot(aes(fill = !!sym(avar))) +
          geom_sf(color = NA, fill = "grey", size = 0.2, data = shp, inherit.aes = F) +
          geom_sf(color = "grey30", size = 0.1) +
          labs(title = Item) +
          scale_fill_gradientn(
            name = aname,
            limits = limits,
            oob = squish,
            colors = apal,
            breaks = breaks_extended(n = 5),
            guide = guide_coloursteps(title.position = "top")
          ) +
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
          ggplot(aes(y = !!sym(avar))) +
          geom_histogram(aes(x = after_stat(..count..) / max(after_stat(..count..)), fill = ..y..), bins = 20, color = "black", show.legend = F) +
          # geom_vline(xintercept = median(data$cor), linetype = 'longdash', linewidth = 1) +
          # geom_text_repel(aes(y = 1.4, label = round_pettry(cor)), data = . %>% st_drop_geometry() %>% summarise(cor = round(median(cor), 2)), nudge_x = 0.2, family = "Roboto Condensed", size = 5, min.segment.length = Inf) +
          scale_fill_gradientn(
            name = aname,
            limits = limits,
            oob = squish,
            colors = apal
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

    ggsave(str_c("figures/background_map_", avar, ".pdf"), patch, width = 1, height = 1.3, scale = 20)
  },
  .progress = T
)

patch <- bind_rows(
  datam %>%
    slice_max(Area_harvested, n = 3, by = Item),
  datam %>%
    slice_min(Area_harvested, n = 3, by = Item)
) %>%
  select(FAO_CODE, Item, ADM_NAME, ADM_CODE) %>%
  inner_join(data) %>%
  nest(.by = Item) %>%
  arrange(factor(Item, order)) %>%
  pmap(function(Item, data) {
    data %>%
      mutate(ISO = countrycode::countrycode(as.integer(ADM_CODE), "gaul", "iso3c", custom_match = c("147295" = "CHN", "147296" = "TWN"))) %>%
      ggplot(aes(x = year, y = log(Yield), color = ADM_NAME)) +
      geom_point(show.legend = F, size = 3) +
      geom_line(linewidth = 1, show.legend = F) +
      geom_label_repel(aes(label = ISO), data = . %>% slice_max(year, n = 1, by = ADM_NAME), show.legend = F, direction = "y", xlim = c(2019, Inf), max.overlaps = Inf, arrow = arrow, family = "Roboto Condensed", size = 5, min.segment.length = Inf) +
      scale_color_manual(values = WrensBookshelf::WB_brewer("MoreThanALittle")) +
      scale_x_continuous("Year", limits = c(NA, 2020), breaks = c(2003, 2011, 2019)) +
      labs(title = Item) +
      theme_half_open(font_size = 16, font_family = "Roboto Condensed") +
      background_grid(major = "x", minor = "x") +
      coord_cartesian(clip = "off") +
      theme(plot.title = element_text(
        size = 25, hjust = 0.5,
        margin = margin(0.5, 0, 0.5, 0, "lines")
      ))
  }) %>%
  wrap_plots(ncol = 3)

ggsave("figures/background_line.pdf", patch, width = 1, height = 1.3, scale = 16)

ozone_trend <- read_rds("data/outputs/ozone/extracted_aot40.rds") %>%
  inner_join(season %>% filter(var == "aot40")) %>%
  mutate(aot40 = map2(aot40, prop, function(adata, aprop) {
    adata %>%
      filter(prop == aprop) %>%
      mutate(aot40 = aot40 / 1e3) %>%
      nest(.by = c(ADM_CODE, ADM_NAME, FAO_CODE), .key = "fdata") %>%
      mutate(aot40 = map(fdata, function(afdata) {
        feols(aot40 ~ year, afdata) %>%
          tidy()
      }), .keep = "unused")
  }), .keep = "unused") %>%
  select(-c(var)) %>%
  unnest() %>%
  unnest() %>%
  filter(term == "year") %>%
  mutate(sig = fifelse(p.value < 0.1, "grey30", "grey90")) %>%
  nest(.by = Item) %>%
  arrange(factor(Item, order))

pt1 <- pmap(ozone_trend, function(Item, data) {
  data <- data %>%
    inner_join(shp, ., by = join_by(ADM_CODE, ADM_NAME, FAO_CODE))

  plot_data <- overlap_CN(data)

  limits <- quantile(data %>% pull(estimate), c(0.1, 0.9)) %>% symmetric_limits()

  p1 <- plot_data %>%
    ggplot(aes(fill = estimate, color = sig)) +
    geom_sf(color = NA, fill = "grey", size = 0.2, data = shp, inherit.aes = F) +
    geom_sf(size = 0.1) +
    labs(title = Item) +
    scale_fill_gradientn(
      name = TeX("AOT40 (ppm h $yr^{-1}$)"),
      limits = limits,
      oob = squish,
      colors = ozone_trend_pal,
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
      limits = limits,
      oob = squish,
      colors = ozone_trend_pal
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

ggsave(str_c("figures/background_map_trend_aot40.pdf"), pt1, width = 1, height = 1.3, scale = 20)

tmax_trend <- read_rds("data/outputs/tmax/extracted_tmax.rds") %>%
  inner_join(season %>% filter(var == "tmax")) %>%
  mutate(tmax = map2(tmax, prop, function(adata, aprop) {
    adata %>%
      nest(.by = c(ADM_CODE, ADM_NAME, FAO_CODE), .key = "fdata") %>%
      mutate(tmax = map(fdata, function(afdata) {
        feols(tmax ~ year, afdata) %>%
          tidy()
      }), .keep = "unused")
  }), .keep = "unused") %>%
  select(-c(var)) %>%
  unnest() %>%
  unnest() %>%
  filter(term == "year") %>%
  mutate(sig = fifelse(p.value < 0.1, "grey30", "grey90")) %>%
  nest(.by = Item) %>%
  arrange(factor(Item, order))

pt2 <- pmap(tmax_trend, function(Item, data) {
  data <- data %>%
    inner_join(shp, ., by = join_by(ADM_CODE, ADM_NAME, FAO_CODE))

  plot_data <- overlap_CN(data)

  limits <- quantile(data %>% pull(estimate), c(0.1, 0.9)) %>% symmetric_limits()

  p1 <- plot_data %>%
    ggplot(aes(fill = estimate, color = sig)) +
    geom_sf(color = NA, fill = "grey", size = 0.2, data = shp, inherit.aes = F) +
    geom_sf(size = 0.1) +
    labs(title = Item) +
    scale_fill_gradientn(
      name = TeX("Daily maximum temperature trend ($\\degree C$ $yr^{-1}$)"),
      limits = limits,
      oob = squish,
      colors = temperature_trend_pal,
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
      limits = limits,
      oob = squish,
      colors = temperature_trend_pal
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

ggsave(str_c("figures/background_map_trend_tmax.pdf"), pt2, width = 1, height = 1.3, scale = 20)
