source("script/loadpackages.R")
source("script/loadfunctions.R")

impacts <- qread("data/outputs/impact_yield_future.qs")
impacts_ <- qread("data/outputs/impact_yield_future_no.qs")

# maps --------------------------------------------------------------------

shp <- get_shp() %>%
  st_simplify(dTolerance = 0.1) %>%
  filter(!ADM_NAME %in% c("Antarctica"))

impacts$aot40 %>%
  unnest_wider(data) %>%
  select(-fortrend) %>%
  unnest() %>%
  mutate(sig = ifelse(.lower * .upper >= 0, "grey30", "grey90")) %>%
  nest(.by = c(Item, indx)) %>%
  arrange(factor(Item, get_order())) %>%
  nest(.by = indx, .key = "ddata") %>%
  pwalk(function(indx, ddata) {
    p <- pmap(ddata, function(Item, data) {
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
    
    ggsave(str_c("figures/impacts_map_aot40_", indx, ".pdf"), p, width = 1, height = 1.3, scale = 20)
  })

impacts$tmax %>%
  unnest_wider(data) %>%
  select(-fortrend) %>%
  unnest() %>%
  mutate(sig = ifelse(.lower * .upper >= 0, "grey30", "grey90")) %>%
  nest(.by = c(Item, indx)) %>%
  arrange(factor(Item, get_order())) %>%
  nest(.by = indx, .key = "ddata") %>%
  pwalk(function(indx, ddata) {
    p <- pmap(ddata, function(Item, data) {
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
          colors = rev(temperature_trend_pal),
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
          colors = rev(temperature_trend_pal)
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
    
    ggsave(str_c("figures/impacts_map_tmax_", indx, ".pdf"), p, width = 1, height = 1.3, scale = 20)
  })

list(
  impacts$aot40,
  impacts_$aot40
) %>%
  map2_dfr(c("Full", "No"), function(adata, aname) {
    adata %>%
      unnest_wider(data) %>%
      select(-fortrend) %>%
      unnest() %>%
      filter(.width == uncertainty[2]) %>%
      select(-c(.width, .point, .interval)) %>%
      mutate(name = aname)
  }) %>%
  pivot_longer(c(value, .upper, .lower),
               names_to = "type"
  ) %>%
  pivot_wider() %>%
  mutate(delta = Full - No, .keep = "unused") %>%
  pivot_wider(values_from = delta, names_from = type) %>%
  mutate(sig = ifelse(.lower * .upper >= 0, "grey30", "grey90")) %>%
  nest(.by = c(Item, indx)) %>%
  arrange(factor(Item, get_order())) %>%
  nest(.by = indx, .key = "ddata") %>%
  pwalk(function(indx, ddata) {
    p <- pmap(ddata, function(Item, data) {
      data <- data %>%
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
          name = "Change in ozone impact",
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
    
    ggsave(str_c("figures/impacts_delta_map_aot40_", indx, ".pdf"), p, width = 1, height = 1.3, scale = 20)
  })

list(
  impacts$tmax,
  impacts_$tmax
) %>%
  map2_dfr(c("Full", "No"), function(adata, aname) {
    adata %>%
      unnest_wider(data) %>%
      select(-fortrend) %>%
      unnest() %>%
      filter(.width == uncertainty[2]) %>%
      select(-c(.width, .point, .interval)) %>%
      mutate(name = aname)
  }) %>%
  pivot_longer(c(value, .upper, .lower),
               names_to = "type"
  ) %>%
  pivot_wider() %>%
  mutate(delta = Full - No, .keep = "unused") %>%
  pivot_wider(values_from = delta, names_from = type) %>%
  mutate(sig = ifelse(.lower * .upper >= 0, "grey30", "grey90")) %>%
  nest(.by = c(Item, indx)) %>%
  arrange(factor(Item, get_order())) %>%
  nest(.by = indx, .key = "ddata") %>%
  pwalk(function(indx, ddata) {
    p <- pmap(ddata, function(Item, data) {
      data <- data %>%
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
          name = "Change in temperature impact",
          limits = limits,
          label = label_percent(),
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
        ggplot(aes(y = value)) +
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
    
    ggsave(str_c("figures/impacts_delta_map_tmax_", indx, ".pdf"), p, width = 1, height = 1.3, scale = 20)
  })

# overall, box ------------------------------------------------------------

# po = bind_rows(
#   impacts$aot40 %>%
#     unnest_wider(data) %>%
#     select(-formap) %>%
#     unnest() %>%
#     summarise(across(c(value, .lower, .upper), mean), .by = c(Item, indx, .width)) %>%
#     mutate(type = 'Full model'),
#   impacts_$aot40 %>%
#     unnest_wider(data) %>%
#     select(-formap) %>%
#     unnest() %>%
#     summarise(across(c(value, .lower, .upper), mean), .by = c(Item, indx, .width)) %>%
#     mutate(type = 'Full model without temperature')
# ) %>%
#   ggplot(aes(x = indx, y = value, color = type)) +
#   facet_wrap(~ factor(Item, get_order()), scales = 'free_y') +
#   geom_hline(yintercept = 0, linetype = 'longdash', color = 'grey') +
#   geom_pointinterval(aes(ymax = .upper, ymin = .lower), position = position_dodge(width = 0.5)) +
#   scale_y_percent('Change in yield') +
#   scale_color_manual(values = c('gray30', '#89AEE5')) +
#   theme_half_open(font_size = 16, font_family = "Roboto Condensed") +
#   background_grid(major = 'x') +
#   theme(strip.background = element_blank(),
#         legend.title = element_blank(),
#         legend.position = c(0.6, 0.1),
#         axis.title.x = element_blank())
#
# ggsave('figures/impacts_pointinterval_aot40.pdf', po, width = 2, height = 1, scale = 7)
#
# pt = bind_rows(
#   impacts$tmax %>%
#     unnest_wider(data) %>%
#     select(-formap) %>%
#     unnest() %>%
#     summarise(across(c(value, .lower, .upper), mean), .by = c(Item, indx, .width)) %>%
#     mutate(type = 'Full model'),
#   impacts_$tmax %>%
#     unnest_wider(data) %>%
#     select(-formap) %>%
#     unnest() %>%
#     summarise(across(c(value, .lower, .upper), mean), .by = c(Item, indx, .width)) %>%
#     mutate(type = 'Full model without ozone')
# ) %>%
#   ggplot(aes(x = indx, y = value, color = type)) +
#   facet_wrap(~ factor(Item, get_order()), scales = 'free_y') +
#   geom_hline(yintercept = 0, linetype = 'longdash', color = 'grey') +
#   geom_pointinterval(aes(ymax = .upper, ymin = .lower), position = position_dodge(width = 0.5)) +
#   scale_y_percent('Change in yield') +
#   scale_color_manual(values = c('gray30', "#C15555")) +
#   theme_half_open(font_size = 16, font_family = "Roboto Condensed") +
#   background_grid(major = 'x') +
#   theme(strip.background = element_blank(),
#         legend.title = element_blank(),
#         legend.position = c(0.6, 0.1),
#         axis.title.x = element_blank())
#
# ggsave('figures/impacts_pointinterval_tmax.pdf', pt, width = 2, height = 1, scale = 7)


pp <- bind_rows(
  impacts$aot40 %>%
    unnest_wider(data) %>%
    select(-formap) %>%
    unnest() %>%
    summarise(across(c(value, .lower, .upper), mean), .by = c(Item, indx, .width)) %>%
    mutate(
      var = "Ozone impact",
      type = "Full model"
    ),
  impacts_$aot40 %>%
    unnest_wider(data) %>%
    select(-formap) %>%
    unnest() %>%
    summarise(across(c(value, .lower, .upper), mean), .by = c(Item, indx, .width)) %>%
    mutate(
      var = "Ozone impact",
      type = "Full model without temperature"
    ),
  impacts$tmax %>%
    unnest_wider(data) %>%
    select(-formap) %>%
    unnest() %>%
    summarise(across(c(value, .lower, .upper), mean), .by = c(Item, indx, .width)) %>%
    mutate(
      var = "Temperature impact",
      type = "Full model"
    ),
  impacts_$tmax %>%
    unnest_wider(data) %>%
    select(-formap) %>%
    unnest() %>%
    summarise(across(c(value, .lower, .upper), mean), .by = c(Item, indx, .width)) %>%
    mutate(
      var = "Temperature impact",
      type = "Full model without ozone"
    )
) %>%
  ggplot(aes(x = indx, y = value, shape = var, color = type)) +
  facet_wrap(~ factor(Item, get_order()), ncol = 4, scales = "free_y") +
  geom_hline(yintercept = 0, linetype = "longdash", color = "grey") +
  geom_pointinterval(aes(ymax = .upper, ymin = .lower), position = position_dodge(width = 0.7)) +
  scale_y_percent("Change in yield") +
  scale_color_manual(values = c("gray30", "#C15555", "#89AEE5")) +
  scale_shape_manual(values = c(16, 15)) +
  theme_half_open(font_size = 16, font_family = "Roboto Condensed") +
  background_grid(major = "x") +
  theme(
    strip.background = element_blank(),
    legend.title = element_blank(),
    legend.position = c(0.5, 0.1),
    legend.box = "horizontal",
    axis.title.x = element_blank()
  )

ggsave("figures/impacts_pointinterval.pdf", pp, width = 1.2, height = 1, scale = 10)

# overall, text ------------------------------------------------------------

a = 
  bind_rows(
    impacts$aot40 %>%
      unnest_wider(data) %>%
      select(-formap) %>%
      unnest() %>%
      summarise(across(c(value, .lower, .upper), mean), .by = c(Item, indx, .width)) %>%
      mutate(
        var = "Ozone impact",
        type = "Full model"
      ),
    impacts_$aot40 %>%
      unnest_wider(data) %>%
      select(-formap) %>%
      unnest() %>%
      summarise(across(c(value, .lower, .upper), mean), .by = c(Item, indx, .width)) %>%
      mutate(
        var = "Ozone impact",
        type = "Full model without temperature"
      ),
    impacts$tmax %>%
      unnest_wider(data) %>%
      select(-formap) %>%
      unnest() %>%
      summarise(across(c(value, .lower, .upper), mean), .by = c(Item, indx, .width)) %>%
      mutate(
        var = "Temperature impact",
        type = "Full model"
      ),
    impacts_$tmax %>%
      unnest_wider(data) %>%
      select(-formap) %>%
      unnest() %>%
      summarise(across(c(value, .lower, .upper), mean), .by = c(Item, indx, .width)) %>%
      mutate(
        var = "Temperature impact",
        type = "Full model without ozone"
      )
  )

a %>% 
  filter(Item %in% c('Maize', 'Wheat'), indx == '2045-2049', type == 'Full model', .width == 0.95) %>% 
  mutate(across(c(value, .lower, .upper), ~ round_pettry(.x * 100))) %>% 
  arrange(Item) %>% 
  mutate(text = str_c('\\qty{', value, '}{\\percent} (\\qty{95}{\\percent} CI: \\qtyrange{', .lower, '}{', .upper, '}{\\percent})')) %>% 
  pull() %>% 
  walk(writeLines)

a %>% 
  filter(Item %in% c('Maize', 'Soybean'), indx == '2045-2049', var == 'Temperature impact', .width == 0.95) %>% 
  summarise(across(c(value, .lower, .upper), ~ round_pettry(-diff(.x) * 100)), .by = Item) %>% 
  mutate(text = str_c('\\qty{', value, '}{\\percent} (\\qty{95}{\\percent} CI: \\qtyrange{', .upper, '}{', .lower, '}{\\percent})')) %>% 
  pull() %>% 
  walk(writeLines)











