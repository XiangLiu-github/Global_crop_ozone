source("script/loadpackages.R")
source("script/loadfunctions.R")

# must in largemem

coefs <- qread("data/outputs/boot_coefs.qs")

region <- get_shp()

shp <- region %>%
  st_simplify(dTolerance = 0.1) %>%
  filter(!ADM_NAME %in% c("Antarctica"))

weights <- qread("data/inputs/cropland/weights.qs") %>%
  mutate(weights = map(weights, function(aweight) {
    aweight %>%
      `names<-`("weights") %>%
      as.data.frame(xy = T) %>%
      filter(weights != 0)
  }, .progress = T))

extracted_weights <- qread("data/inputs/cropland/weights.qs") %>%
  mutate(weights = map(weights, function(aweight) {
    region %>%
      st_drop_geometry() %>%
      mutate(harvest_area = exact_extract(aweight, region, "sum", progress = F))
  }, .progress = T))

# ozone ---------------------------------------------------------------------

aot40_hist <- read_rds("data/outputs/ozone/extracted_aot40.rds") %>% unnest()
aot40_nat <- read_rds("data/outputs/ozone/extracted_aot40_nat.rds") %>%
  unnest() %>%
  unnest() %>%
  rename_with(~ str_c(.x, "_nat"), starts_with("aot40"))

# nat is counterfactual
aot40 <- inner_join(aot40_hist, aot40_nat) %>%
  mutate(aot40 = aot40 - aot40_nat, .keep = "unused") %>%
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
      fgroup_by(model, ADM_CODE, ADM_NAME, FAO_CODE, ISO, NAME, name, harvest_area) %>%
      fsummarise(value = mean(value)) %>%
      group_by(ADM_CODE, ADM_NAME, FAO_CODE, ISO, NAME, harvest_area) %>%
      median_qi(value, .width = c(0.66, 0.9))

    # weighted mean by region is needed
    fortrend <- result %>%
      fgroup_by(year, name) %>%
      fsummarise(value = weighted.mean(value, harvest_area)) %>%
      group_by(year) %>%
      median_qi(value, .width = c(0.66, 0.9))

    lst(formap, fortrend)
  }, .progress = T), .keep = "unused")

# aot40 %>%
#   unnest_wider(data) %>%
#   select(-formap) %>%
#   unnest() %>%
#   ggplot(aes(x = year, y = value)) +
#   facet_wrap(~ Item, scales = 'free_y') +
#   geom_lineribbon(aes(ymax = .upper, ymin = .lower)) +
#   scale_y_percent() +
#   scale_fill_manual(values = ozone_pal[c(1,3)])

# should calculate at pixel level!!!! ozone is fine, because it is not nonlinear effect
# # tmax --------------------------------------------------------------------
#
# tmax_hist = read_rds('data/outputs/tmax/extracted_tmax.rds') %>% unnest()
# tmax_nat = read_rds('data/outputs/tmax/extracted_tmax_nat.rds') %>% unnest() %>% unnest() %>% rename_with(~ str_c(.x, '_nat'), starts_with('tmax'))
#
# tmax = inner_join(tmax_hist, tmax_nat) %>%
#   mutate(tmax = tmax - tmax_nat, tmax2 = tmax2 - tmax2_nat, tmax3 = tmax3 - tmax3_nat, .keep = 'unused') %>%
#   nest(.by = Item)
#
# tmax = list(coefs, tmax, extracted_weights) %>%
#   reduce(inner_join) %>%
#   mutate(data = pmap(list(coefs, data, weights), function(acoef, adata, aweight){
#
#     acoef = acoef %>% select(starts_with('tmax'))
#
#     rel_X = adata %>%
#       select(starts_with('tmax'))
#
#     result = tcrossprod(as.matrix(rel_X), as.matrix(acoef)) %>%
#       expm1() %>%
#       as_tibble() %>%
#       bind_cols(adata %>% select(-starts_with('tmax'))) %>%
#       pivot_longer(num_range('V', 1:1e3)) %>%
#       inner_join(aweight, by = join_by(ADM_CODE, ADM_NAME, FAO_CODE))
#
#     formap = result %>%
#       fgroup_by(model, ADM_CODE, ADM_NAME, FAO_CODE, ISO, NAME, name, harvest_area) %>%
#       fsummarise(value = mean(value)) %>%
#       group_by(ADM_CODE, ADM_NAME, FAO_CODE, ISO, NAME, harvest_area) %>%
#       median_qi(value, .width = c(0.66, 0.9))
#
#     # weighted mean by region is needed
#     fortrend = result %>%
#       fgroup_by(year, name, model) %>%
#       fsummarise(value = weighted.mean(value, harvest_area)) %>%
#       group_by(year) %>%
#       median_qi(value, .width = c(0.66, 0.9))
#
#     lst(formap, fortrend)
#
#   }, .progress = T), .keep = 'unused')
#
# tmax %>%
#   unnest_wider(data) %>%
# select(-formap) %>%
# unnest() %>%
# ggplot(aes(x = year, y = value)) +
# facet_wrap(~ Item, scales = 'free_y') +
# geom_lineribbon(aes(ymax = .upper, ymin = .lower)) +
# scale_y_percent() +
# scale_fill_manual(values = temperature_pal[c(6,8)])
#
# p = tmax %>%
#   unnest_wider(data) %>%
#   select(-fortrend) %>%
#   arrange(factor(Item, get_order())) %>%
#   pmap(function(Item, formap){
#
#     formap = formap %>%
#       filter(.width == 0.66) %>%
#       inner_join(shp, ., by = join_by(ADM_CODE, ADM_NAME, FAO_CODE))
#
#     limits = quantile(formap %>% pull(value), c(0.1, 0.9)) %>%
#       ggpmisc::symmetric_limits()
#
#     formap %>%
#       ggplot(aes(fill = value)) +
#       geom_sf(color = NA, fill = 'grey', size = 0.2, data = shp, inherit.aes = F) +
#       geom_sf(color = 'grey30', size = 0.1) +
#       labs(title = Item) +
#       scale_fill_gradientn(name = 'T',
#                            limits = limits,
#                            label = label_percent(),
#                            oob = squish,
#                            colors = rev(temperature_pal),
#                            breaks = breaks_extended(n = 5),
#                            guide = guide_coloursteps(title.position = "top")) +
#       theme_void(base_family = "Roboto Condensed", base_size = 16) +
#       theme(plot.title = element_text(size = 25, hjust = 0.5,
#                                       margin = margin(0.5,0,0.5,0, "lines")),
#             legend.position = 'bottom',
#             legend.direction = 'horizontal',
#             legend.key.width = unit(5, 'lines'),
#             legend.key.height = unit(0.5, 'lines'),
#             legend.title.align = 0.5) +
#       coord_sf(expand = F, crs = "+proj=robin")
#
#   }) %>%
#   wrap_plots()
#
# ggsave('test.pdf', p, width = 2, height = 1, scale = 20)
#

# tmax new ----------------------------------------------------------------

tmax_hist <- qread("data/outputs/tmax/gs_tmax.qs")
tmax_nat <- qread("data/outputs/tmax/gs_tmax_nat.qs")

season <- get_season() %>%
  filter(var == "tmax") %>%
  select(Item, prop)

tmax_hist <- inner_join(season, tmax_hist) %>%
  mutate(tmax = map2(prop, tmax, function(aprop, atmax) {
    aprop <- str_c("\\.", aprop)

    atmax <- atmax %>%
      subset(names(atmax)[str_detect(names(atmax), aprop)])

    names(atmax) <- str_remove(names(atmax), aprop)

    return(atmax)
  })) %>%
  select(-prop)

tmax <- tmax_nat %>%
  unnest() %>%
  select(-prop) %>%
  rename(nat = tmax) %>%
  inner_join(tmax_hist %>% rename(hist = tmax)) %>%
  nest(.by = Item)

library(dtplyr)
tmax <- tmax %>%
  inner_join(weights) %>%
  inner_join(coefs) %>%
  mutate(data = pmap(list(data, weights, coefs), function(adata, aweight, acoef) {
    acoef <- acoef %>% select(starts_with("tmax"))

    res <- adata %>%
      # head(1) %>%
      mutate(temp = map2(nat, hist, function(anat, ahist) {
        stopifnot(names(ahist) == names(anat))

        df <- (ahist - anat) %>%
          `names<-`(names(ahist)) %>%
          as.data.frame(xy = T) %>%
          pivot_longer(-c(x, y),
            names_to = c("name", "year"), names_sep = "\\.", names_transform = list(year = as.integer),
            values_drop_na = T
          ) %>%
          pivot_wider() %>%
          inner_join(aweight, by = c("x", "y"))

        rel_X <- df %>%
          select(starts_with("tmax"))

        result <- tcrossprod(as.matrix(rel_X), as.matrix(acoef)) %>%
          expm1() %>%
          as.data.table() %>%
          lazy_dt() %>%
          cbind(df %>% select(-starts_with("tmax"))) %>%
          pivot_longer(num_range("V", 1:1e3)) %>%
          as.data.table()

        formap <- result %>%
          lazy_dt() %>%
          group_by(x, y, name, weights) %>%
          summarise(value = mean(value), .groups = "drop") %>%
          as_tibble()

        fortrend <- result %>%
          lazy_dt() %>%
          group_by(year, name) %>%
          summarise(value = weighted.mean(value, weights), .groups = "drop") %>%
          as_tibble()

        lst(formap, fortrend)
      }, .progress = T), .keep = "unused") %>%
      unnest_wider(temp)

    formap <- res %>%
      select(model, formap) %>%
      unnest() %>%
      lazy_dt() %>%
      group_by(x, y) %>%
      summarise(value = median(value), .groups = "drop") %>%
      as_tibble()

    fortrend <- res %>%
      select(model, fortrend) %>%
      unnest() %>%
      group_by(year) %>%
      median_qi(value, .width = c(0.66, 0.9))

    lst(formap, fortrend)
  }, .progress = T), .keep = "unused")

qsave(lst(aot40, tmax), "data/outputs/impact_yield.qs")

a$tmax %>%
  unnest_wider(data) %>%
  select(-formap) %>%
  unnest() %>%
  ggplot(aes(x = year, y = value)) +
  facet_wrap(~Item, scales = "free_y") +
  geom_lineribbon(aes(ymax = .upper, ymin = .lower)) +
  scale_y_percent() +
  scale_fill_manual(values = temperature_pal[c(6, 8)])
