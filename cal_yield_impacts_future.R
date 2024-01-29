source("script/loadpackages.R")
source("script/loadfunctions.R")

coefs <- qread("data/outputs/boot_coefs.qs")
coefs_noaot40 <- qread("data/outputs/boot_coefs_noaot40.qs")
coefs_notmax <- qread("data/outputs/boot_coefs_notmax.qs")

region <- get_shp()

shp <- region %>%
  st_simplify(dTolerance = 0.1) %>%
  filter(!ADM_NAME %in% c("Antarctica"))

# weights = qread('data/inputs/cropland/weights.qs') %>%
#   mutate(weights = map(weights, function(aweight){
#
#     aweight %>%
#       `names<-`('weights') %>%
#       as.data.frame(xy = T) %>%
#       filter(weights != 0)
#
#   }, .progress = T))

extracted_weights <- qread("data/inputs/cropland/weights.qs") %>%
  mutate(weights = map(weights, function(aweight) {
    region %>%
      st_drop_geometry() %>%
      mutate(harvest_area = exact_extract(aweight, region, "sum", progress = F))
  }, .progress = T))

aot40_hist <- read_rds("data/outputs/ozone/extracted_aot40.rds") %>% unnest()
aot40_fut <- read_rds("data/outputs/ozone/extracted_aot40_future.rds") %>%
  unnest() %>%
  unnest() %>%
  rename_with(~ str_c(.x, "_fut"), starts_with("aot40"))

tmax_hist <- read_rds("data/outputs/tmax/extracted_tmax.rds") %>% unnest()
tmax_fut <- read_rds("data/outputs/tmax/extracted_tmax_future.rds") %>%
  unnest() %>%
  unnest() %>%
  rename_with(~ str_c(.x, "_fut"), starts_with("tmax"))

# ozone ---------------------------------------------------------------------

aot40 <- inner_join(aot40_hist, aot40_fut) %>%
  mutate(aot40 = aot40_fut - aot40, .keep = "unused") %>%
  nest(.by = Item)

plan(multisession, workers = 3)

aot40 <- list(coefs, aot40, extracted_weights) %>%
  reduce(inner_join) %>%
  mutate(data = future_pmap(list(coefs, data, weights), function(acoef, adata, aweight) {
    adata <- adata %>%
      filter(year %in% 2005:2009)

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
      fgroup_by(model, indx, ADM_CODE, ADM_NAME, FAO_CODE, ISO, NAME, name, harvest_area) %>%
      fsummarise(value = mean(value)) %>%
      group_by(indx, ADM_CODE, ADM_NAME, FAO_CODE, ISO, NAME, harvest_area) %>%
      median_qi(value, .width = uncertainty)

    # weighted mean by region is needed
    fortrend <- result %>%
      fgroup_by(indx, year, name, model) %>%
      fsummarise(value = weighted.mean(value, harvest_area)) %>%
      group_by(indx, year) %>%
      median_qi(value, .width = uncertainty)

    lst(formap, fortrend)
  }, .progress = T), .keep = "unused")

plan(sequential)

# aot40 %>%
#   unnest_wider(data) %>%
#   select(-formap) %>%
#   unnest() %>%
#   ggplot(aes(x = year, y = value)) +
#   facet_wrap(~ Item, scales = 'free_y') +
#   geom_lineribbon(aes(ymax = .upper, ymin = .lower)) +
#   scale_y_percent() +
#   scale_fill_manual(values = ozone_pal[c(1,3)])

# tmax --------------------------------------------------------------------

tmax <- inner_join(tmax_hist, tmax_fut) %>%
  mutate(tmax = tmax_fut - tmax, tmax2 = tmax2_fut - tmax2, tmax3 = tmax3_fut - tmax3, .keep = "unused") %>%
  nest(.by = Item)

plan(multisession, workers = 3)

tmax <- list(coefs, tmax, extracted_weights) %>%
  reduce(inner_join) %>%
  mutate(data = future_pmap(list(coefs, data, weights), function(acoef, adata, aweight) {
    adata <- adata %>%
      filter(year %in% 2005:2009)

    acoef <- acoef %>% select(starts_with("tmax"))

    rel_X <- adata %>%
      select(starts_with("tmax"))

    result <- tcrossprod(as.matrix(rel_X), as.matrix(acoef)) %>%
      expm1() %>%
      as_tibble() %>%
      bind_cols(adata %>% select(-starts_with("tmax"))) %>%
      pivot_longer(num_range("V", 1:1e3)) %>%
      inner_join(aweight, by = join_by(ADM_CODE, ADM_NAME, FAO_CODE))

    formap <- result %>%
      fgroup_by(model, indx, ADM_CODE, ADM_NAME, FAO_CODE, ISO, NAME, name, harvest_area) %>%
      fsummarise(value = mean(value)) %>%
      group_by(indx, ADM_CODE, ADM_NAME, FAO_CODE, ISO, NAME, harvest_area) %>%
      median_qi(value, .width = uncertainty)

    # weighted mean by region is needed
    fortrend <- result %>%
      fgroup_by(indx, year, name, model) %>%
      fsummarise(value = weighted.mean(value, harvest_area)) %>%
      group_by(indx, year) %>%
      median_qi(value, .width = uncertainty)

    lst(formap, fortrend)
  }, .progress = T), .keep = "unused")

plan(sequential)

qsave(lst(aot40, tmax), "data/outputs/impact_yield_future.qs")

# tmax %>%
#   unnest_wider(data) %>%
#   select(-formap) %>%
#   unnest() %>%
#   filter(indx == '2095-2099', year %in% 2005:2009) %>%
#   ggplot(aes(x = year, y = value)) +
#   facet_wrap(~ Item, scales = 'free_y') +
#   geom_lineribbon(aes(ymax = .upper, ymin = .lower)) +
#   scale_y_percent() +
#   scale_fill_manual(values = temperature_pal[c(6,8)])
#
# p = tmax %>%
#   unnest_wider(data) %>%
#   select(-fortrend) %>%
#   arrange(factor(Item, get_order())) %>%
#   pmap(function(Item, formap){
#
#     formap = formap %>%
#       filter(indx == '2095-2099', .width == 0.66) %>%
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

# ozone no tmax ---------------------------------------------------------------------

aot40 <- inner_join(aot40_hist, aot40_fut) %>%
  mutate(aot40 = aot40_fut - aot40, .keep = "unused") %>%
  nest(.by = Item)

plan(multisession, workers = 3)

aot40 <- list(coefs_notmax, aot40, extracted_weights) %>%
  reduce(inner_join) %>%
  mutate(data = future_pmap(list(coefs, data, weights), function(acoef, adata, aweight) {
    adata <- adata %>%
      filter(year %in% 2005:2009)

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
      fgroup_by(model, indx, ADM_CODE, ADM_NAME, FAO_CODE, ISO, NAME, name, harvest_area) %>%
      fsummarise(value = mean(value)) %>%
      group_by(indx, ADM_CODE, ADM_NAME, FAO_CODE, ISO, NAME, harvest_area) %>%
      median_qi(value, .width = uncertainty)

    # weighted mean by region is needed
    fortrend <- result %>%
      fgroup_by(indx, year, name, model) %>%
      fsummarise(value = weighted.mean(value, harvest_area)) %>%
      group_by(indx, year) %>%
      median_qi(value, .width = uncertainty)

    lst(formap, fortrend)
  }, .progress = T), .keep = "unused")

plan(sequential)

# aot40 %>%
#   unnest_wider(data) %>%
#   select(-formap) %>%
#   unnest() %>%
#   ggplot(aes(x = year, y = value)) +
#   facet_wrap(~ Item, scales = 'free_y') +
#   geom_lineribbon(aes(ymax = .upper, ymin = .lower)) +
#   scale_y_percent() +
#   scale_fill_manual(values = ozone_pal[c(1,3)])

# tmax no ozone --------------------------------------------------------------------

tmax <- inner_join(tmax_hist, tmax_fut) %>%
  mutate(tmax = tmax_fut - tmax, tmax2 = tmax2_fut - tmax2, tmax3 = tmax3_fut - tmax3, .keep = "unused") %>%
  nest(.by = Item)

plan(multisession, workers = 3)

tmax <- list(coefs_noaot40, tmax, extracted_weights) %>%
  reduce(inner_join) %>%
  mutate(data = future_pmap(list(coefs, data, weights), function(acoef, adata, aweight) {
    adata <- adata %>%
      filter(year %in% 2005:2009)

    acoef <- acoef %>% select(starts_with("tmax"))

    rel_X <- adata %>%
      select(starts_with("tmax"))

    result <- tcrossprod(as.matrix(rel_X), as.matrix(acoef)) %>%
      expm1() %>%
      as_tibble() %>%
      bind_cols(adata %>% select(-starts_with("tmax"))) %>%
      pivot_longer(num_range("V", 1:1e3)) %>%
      inner_join(aweight, by = join_by(ADM_CODE, ADM_NAME, FAO_CODE))

    formap <- result %>%
      fgroup_by(model, indx, ADM_CODE, ADM_NAME, FAO_CODE, ISO, NAME, name, harvest_area) %>%
      fsummarise(value = mean(value)) %>%
      group_by(indx, ADM_CODE, ADM_NAME, FAO_CODE, ISO, NAME, harvest_area) %>%
      median_qi(value, .width = uncertainty)

    # weighted mean by region is needed
    fortrend <- result %>%
      fgroup_by(indx, year, name, model) %>%
      fsummarise(value = weighted.mean(value, harvest_area)) %>%
      group_by(indx, year) %>%
      median_qi(value, .width = uncertainty)

    lst(formap, fortrend)
  }, .progress = T), .keep = "unused")

plan(sequential)

qsave(lst(aot40, tmax), "data/outputs/impact_yield_future_no.qs")

# # tmax new ----------------------------------------------------------------
#
# tmax_hist = qread('data/outputs/tmax/gs_tmax.qs')
# tmax_fut = qread('data/outputs/tmax/gs_tmax_future.qs')
#
# season = get_season() %>%
#   filter(var == 'tmax') %>%
#   select(Item, prop)
#
# tmax_hist = inner_join(season, tmax_hist) %>%
#   mutate(tmax = map2(prop, tmax, function(aprop, atmax){
#
#     aprop = str_c('\\.', aprop)
#
#     atmax = atmax %>%
#       subset(names(atmax)[str_detect(names(atmax), aprop)])
#
#     names(atmax) <- str_remove(names(atmax), aprop)
#
#     atmax = atmax %>%
#       subset(names(atmax)[str_sub(names(atmax), start = -4) %in% 2005:2010])
#
#     return(atmax)
#
#   })) %>%
#   select(-prop)
#
# tmax = tmax_fut %>%
#   unnest() %>%
#   select(-prop) %>%
#   rename(fut = tmax) %>%
#   inner_join(tmax_hist %>% rename(hist = tmax)) %>%
#   nest(.by = Item)
#
# library(dtplyr)
# tmax = tmax %>%
#   inner_join(weights) %>%
#   inner_join(coefs) %>%
#   mutate(data = pmap(list(data, weights, coefs), function(adata, aweight, acoef){
#
#     acoef = acoef %>% select(starts_with('tmax'))
#
#     res = adata %>%
#       # head(1) %>%
#       mutate(temp = map2(fut, hist, function(afut, ahist){
#
#         stopifnot(names(ahist) == names(afut))
#
#         df = (ahist - afut) %>%
#           `names<-`(names(ahist)) %>%
#           as.data.frame(xy = T) %>%
#           pivot_longer(-c(x, y),
#                        names_to = c('name', 'year'), names_sep = '\\.', names_transform = list(year = as.integer),
#                        values_drop_na = T) %>%
#           pivot_wider() %>%
#           inner_join(aweight, by = c('x', 'y'))
#
#         rel_X = df %>%
#           select(starts_with('tmax'))
#
#         result = tcrossprod(as.matrix(rel_X), as.matrix(acoef)) %>%
#           expm1() %>%
#           as.data.table() %>%
#           lazy_dt() %>%
#           cbind(df %>% select(-starts_with('tmax'))) %>%
#           pivot_longer(num_range('V', 1:1e3)) %>%
#           as.data.table()
#
#         formap = result %>%
#           lazy_dt() %>%
#           group_by(x, y, name, weights) %>%
#           summarise(value = mean(value), .groups = 'drop') %>%
#           as_tibble()
#
#         fortrend = result %>%
#           lazy_dt() %>%
#           group_by(year, name) %>%
#           summarise(value = weighted.mean(value, weights), .groups = 'drop') %>%
#           as_tibble()
#
#         lst(formap, fortrend)
#
#       }, .progress = T), .keep = 'unused') %>%
#       unnest_wider(temp)
#
#     formap = res %>%
#       select(model, formap) %>%
#       unnest() %>%
#       lazy_dt() %>%
#       group_by(x, y) %>%
#       summarise(value = median(value), .groups = 'drop') %>%
#       as_tibble()
#
#     fortrend = res %>%
#       select(model, fortrend) %>%
#       unnest() %>%
#       group_by(year) %>%
#       median_qi(value, .width = c(0.66, 0.9))
#
#     lst(formap, fortrend)
#
#   }, .progress = T), .keep = 'unused')
#
