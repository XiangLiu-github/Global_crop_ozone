source("script/loadPackages.R")
source("script/loadFunctions.R")

data <- read_rds('data/inputs/US_ozone.rds') %>%
  mutate(O3 = rollingo3)

model_aot40 <- bam(AOT40 ~ s(O3), data = data, family = tw(), discrete = T)

folds <- vfold_cv(data)

plan(multisession, workers = 5)

cv_aot40 <- future_map_dfr(folds$splits, function(asplit) {
  model <- bam(AOT40 ~ s(O3),
               data = asplit %>% analysis(),
               family = tw(),
               discrete = T
  )
  
  tibble(
    predicted = exp(predict(model, asplit %>% assessment())),
    raw = asplit %>% assessment() %>% pull(AOT40)
  )
}, .progress = T)

plan(sequential)

p1 <- data %>%
  ggplot(aes(x = O3, y = AOT40)) +
  geom_pointdensity(adjust = 4, method = "default") +
  scale_color_viridis(
    trans = "log10",
    name = "n_neighbors",
    limits = c(1, 8e2), oob = scales::squish,
    option = "magma"
  ) +
  xlab(TeX("MDA8 (ppb)")) +
  ylab(TeX("AOT40 (ppb h)")) +
  scale_y_continuous(limits = c(0, 23e3)) +
  scale_x_continuous(limits = c(0, 150)) +
  geom_xsidehistogram(bins = 60) +
  geom_ysidehistogram(bins = 60) +
  ggside(x.pos = "top", y.pos = "right", collapse = "all") +
  scale_xsidey_continuous(labels = NULL) +
  scale_ysidex_continuous(labels = NULL) +
  geom_line(aes(x = O3, y = AOT40_pred),
            inherit.aes = F,
            data = tibble(O3 = seq(0, 150, by = 0.1)) %>%
              mutate(AOT40_pred = exp(predict(model_aot40, .))),
            size = 1
  ) +
  theme_half_open(18, font_family = "Roboto Condensed") +
  background_grid() +
  theme(legend.position = c(0.2, 0.6))

p2 <- cv_aot40 %>%
  ggplot(aes(x = raw, y = predicted)) +
  geom_pointdensity(adjust = 4, method = "default") +
  scale_color_viridis(
    trans = "log10",
    name = "n_neighbors",
    limits = c(1, 1e2), oob = scales::squish,
    option = "magma"
  ) +
  geom_abline(slope = 1, intercept = 0, size = 1, color = "blue") +
  geom_smooth(method = "lm", se = F, color = "red", size = 1) +
  xlab(TeX("Observation AOT40 (ppb h)")) +
  ylab(TeX("Prediction AOT40 (ppb h)")) +
  stat_poly_eq(
    aes(label = paste(after_stat(rr.label),
                      after_stat(n.label),
                      sep = "*\", \"*"
    )),
    method = "lm", family = "Roboto Condensed", size = 5
  ) +
  scale_x_continuous(limits = c(0, 23e3)) +
  scale_y_continuous(limits = c(0, 23e3)) +
  geom_xsidehistogram(bins = 60) +
  geom_ysidehistogram(bins = 60) +
  ggside(x.pos = "top", y.pos = "right", collapse = "all") +
  scale_xsidey_continuous(labels = NULL, limits = c(0, 2e4)) +
  scale_ysidex_continuous(labels = NULL, limits = c(0, 2e4)) +
  theme_half_open(18, font_family = "Roboto Condensed") +
  background_grid() +
  theme(legend.position = c(0.6, 0.2))

patch <-
  (p1 + p2) +
  # (p3 + p4) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 30, family = "Roboto Condensed"))

ggsave("figures/gam_performance.pdf", patch, width = 2, height = 1, scale = 8)
