source("script/loadpackages.R")
source("script/loadfunctions.R")

season <- get_season() %>%
  select(Item, var, prop)
order <- get_order()

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

cvdata <- data %>%
  mutate(cor = map(fdata, function(afdata) {
    defdata <- afdata %>%
      select(logYield, starts_with("aot40"), starts_with("tmax"), starts_with("AOD"), starts_with("sm"), starts_with("cloud")) %>%
      names() %>%
      str_c(collapse = ", ") %>%
      str_c("c(", ., ") ~ 1 | ADM_CODE[year + year ^ 2]") %>%
      as.formula() %>%
      feols(data = afdata) %>%
      map_dfc("residuals") %>%
      purrr::set_names(~ str_remove(.x, "lhs: ")) %>%
      mutate(ADM_CODE = afdata$ADM_CODE)

    set.seed(2023)

    bootstraps <- group_vfold_cv(defdata, ADM_CODE, v = 10)

    map_dfr(bootstraps$splits, function(asplit) {
      aanalysis <- analysis(asplit)
      aassessment <- assessment(asplit)

      list(
        fml_full = logYield ~ aot40 + tmax + tmax2 + tmax3 + AOD + AOD2 + AOD3 + sm + sm2 + sm3 + cloud + cloud2 + cloud3,
        fml_no_ozone = logYield ~ tmax + tmax2 + tmax3 + AOD + AOD2 + AOD3 + sm + sm2 + sm3 + cloud + cloud2 + cloud3,
        fml_no_tmax = logYield ~ aot40 + AOD + AOD2 + AOD3 + sm + sm2 + sm3 + cloud + cloud2 + cloud3
        # fml_no_sm = logYield ~ aot40 + tmax + tmax2 + tmax3 + AOD + AOD2 + AOD3 + cloud + cloud2 + cloud3,
        # fml_no_opacity = logYield ~ aot40 + tmax + tmax2 + tmax3 + sm + sm2 + sm3
      ) %>%
        map(feols, data = aanalysis) %>%
        map_dfc(predict, newdata = aassessment) %>%
        mutate(obs = aassessment$logYield)
    }) %>%
      pivot_longer(-obs) %>%
      summarise(
        cor = cor(obs, value),
        rmse = rmse_vec(obs, value), .by = name
      )
  }), .keep = "unused") %>%
  unnest()

cvdata %>%
  ggplot(aes(x = name, y = cor)) +
  facet_wrap(~ factor(Item, order)) +
  geom_col()

# this should compare with a model that only has fixed effects
cvdata %>%
  select(-cor) %>%
  pivot_wider(values_from = rmse) %>%
  # mutate(across(starts_with('fml_no'), ~ (fml_full - .x) / fml_full)) %>%
  # select(-fml_full) %>%
  pivot_longer(-Item) %>%
  ggplot(aes(x = factor(Item, order), y = value, fill = name)) +
  geom_col(position = position_dodge())

saved <- qread("data/outputs/oob.qs") %>%
  unnest() %>%
  group_by(Item) %>%
  median_qi(Testing, .width = uncertainty[2]) %>%
  group_by(.lower < 0) %>%
  arrange(desc(Testing), .by_group = T) %>%
  filter(`.lower < 0` == F) %>%
  pull(Item)

p1 <- cvdata %>%
  filter(Item %in% saved) %>% 
  select(-rmse) %>%
  pivot_wider(values_from = cor) %>%
  mutate(across(starts_with("fml_no"), ~ (fml_full - .x) / fml_full)) %>%
  select(-fml_full) %>%
  pivot_longer(-Item) %>%
  ggplot(aes(x = factor(Item, order), y = value, fill = factor(name, labels = c("Full model without ozone", "Full model without temperature")))) +
  geom_hline(yintercept = 0, linetype = "longdash", color = "grey") +
  geom_col(position = position_dodge(), width = 0.7, linewidth = 0.5, color = "black") +
  geom_text_repel(aes(x = factor(Item, order), y = -0.05, label = str_c(round_pettry(delta * 100), "%")), data = . %>% pivot_wider() %>% mutate(delta = fml_no_tmax - fml_no_ozone), inherit.aes = F, size = 6, family = "Roboto Condensed") +
  scale_fill_manual(values = c("#008500", "#FF8500")) +
  # ggbreak::scale_y_break(c(0.3, 0.8)) +
  scale_y_percent("Relative change in within-r") +
  theme_half_open(font_size = 16, font_family = "Roboto Condensed") +
  background_grid(major = "x") +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.justification = 0.5,
    axis.line.y.right = element_blank(),
    axis.title.y.right = element_blank(),
    axis.text.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.title.x = element_blank()
  )

ggsave("figures/cv_ozone_temperature_term.pdf", p1, width = 5, height = 2, scale = 3)

# pdftools::pdf_subset("figures/cv_ozone_temperature_term.pdf", 2, "figures/temp.pdf")

# file.rename("figures/temp.pdf", "figures/cv_ozone_temperature_term.pdf")
