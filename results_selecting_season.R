source("script/loadpackages.R")
source("script/loadfunctions.R")

data <- get_data() %>%
  mutate(logYield = log(Yield)) %>%
  nest(.by = c(Item), .key = "fdata")

plan(multisession, workers = 5)

cvdata <- data %>%
  mutate(fdata = future_map(fdata, function(afdata) {
    afdata <- afdata %>%
      pivot_wider(
        names_from = prop,
        values_from = c(aot40, starts_with("tmax"), starts_with("AOD"), starts_with("sm"), starts_with("cloud"))
      ) %>%
      drop_na() # this will drop some observations

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

    expand_grid(
      var = list("aot40", c("tmax", "tmax2", "tmax3"), c("AOD", "AOD2", "AOD3"), c("sm", "sm2", "sm3"), c("cloud", "cloud2", "cloud3")),
      prop = c("0_2", "0_4", "0_6", "0_8", "0_10", "2_4", "2_6", "2_8", "2_10", "4_6", "4_8", "4_10", "6_8", "6_10", "8_10")
    ) %>%
      mutate(cors = map2(var, prop, function(avar, aprop) {
        avar <- str_c(avar, "_", aprop)

        avar <- str_c(avar, collapse = "+")

        map_dfr(bootstraps$splits, function(asplit) {
          analysisd <- analysis(asplit)
          assessmentd <- assessment(asplit)

          model <- model <- feols(as.formula(str_c("logYield ~ ", avar)), data = analysisd)

          # tibble(
          #   Training = cor(analysisd$logYield, predict(model, analysisd)),
          #   Testing = cor(assessmentd$logYield, predict(model, assessmentd))
          # )

          tibble(obs = assessmentd$logYield, prd = predict(model, assessmentd))
        }, .id = "id") %>%
          summarise(cor = cor(obs, prd))
      })) %>%
      mutate(var = map_chr(var, ~ unique(str_remove(.x, "[23]")))) %>%
      unnest(cols = c(cors))
  }, .progress = T, .options = furrr_options(seed = T)))

plan(sequential)

qsave(cvdata, "data/outputs/oob_seasons.qs")

cvdata <- qread("data/outputs/oob_seasons.qs")

cvdata %>%
  unnest() %>%
  # group_by(Item, var, prop) %>%
  # median_qi(Testing, .width = 0.9) %>%
  slice_max(cor, n = 1, by = c(Item, var))


p <- cvdata %>%
  unnest() %>%
  ggplot(aes(x = factor(var, labels = c("AOD", "AOT40", "COD", "W", "T")), y = factor(prop, c("0_2", "0_4", "0_6", "0_8", "0_10", "2_4", "2_6", "2_8", "2_10", "4_6", "4_8", "4_10", "6_8", "6_10", "8_10")))) +
  facet_wrap(~ factor(Item, get_order()), ncol = 3) +
  geom_tile(aes(fill = cor), color = "grey90") +
  geom_text(aes(label = round_pettry(cor)), data = . %>% slice_max(cor, n = 1, by = c(Item, var)), family = "Roboto Condensed") +
  scale_x_discrete("Variable", expand = expansion()) +
  scale_y_discrete("Season candidate", expand = expansion()) +
  scale_fill_gradientn(
    name = "Within r",
    colours = NatParksPalettes::natparks.pals("Arches", n = 11),
    limits = symmetric_limits,
    breaks = breaks_extended(n = 5),
    guide = guide_coloursteps(title.position = "top")
  ) +
  theme_half_open(font_size = 16, font_family = "Roboto Condensed") +
  theme(
    strip.background = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.width = unit(7, "lines"),
    legend.key.height = unit(0.5, "lines"),
    legend.justification = 0.5,
    legend.title.align = 0.5
  )

ggsave("figures/selected_seasons.pdf", p, width = 1, height = 2, scale = 10)
