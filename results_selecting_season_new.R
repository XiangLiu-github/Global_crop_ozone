source("script/loadpackages.R")
source("script/loadfunctions.R")

data <- get_data() %>%
  mutate(logYield = log(Yield)) %>%
  nest(.by = c(Item), .key = "fdata")

props <- c("0_2", "0_4", "0_6", "0_8", "0_10", "2_4", "2_6", "2_8", "2_10", "4_6", "4_8", "4_10", "6_8", "6_10", "8_10")

if (F) {
  plan(multisession, workers = 10)

  fmls <- expand_grid(aot40 = props, tmax = props, AOD = props, sm = props, cloud = props) %>%
    mutate(fml = future_pmap(pick(1:5), function(aot40, tmax, AOD, sm, cloud) {
      x1 <- str_c("aot40.", aot40)
      x2 <- map_chr(c("", "2", "3"), ~ str_c("tmax", .x, ".", tmax))
      x3 <- map_chr(c("", "2", "3"), ~ str_c("AOD", .x, ".", AOD))
      x4 <- map_chr(c("", "2", "3"), ~ str_c("sm", .x, ".", sm))
      x5 <- map_chr(c("", "2", "3"), ~ str_c("cloud", .x, ".", cloud))

      vars <- str_c(c(x1, x2, x3, x4, x5), collapse = "+")

      as.formula(str_c("logYield ~ ", vars))
    }, .progress = T))

  plan(sequential)

  qsave(fmls, "data/outputs/fmls.qs")
} else {
  fmls <- qread("data/outputs/fmls.qs")
}

fmls <- fmls %>% filter(AOD == "0_10", sm == "0_10", cloud == "0_10")

# fmls = expand_grid(aot40 = props, tmax = props, AOD = props, sm = props, cloud = props) %>%
#   mutate(x1 = stringi::stri_join('aot40.', aot40),
#          x20 = stringi::stri_join('tmax.', tmax),
#          x22 = stringi::stri_join('tmax2.', tmax),
#          x23 = stringi::stri_join('tmax3.', tmax),
#          x30 = stringi::stri_join('AOD.', AOD),
#          x32 = stringi::stri_join('AOD2.', AOD),
#          x33 = stringi::stri_join('AOD3.', AOD),
#          x40 = stringi::stri_join('sm.', sm),
#          x42 = stringi::stri_join('sm2.', sm),
#          x43 = stringi::stri_join('sm3.', sm),
#          x50 = stringi::stri_join('cloud.', cloud),
#          x52 = stringi::stri_join('cloud2.', cloud),
#          x53 = stringi::stri_join('cloud3.', cloud),
#          fml = stringi::stri_join('logYield ~ ', stringi::stri_join(c(x1, x20, x22, x23, x30, x32, x33, x40, x42, x43, x50, x52, x53), collapse = '+'))) %>%
#   select(-starts_with('x'))

plan(multisession, workers = 5)

cvdata <- data %>%
  mutate(fdata = future_map(fdata, function(afdata) {
    afdata <- afdata %>%
      pivot_wider(
        names_from = prop, names_sep = ".",
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

    fmls %>%
      mutate(cors = map_dbl(fml, function(afml) {
        map_dfr(bootstraps$splits, function(asplit) {
          analysisd <- analysis(asplit)
          assessmentd <- assessment(asplit)

          model <- model <- feols(afml, data = analysisd)

          tibble(obs = assessmentd$logYield, prd = predict(model, assessmentd))
        }, .id = "id") %>%
          summarise(cor = cor(obs, prd)) %>%
          deframe()
      }, .progress = F), .keep = "unused")
  }, .progress = T, .options = furrr_options(seed = T)))

plan(sequential)

qsave(cvdata, "data/outputs/oob_seasons.qs")

cvdata <- qread("data/outputs/oob_seasons.qs")

cvdata %>%
  unnest() %>%
  # group_by(Item, var, prop) %>%
  # median_qi(Testing, .width = 0.9) %>%
  slice_max(cors, n = 1, by = c(Item))

p <- cvdata %>%
  unnest() %>%
  mutate(rank = rank(desc(cors)), .by = Item) %>%
  # filter(Item %in% get_order()[1:6]) %>%
  ggplot(aes(
    x = factor(aot40, 
               c("0_2", "0_4", "0_6", "0_8", "0_10", "2_4", "2_6", "2_8", "2_10", "4_6", "4_8", "4_10", "6_8", "6_10", "8_10"), 
               c("00-02", "00-04", "00-06", "00-08", "00-10", "02-04", "02-06", "02-08", "02-10", "04-06", "04-08", "04-10", "06-08", "06-10", "08-10")),
    y = factor(tmax, 
               c("0_2", "0_4", "0_6", "0_8", "0_10", "2_4", "2_6", "2_8", "2_10", "4_6", "4_8", "4_10", "6_8", "6_10", "8_10"),
               c("00-02", "00-04", "00-06", "00-08", "00-10", "02-04", "02-06", "02-08", "02-10", "04-06", "04-08", "04-10", "06-08", "06-10", "08-10"))
  )) +
  facet_wrap(~ factor(Item, get_order()), ncol = 3) +
  geom_tile(aes(fill = rank), color = "grey30", linewidth = 0.5) +
  ggstar::geom_star(data = . %>% slice_max(cors, n = 1, by = c(Item)), size = 2, fill = 'white', color = NA, starshape = 1) +
  geom_text_repel(aes(label = round_pettry(cors)), data = . %>% slice_max(cors, n = 1, by = c(Item)), family = "Roboto Condensed", size = 5, color = "white") +
  scale_x_discrete("AOT40 season candidates", expand = expansion()) +
  scale_y_discrete("Temperature season candidates", expand = expansion()) +
  scale_fill_gradientn(
    name = "Rank of within-r",
    colours = NatParksPalettes::natparks.pals("Arches", n = 10, direction = -1)[1:5],
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
    legend.title.align = 0.5,
    axis.text.x.bottom = element_text(angle = 90, vjust = 0.5)
  )

# ggsave("selected_seasons_mini.pdf", p, width = 3, height = 2.5, scale = 3)
# library(magick)
# image_read_pdf('selected_seasons_mini.pdf', density = 100) %>%
  # image_write(path = 'selected_seasons_mini.png', format = "png")
ggsave("figures/selected_seasons.pdf", p, width = 1, height = 2, scale = 10)
