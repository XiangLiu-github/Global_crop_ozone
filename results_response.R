source("script/loadpackages.R")
source("script/loadfunctions.R")

quantile_df <- function(x, probs = c(0, 0.5, 1)) {
  tibble(
    val = quantile(x, probs, na.rm = TRUE),
    quant = probs
  )
}

coefs <- qread("data/outputs/boot_coefs.qs")

season <- get_season() %>%
  select(Item, var, prop)
order <- get_order()

data <- get_data() %>%
  mutate(logYield = log(Yield), aot40 = aot40 / 1e3)

data <- data %>%
  pivot_longer(c(aot40, starts_with("tmax"), starts_with("AOD"), starts_with("sm"), starts_with("cloud"))) %>%
  mutate(var = str_remove(name, "[23]")) %>%
  inner_join(season, by = join_by(Item, prop, var)) %>%
  select(-c(prop, var)) %>%
  pivot_wider() %>%
  drop_na()

datanested <- data %>%
  nest(.by = Item, .key = "fdata") %>%
  arrange(factor(Item, order))

pwalk(
  tibble(
    avar = c("aot40", "AOD", "tmax", "sm", "cloud"),
    aby = c(0.1, 0.01, 0.1, 0.1, 0.01),
    aname = list("AOT40 (ppm h)", "AOD", TeX("Daily maximum temperature ($\\degree C$)"), "Soil moisture (%)", "COD")
  ),
  function(avar, aby, aname) {
    ranges <- data %>%
      reframe(quantile_df(!!sym(avar)), .by = Item) %>%
      pivot_wider(names_from = quant, values_from = val, names_prefix = "q")

    results <- inner_join(coefs, ranges, by = join_by(Item)) %>%
      mutate(data = pmap(list(coefs, q0, q0.5, q1), function(acoefs, aq0, aq0.5, aq1) {
        acoefs <- acoefs %>%
          select(starts_with(avar))

        seqx <- seq(aq0, aq1, by = aby)

        get_X <- function(avar) {
          if (avar == "aot40") {
            (tibble(var1 = seqx) - tibble(var1 = rep(aq0.5, times = length(seqx)))) * 1e3
          } else {
            tibble(var1 = seqx, var2 = var1^2, var3 = var1^3) - tibble(var1 = rep(aq0.5, times = length(seqx)), var2 = var1^2, var3 = var1^3)
          }
        }

        X <- get_X(avar)

        tcrossprod(as.matrix(X), as.matrix(acoefs)) %>%
          rowQuantiles(probs = c(0.025, 0.2, 0.5, 0.8, 0.975)) %>%
          expm1() %>%
          as_tibble() %>%
          mutate(x = seqx)
      }), .keep = "unused")

    patch <- inner_join(datanested, results, by = join_by(Item)) %>%
      inner_join(color_tbl, by = join_by(Item)) %>%
      pmap(function(Item, data, fdata, color) {
        ggplot(aes(x = x, y = `50%`), data = data) +
          geom_hline(yintercept = 0, linetype = "longdash", color = "grey") +
          geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2, fill = color, lineend = "round", color = color) +
          geom_ribbon(aes(ymin = `20%`, ymax = `80%`), alpha = 0.6, fill = color, lineend = "round", color = color) +
          ggborderline::geom_borderline(linewidth = 1.5, lineend = "round", bordercolour = "black", color = color) +
          geom_xsidehistogram(aes(x = !!sym(avar)), bins = 50, data = fdata, inherit.aes = F, fill = color, color = "black") +
          scale_y_percent(name = NULL, expand = expansion()) +
          scale_x_continuous(name = aname, expand = expansion(mult = 0.02)) +
          scale_xsidey_continuous(labels = NULL, breaks = NULL, expand = expansion(mult = 0.02)) +
          labs(title = Item) +
          theme_half_open(font_size = 16, font_family = "Roboto Condensed") +
          background_grid(major = "x") +
          theme(
            plot.title = element_text(
              size = 25, hjust = 0.5, colour = color,
              margin = margin(0.5, 0, 0.5, 0, "lines")
            ),
            ggside.panel.scale = 0.2
          ) +
          ggside(x.pos = "bottom")
      }) %>%
      wrap_plots(ncol = 3)

    ggsave(str_c("figures/response_", avar, ".pdf"), patch, width = 1, height = 1.3, scale = 16)
  },
  .progress = T
)
