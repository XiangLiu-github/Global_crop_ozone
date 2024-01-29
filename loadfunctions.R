is.Raster <- function(x) {
  return((
    class(x)[1] == "RasterLayer" ||
      class(x)[1] == "RasterBrick" ||
      class(x)[1] == "RasterStack" || class(x)[1] == "SpatRaster"
  ))
}

overlap_CN <- function(sdata) {
  CN <- st_read("https://geo.datav.aliyun.com/areas_v3/bound/100000.json", quiet = T) %>%
    pull(geometry)

  bind_rows(
    sdata,
    sdata %>%
      filter(ADM_CODE == 147295) %>%
      mutate(wkb_geometry = CN)
  )
}

get_shp <- function() {
  read_rds("data/inputs/shp/GAUL.rds") %>%
    mutate(
      ISO = countrycode::countrycode(as.integer(ADM_CODE), "gaul", "iso3c", custom_match = c("147295" = "CHN", "147296" = "TWN")),
      NAME = countrycode::countrycode(as.integer(ADM_CODE), "gaul", "country.name", custom_match = c("147295" = "China", "147296" = "Taiwan"))
    )
}

get_GTAP_region <- function() {
  read_xlsx("data/inputs/GTAP_aux/GTAP Data Base 10 Regions.xlsx") %>%
    fill(GTAP_CODE, GTAP_NAME) %>%
    inner_join(read_xlsx("data/inputs/GTAP_aux/GTAP_modeling_regions.xlsx") %>%
      mutate(GTAP_CODE = str_to_upper(GTAP_CODE))) %>%
    mutate(NAME = countrycode::countrycode(MEMBER_CODE, "iso3c", "country.name"))
}

get_FAO_groups <- function() {
  temp <- tibble(
    GROUP_NAME = c(
      "Eastern Africa", "Middle Africa", "Northern Africa", "Southern Africa",
      "Western Africa", "Northern America", "Central America", "Caribbean",
      "South America", "Central Asia", "Eastern Asia", "Southern Asia",
      "South-Eastern Asia", "Western Asia", "Eastern Europe", "Northern Europe",
      "Southern Europe", "Western Europe", "Oceania"
    ),
    GROUP_CODE = str_c("FAO", 1:length(GROUP_NAME))
  )

  read_xls("data/inputs/shp/G2015_GroupingsM49.xls") %>%
    select(-FAOSTAT_GROUP_CODE) %>%
    filter(GROUP_NAME %in% temp$GROUP_NAME) %>%
    separate_rows(ADM_CODE, sep = " \\+ ") %>%
    mutate(
      across(everything(), as.character),
      TERRITORY_NAME = case_when(
        ADM_CODE == 147295 ~ "China, mainland",
        ADM_CODE == 147296 ~ "China, Taiwan Province of",
        TRUE ~ TERRITORY_NAME
      ),
      FAO_CODE = case_when(
        ADM_CODE == 147295 ~ "41",
        ADM_CODE == 147296 ~ "214",
        TRUE ~ FAO_CODE
      ),
      ISO = countrycode::countrycode(ADM_CODE, "gaul", "iso3c", custom_match = c("147295" = "CHN", "147296" = "TWN")),
      NAME = countrycode::countrycode(ADM_CODE, "gaul", "country.name", custom_match = c("147295" = "China", "147296" = "Taiwan"))
    ) %>%
    rename(ADM_NAME = TERRITORY_NAME) %>%
    inner_join(temp)
}

uncertainty = c(0.66, 0.95)

color_tbl <- tibble(
  Item = c("Sugarbeet", "Maize", "Soybean", "Sunflower", "Barley", "Rye", "Rapeseed", "Wheat", "Sorghum", "Pea", "Millet", "Potato", "Bean", "Groundnut", "Rice", "Cotton", "Cassava", "Sugarcane"),
  color = c("#F2C4E2", "#FFD700", "#8B4513", "#FFA500", "#FFA07A", "#CD853F", "#FF69B4", "#0420BF", "#8B0000", "#00FF7F", "#8B008B", "#D2691E", "#A0A603", "#BD4D24", "#253C59", "#1E90FF", "#D2B48C", "#228B22")
)

yield_pal <- NatParksPalettes::natparks.pals("Olympic", n = 10, direction = -1)
ozone_pal <- NatParksPalettes::natparks.pals("Arches2", n = 10, direction = -1)
aerosol_pal <- rev(PNWColors::pnw_palette("Mushroom", n = 10))
temperature_pal <- PNWColors::pnw_palette("Bay", n = 10)
moisture_pal <- rev(PNWColors::pnw_palette("Anemone", n = 10))
cloud_pal <- NatParksPalettes::natparks.pals("Glacier", n = 10, direction = -1)

ozone_trend_pal <- MetBrewer::met.brewer("Cassatt2", n = 10, direction = -1)
aerosol_trend_pal <- NatParksPalettes::natparks.pals("Acadia", n = 10)
temperature_trend_pal <- MetBrewer::met.brewer("Benedictus", n = 10, direction = -1)

detrend <- function(value, year) {
  lm(value ~ poly(year, 2), data = NULL)$residuals
}

round_pettry <- function(value) sprintf("%0.2f", round(value, digits = 2))

get_data <- function() {
  read_rds("data/tidied.rds")
  # mutate(across(c(aot40, starts_with('AOD'), starts_with('tmean'), starts_with('sm'), starts_with('cloud')), DescTools::Winsorize, probs = c(0.01, 0.99)), .by = Item)
}

get_season <- function() {
  cvdata <- qread("data/outputs/oob_seasons.qs")

  summed <- cvdata %>%
    unnest() %>%
    pivot_longer(-c(Item, cors), names_to = "var", values_to = "prop") %>%
    slice_max(cors, n = 1, by = c(Item, var))

  return(summed)
}

get_order <- function() {
  cvdata <- qread("data/outputs/oob.qs")

  summed <- cvdata %>%
    unnest() %>%
    group_by(Item) %>%
    median_qi(Testing, .width = uncertainty[2]) %>%
    group_by(.lower < 0) %>%
    arrange(desc(Testing), .by_group = T) %>%
    ungroup() %>%
    pull(Item)

  return(summed)
}

object_size <- function(aobject) {
  aobject %>%
    object.size() %>%
    print(unit = "auto")
}

if (Sys.info()[1] == "Windows") {
  qn <- 10
} else {
  qn <- 3
}

cal_significance <- function(model1, model2, type) {
  # model1 is full model
  # model2 is nested model

  if (type == "tmax") {
    coef_model1 <- coef(model1)[c("tmax", "tmax2", "tmax3")]
    cov_model1 <- vcov(model1)[c("tmax", "tmax2", "tmax3"), c("tmax", "tmax2", "tmax3")]
    coef_model2 <- coef(model2)[c("tmax", "tmax2", "tmax3")]
    cov_model2 <- vcov(model2)[c("tmax", "tmax2", "tmax3"), c("tmax", "tmax2", "tmax3")]

    test_statistic <- t(coef_model1 - coef_model2) %*% solve(cov_model1 + cov_model2) %*% (coef_model1 - coef_model2)

    df1 <- nrow(cov_model1)
    df2 <- nrow(cov_model2)
    df <- df1 + df2
    p_value <- 2 * pt(abs(test_statistic), df, lower.tail = FALSE)

    return(as.vector(p_value))
  } else if (type == "ozone") {
    coef_model1 <- coef(model1)["aot40"]
    cov_model1 <- vcov(model1)["aot40", "aot40"]
    coef_model2 <- coef(model2)["aot40"]
    cov_model2 <- vcov(model2)["aot40", "aot40"]

    test_statistic <- (coef_model1 - coef_model2) / sqrt(cov_model1 + cov_model2)
    p_value <- 2 * pt(abs(test_statistic), df = Inf, lower.tail = FALSE)

    return(p_value)
  }
}

arrow <-
  arrow(
    length = unit(0.015, "npc"),
    ends = "last",
    type = "open"
  )

relpred <- function(object, newdata, baseline = NULL, level = 0.90, x = NULL) {
  # object: lm or fit_feols fitted model
  # newdata: data.frame (or matrix) of X variables with names that correspond to fitted model (names that don't fit are ignored)
  # baseline: list with same length as the number of columns in newdata (or, eventually, matching names). without this, baseline is assumed to be zero.
  # object <- lm(y ~ x + I(x^2) + group, data = data)
  # object <- feols(y ~ x + I(x^2) | group, data = data)
  # newdata <- newdata; level = 0.95
  # baseline <- NULL
  # baseline <- list(x = 1)
  # END DEBUG

  # TODO: Add support for factor variables in newdata later. Requirements
  # - Must either check or convert variables to match levels in object.
  # - Decide how baseline should be chosen for factor variables if baseline is explicitly supplied.
  if (any(sapply(newdata, function(x) !is.numeric(x)))) {
    stop("For now, only numeric variables are supported in newdata.")
  }

  if (!is.null(baseline) & ncol(newdata) != length(baseline)) {
    stop("baseline vector length must be equal to number of columns in newdata.")
  }

  # Populate a data.frame with all variables used in object estimation and return a design matrix that matches that object.
  newdata_filled <- fill_missing_vars(object, newdata)

  # If baseline isn't passed, assume zero
  if (is.null(baseline)) {
    baseline_df <- newdata
    baseline_df[] <- 0
  } else {
    baseline_df <- as.data.frame(baseline)
  }

  baseline_filled <- fill_missing_vars(object, baseline_df)

  # Subtract baseline from newdata to get X
  X <- as.matrix(newdata_filled - baseline_filled)

  # Compute coefficients and standard errors ----
  B <- as.numeric(coef(object))

  # Compute degrees of freedom
  if (inherits(object, "fixest")) {
    # df <- attributes(vcov(object, attr = T))$G
    df <- degrees_freedom(object, type = "t")
  } else {
    df <- object$df.residual
  }

  # Drop any variables in X that are not the fitted model object
  # This will handle a variety of issues, including a missing intercept term
  X <- X[, names(coef(object))]


  fit <- data.frame(fit = X %*% B)
  sig <- vcov(object)
  se <- apply(X, MARGIN = 1, FUN = get_se, sig = sig)

  t_val <- qt((1 - level) / 2 + level, df = df)
  fit$lwr <- fit$fit - t_val * se
  fit$upr <- fit$fit + t_val * se

  if (is.null(x)) {
    fit
  } else {
    fit %>% mutate(x = x)
  }
}


fill_missing_vars <- function(object, X) {
  # Populate a data.frame with all variables used in object estimation and return a design matrix with the same coefficients.
  # Factor variables filled in using the first factor level, numeric filled with zeroes

  orig_vars <- all.vars(formula(object))[-1]

  # Add any factor variables that are in the RHS but not newdata with the baseline value
  # TODO: How to add back a missing factor variable for feols? We can do it with lm using object$xlevels
  # TODO: If the user passes, e.g., factor(group) this will not work. Not sure how to resolve.
  fact_vars <- names(object$xlevels)
  for (f in fact_vars) {
    if (!(f %in% names(X))) {
      X[[f]] <- factor(object$xlevels[[f]][1],
        levels = object$xlevels[[f]]
      )
    } else {
      # This is risky, but will ensure levels in object match those in X...
      X[[f]] <- factor(X[[f]],
        levels = object$xlevels[[f]]
      )
    }
  }

  # If there are any other variables in the formula but not in X, add them with value = 0
  for (v in orig_vars) {
    if (!(v %in% names(X))) {
      X[[v]] <- 0
    }
  }

  # Extract terms object, removing response variable
  tt <- terms(object)
  Terms <- delete.response(tt)

  as.data.frame(model.matrix(Terms, data = X))
}

get_se <- function(r, sig) {
  # Compute linear combination, helper function for predict_partial
  # Given numeric vector r (the constants) and vcov sig, compute SE
  r <- matrix(r, nrow = 1)
  sqrt(r %*% sig %*% t(r))
}
