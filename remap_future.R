source("script/loadpackages.R")
source("script/loadfunctions.R")

# tmax --------------------------------------------------------------------

tmax <- tibble(
  files_ = list.files("D:/Data/CMIP6/", "^tasmax"),
  files = list.files("D:/Data/CMIP6/", "^tasmax", full.names = F)
) %>%
  separate(files_, c(NA, NA, "model", "scenario", "realization", NA, NA), "_") %>%
  filter(scenario != "historical") %>%
  mutate(model = str_c(model, "_", realization)) %>%
  reframe(files = str_c(files, collapse = " "), .by = c(model, scenario)) %>%
  mutate(cmd = str_c("cdo -P 5 -subc,273.15 -ifthen gridinfo.nc -remapbil,gridinfo.nc -mergetime [ ", files, " ] remapped/tasmax_", model, "_", scenario, ".nc"))

walk(tmax$cmd, print)


# ozone -------------------------------------------------------------------

ozone <- tibble(
  files_ = list.files("D:/Data/CMIP6/", "^o3"),
  files = list.files("D:/Data/CMIP6/", "^o3", full.names = F)
) %>%
  separate(files_, c(NA, NA, "model", "scenario", "realization", NA, "period"), "_") %>%
  mutate(model = str_c(model, "_", realization)) %>%
  reframe(files = str_c(files, collapse = " "), .by = c(model, scenario)) %>%
  mutate(cmd = str_c("cdo -P 5 -ifthen gridinfo.nc -remapbil,gridinfo.nc -sellevidx,1 -selvar,o3 -mergetime -apply,-selyear,2004/2010,2044/2050,2094/2100 [ ", files, " ] remapped/o3_", model, "_", scenario, ".nc"))

walk(ozone$cmd, print)

# scale
# str_c('cdo div o3_', unique(ozone$model), '_ssp370.nc o3_', unique(ozone$model), '_historical.nc o3_', unique(ozone$model), '_scale.nc') %>%
#   walk(print)
