source("script/loadpackages.R")
source("script/loadfunctions.R")

temp <- raster("data/inputs/crop_calendar/mai_rf_ggcmi_crop_calendar_phase3_v1.01.nc4")

data <- tibble(
  Item = c(
    "Barley", "Bean", "Cassava", "Cotton", "Groundnut", "Maize",
    "Millet", "Pea", "Potato", "Rye", "Rapeseed", "Sorghum", "Soybean",
    "Sugarbeet", "Sugarcane", "Sunflower", "Wheat", "Rice"
  ),
  weights = map(Item, function(aitem) {
    aitem <- str_to_lower(aitem)
    raster(str_c("/vsizip/data/inputs/cropland/Monfreda/", aitem, "_HarvAreaYield_Geotiff.zip/", aitem, "_HarvAreaYield_Geotiff/", aitem, "_HarvestedAreaHectares.tif")) %>%
      exact_resample(temp, "sum")
  })
)

qsave(data, "data/inputs/cropland/weights.qs")
