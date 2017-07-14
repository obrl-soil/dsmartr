# modify dsmart demo data for dsmartr

library(sp)
library(sf)
library(raster)
library(tidyverse)

load(file.path(getwd(), 'data', 'dsT_composition.rda'))
load(file.path(getwd(), 'data', 'dsT_polygons.rda'))

dst_comp_wide <- dsT_composition %>%
  select(-mapunit) %>%
  rename(NUID = poly, CLASS = soil_class, PERC= proportion) %>%
  mutate(CLASS = as.character(CLASS)) %>%
  gather(key, value, CLASS, PERC) %>%
  group_by(NUID, key) %>%
  mutate(nn = row_number()) %>%
  ungroup() %>%
  unite(concat, key, nn) %>%
  spread(concat, value)

dst_sf <- st_as_sf(dsT_polygons) %>%
  st_transform(., 3577) %>%
  select(POLY_NO, MAP_CODE) %>%
  rename(NUID = POLY_NO) %>%
  mutate(MAP_CODE = as.character(MAP_CODE)) %>%
  left_join(., dst_comp_wide, by = 'NUID')

save(dst_sf, file= file.path(getwd(), 'data', 'dst_sf_polygons.rda'))
