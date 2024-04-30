

geoconflict_main_weights <- geoconflict_main %>%
  filter(year == 2011) %>%
  st_as_sf(coords = c("lat", "lon"), crs = "WGS84")

#geoconflict_main_renmae <- geoconflict_main %>%
#  st_as_sf(coords = c("lat", "lon"), crs = "WGS84")
  
neighboors <- dnearneigh(st_centroid(geoconflict_main_weights), d1 = 0, d2 = 180)
#neighboors <- poly2nb(geoconflict_main_weights, queen = TRUE)

listwmat <- nb2listw(neighboors, style ="W")

spweightmat <- listw2mat(listwmat)

summary(spweightmat)


fm1 <- ANY_EVENT_ACLED ~ SPEI4pg + GSmain_ext_SPEI4pg

sararbfemod <- spml(fm1, 
                    geoconflict_main, 
                    listw=listwmat, 
                    lag=TRUE, 
                    spatial.error = "none",
                    model="within", 
                    effect="twoways") 

summary(sararbfemod)

res1  <- blmpSDPD(fm1, data=geoconflict_main, 
                  listwmat,
                  index = c("cell","year"),
                  dynamic = TRUE,
                  model = list("sar","sdm","sem","sdem"),
                  effect = "twoways")



# Convert to an sf object
data_sf <- st_as_sf(geoconflict_main, coords = c("lon", "lat"), crs = 4326)  # EPSG:4326 is WGS 84

# Remove duplicates
data_sf <- data_sf[!duplicated(st_coordinates(data_sf)), ]

# Create a distance-based spatial weights matrix
# Convert distance to class `nb` for neig4hbor definition
dist_threshold <- 180  # Convert km to meters if CRS is in meters
nb <- dnearneigh(st_coordinates(data_sf), 0, dist_threshold)

# Convert to a binary spatial weights matrix
w_bin_180 <- nb2listw(nb, style = "B", zero.policy = TRUE)






library(plm)
library(sphet)
library(xtable)
library(lmtest)


# Save or use the spatial weights matrix
#saveRDS(w_bin_180, "W_bin_180.rds")
w_bin_180 <- readRDS("W_bin_180.rds")


# Regression Model I: Fixed Effects with Linear Time Trend
model_fe <- plm(ANY_EVENT_ACLED ~ SPEI4pg + L1_SPEI4pg + L2_SPEI4pg + GSmain_ext_SPEI4pg + L1_GSmain_ext_SPEI4pg + L2_GSmain_ext_SPEI4pg +  elevation_cell + rough_cell + area_cell + use_primary + 
                  dis_river_cell + shared + border + any_mineral + ELF + as.factor(country_largest_share):as.numeric(year),
                data = geoconflict_main, model = "within", index = c("cell", "year"),inst.method = c("am"), effect = "time")
summary(model_fe, vcovHC(model_fe, method = "arellano"))

# Regression Model II: Spatial HAC
model_hac <- sphet(ANY_EVENT_ACLED ~ SPEI4pg + L1_SPEI4pg + L2_SPEI4pg + GSmain_ext_SPEI4pg + L1_GSmain_ext_SPEI4pg + L2_GSmain_ext_SPEI4pg +  elevation_cell + rough_cell + area_cell + use_primary + 
                  dis_river_cell + shared + border + any_mineral + ELF + as.factor(country_largest_share):as.numeric(year),
                   data = data_sf, listw = w_bin_180, model = "HAC", doTest = FALSE)
summary(model_hac)


spreg(ANY_EVENT_ACLED ~ SPEI4pg + L1_SPEI4pg + L2_SPEI4pg + GSmain_ext_SPEI4pg + L1_GSmain_ext_SPEI4pg + L2_GSmain_ext_SPEI4pg +  elevation_cell + rough_cell + area_cell + use_primary + 
        dis_river_cell + shared + border + any_mineral + ELF,
      data = data_sf, listw = w_bin_180, listw2 = NULL, 
      endog = NULL, instruments = NULL, 
      lag.instr = FALSE, initial.value = 0.2, 
      model = c("sarar"),
      het = FALSE, verbose = FALSE, 
      na.action = na.fail,  HAC = FALSE, 
      distance = NULL, type =  c("Epanechnikov"), 
      bandwidth = "variable", Durbin = TRUE)




# Regression Model III: Spatial Durbin Model (SDM) with Fixed Effects
model_sdm <- spml(ANY_EVENT_ACLED ~ SPEI4pg + L1_SPEI4pg + L2_SPEI4pg + GSmain_ext_SPEI4pg + L1_GSmain_ext_SPEI4pg + L2_GSmain_ext_SPEI4pg +  elevation_cell + rough_cell + area_cell + use_primary + 
                    dis_river_cell + shared + border + any_mineral + ELF, 
                  data = data_sf, listw = w_bin_180, model = "within", spatial.model = "none")
summary(model_sdm)


sararbfemod <- spml(fm3, Produc, listw=usalw, lag=TRUE, spatial.error = "b", 
                    model="within", effect="twoways") 
summary(sararbfemod)

# Output results in a readable format
print(xtable(summary(model_fe, vcovHC(model_fe, method = "arellano"))), type = "html")
print(xtable(summary(model_hac)), type = "html")
print(xtable(summary(model_sdm)), type = "html")














pacman::p_load(spatialreg, fixest, splm, stringi, stringr, stringdist, haven, sf, dplyr, fuzzyjoin, 
               comparator, digest, zoomerjoin, ggplot2, tidyr, ggthemes, viridis, 
               fixest, conleyreg, plm, stargazer, magrittr, tidyverse, tmap, spdep, SDPDmod,
               igraph, generics, knitr, kableExtra, formatR,readxl, haven, units)



setwd("~/Desktop/GITHUB/spatial_econ/data/03_assignment/dataset")

raster_Africa <- read_sf("~/Desktop/GITHUB/spatial_econ/data/03_assignment/dataset/raster_Africa.shp")
geoconflict_main <- read_dta("~/Desktop/GITHUB/spatial_econ/data/03_assignment/dataset/geoconflict_main.dta")
intersect_coord <- read_dta("~/Desktop/GITHUB/spatial_econ/data/03_assignment/dataset/intersect_coord.dta")



library(fixest)

geoconflict_main_if <- geoconflict_main %>% 
  filter(!year == 1997)

geoconflict_main_if$year <- as.factor(geoconflict_main_if$year)
geoconflict_main_if$cell <- as.factor(geoconflict_main_if$cell)

mod1 <- feols(ANY_EVENT_ACLED ~ SPEI4pg + L1_SPEI4pg + L2_SPEI4pg + GSmain_ext_SPEI4pg + L1_GSmain_ext_SPEI4pg + L2_GSmain_ext_SPEI4pg + elevation_cell + rough_cell + area_cell + as.factor(use_primary) + 
        dis_river_cell + as.factor(shared) +  as.factor(border) + as.factor(any_mineral) + (ELF) + i(country_largest_share, as.numeric(year)) | as.factor(year), 
      data = geoconflict_main_if, 
      panel.id=c('cell', 'year'))

etable(mod1)



#preparation of spatial weights matrix
# Create spatial neighbors based on distance
geoconflict_main_if <- geoconflict_main %>% 
  filter(!year == 1998)

coords <- cbind(geoconflict_main_if$lon, geoconflict_main_if$lat)
neighbors <- spdep::dnearneigh(coords, 0, 180, longlat = TRUE)

# Create the weights matrix using binary weights
weights <- spdep::nb2mat(neighbors, style = "B", zero.policy = TRUE)

## Model 3
contiguity_matrix <- mOrdNbr(raster, m =1)

## Setting a 180 km cutoff
dist_matrix <- st_distance(df$geometry)
dist_threshold <- set_units(180000, "m")
dist_threshold_matrix <- ifelse(dist_matrix > dist_threshold, 0, 1)

# Remove diagonal elements
diag(dist_threshold_matrix) <- 0

binary_matrix <- contiguity_matrix * dist_threshold_matrix


mod2 <- feols(ANY_EVENT_ACLED ~ SPEI4pg + L1_SPEI4pg + L2_SPEI4pg + GSmain_ext_SPEI4pg + L1_GSmain_ext_SPEI4pg + L2_GSmain_ext_SPEI4pg + W_L2_GSmain_ext_SPEI4pg + W_L1_GSmain_ext_SPEI4pg + W_GSmain_ext_SPEI4pg + W_SPEI4pg + W_L1_SPEI4pg + W_L2_SPEI4pg + elevation_cell +  rough_cell + area_cell + as.factor(use_primary) + 
                dis_river_cell + as.factor(shared) +  as.factor(border) + as.factor(any_mineral)+ ELF + i(country_largest_share, as.numeric(year)) + 
                W_elevation_cell + W_rough_cell + W_area_cell + W_ELF + as.character(W_any_mineral) + as.character(W_shared)  + W_dis_river_cell + as.character(W_use_primary) + i(country_largest_share, as.numeric(year)) | as.factor(year), 
              panel.id=c('cell', 'year'), 
              data = geoconflict_main)

etable(mod2)



### model 3

# spatial weights matrix 

geoconflict_main_weights <- geoconflict_main %>%
  filter(year == 2011) %>%
  st_as_sf(coords = c("lat", "lon"), crs = "WGS84")


neighboors <- dnearneigh(st_centroid(geoconflict_main_weights), d1 = 0, d2 = 180)
#neighboors <- poly2nb(geoconflict_main_weights, queen = TRUE)

listwmat <- nb2listw(neighboors, style ="W", zero.policy = TRUE)

spweightmat <- listw2mat(listwmat)

## Model 3
contiguity_matrix <- mOrdNbr(raster, m =1)

## Setting a 180 km cutoff
dist_matrix <- st_distance(geoconflict_main_weights$geometry)
dist_threshold <- set_units(180000, "m")
dist_threshold_matrix <- ifelse(dist_matrix > dist_threshold, 0, 1)

# Remove diagonal elements
diag(dist_threshold_matrix) <- 0

binary_matrix <- contiguity_matrix * dist_threshold_matrix

#country year fe
geoconflict_main_fe <- geoconflict_main_if %>% 
  mutate(country_year_fe = paste0(country_largest_share, year)) %>%
  arrange(cell)

form1 <- ANY_EVENT_ACLED ~ SPEI4pg + L1_SPEI4pg + GSmain_ext_SPEI4pg + L1_GSmain_ext_SPEI4pg + elevation_cell +  rough_cell + area_cell + as.factor(use_primary) + 
  dis_river_cell + as.factor(shared) +  as.factor(border) + as.factor(any_mineral)+ ELF 

#+ i(country_largest_share, as.numeric(year))
mod3  <- SDPDm(formula = ANY_EVENT_ACLED ~ SPEI4pg, data = as.data.frame(geoconflict_main_fe), 
               W = dist_threshold_matrix,
               index= c("cell","year"),
               effect = "time",
               ldet = "mc", 
               dynamic = TRUE)

mod5<-blmpSDPD(formula = ANY_EVENT_ACLED ~ SPEI4pg, data = as.data.frame(geoconflict_main_fe), 
            W = dist_threshold_matrix,
            index = c("cell","year"),
            model = "sdm", 
            effect = "twoways",
            LYtrans = T,
            dynamic = T,
            tlaginfo = list(ind = NULL, tl = T, stl = T))
summary(mod5)


# Convert data to panel_data format for ease of use
dat <- panel_data((geoconflict_main_fe), id = cell, wave = year)




