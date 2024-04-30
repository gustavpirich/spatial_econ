
# Exercise B - Replication of Table 2 -------------------------------------



# loading packages
pacman::p_load(spatialreg, fixest, splm, stringi, stringr, stringdist, haven, sf, dplyr, fuzzyjoin, 
               comparator, digest, zoomerjoin, ggplot2, tidyr, ggthemes, viridis, 
               fixest, conleyreg, plm, stargazer, magrittr, tidyverse, tmap, spdep, SDPDmod,
               igraph, generics, knitr, kableExtra, formatR,readxl, haven, units)


# change to local directory
raster_Africa <- read_sf("~/Desktop/GITHUB/spatial_econ/data/03_assignment/dataset/raster_Africa.shp")
geoconflict_main <- read_dta("~/Desktop/GITHUB/spatial_econ/data/03_assignment/dataset/geoconflict_main.dta")
intersect_coord <- read_dta("~/Desktop/GITHUB/spatial_econ/data/03_assignment/dataset/intersect_coord.dta")




# preparation of dataset --------------------------------------------------

# exclude year 1997 to obtain same number of observations
geoconflict_main_if <- geoconflict_main %>% 
  filter(!year == 1997)

# format as factor 
geoconflict_main_if$year <- as.factor(geoconflict_main_if$year)
geoconflict_main_if$cell <- as.factor(geoconflict_main_if$cell)


################
# 1st Column OLS

fml1 = ANY_EVENT_ACLED ~ SPEI4pg + L1_SPEI4pg + L2_SPEI4pg + GSmain_ext_SPEI4pg + L1_GSmain_ext_SPEI4pg + L2_GSmain_ext_SPEI4pg + elevation_cell + rough_cell + area_cell + as.factor(use_primary) + 
  dis_river_cell + as.factor(shared) +  as.factor(border) + as.factor(any_mineral) + ELF + i(country_largest_share, as.numeric(year), ref = "Zimbabwe") | as.factor(year)

mod1 <- feols(fml1, 
              data = geoconflict_main_if, 
              panel.id = c('cell', 'year'))

etable(mod1)
################


################
# 2nd Column OLS - missing spatially lagged linear time trends - also column 2 in the table looks like there are some coefficients/standard errors missing? 

fml2 = ANY_EVENT_ACLED ~ SPEI4pg + L1_SPEI4pg + L2_SPEI4pg + GSmain_ext_SPEI4pg + L1_GSmain_ext_SPEI4pg + L2_GSmain_ext_SPEI4pg + W_L2_GSmain_ext_SPEI4pg + W_L1_GSmain_ext_SPEI4pg + W_GSmain_ext_SPEI4pg + W_SPEI4pg + W_L1_SPEI4pg + W_L2_SPEI4pg + elevation_cell + 
  rough_cell + area_cell + as.factor(use_primary) + 
  dis_river_cell + as.factor(shared) +  as.factor(border) + as.factor(any_mineral)+ ELF + W_elevation_cell + W_rough_cell + W_area_cell + W_ELF + as.factor(W_any_mineral) + as.factor(W_shared)  + 
  W_dis_river_cell + as.factor(W_use_primary) + i(as.factor(country_largest_share), as.numeric(year), ref = "Zimbabwe") | as.factor(year)

mod2 <- feols(fml2, 
              data = geoconflict_main_if, 
              panel.id = c('cell', 'year'))

etable(mod2)
################


################
# 3rd Column  -----------------------------------------------------

# generate spatial weights matrix -----------------------------------------
geoconflict_main_weights <- geoconflict_main_if %>%
  st_as_sf(coords = c("lat", "lon"), crs = "WGS84") %>%
  filter(year == 2011) # could be any year

# setting the distance threshold to 180 kilometers
neighboors <- dnearneigh(st_centroid(geoconflict_main_weights), d1 = 0, d2 = 180, use_s2 = TRUE)

# create binary spatial weights matrix
listwmat <- nb2listw(neighboors, style ="B", zero.policy = TRUE)

# specification 3 - we might be missing the country specific linear time trends of the neighbors 
mod3 <- spml(ANY_EVENT_ACLED ~ lag(ANY_EVENT_ACLED) + SPEI4pg + L1_SPEI4pg + L2_SPEI4pg + GSmain_ext_SPEI4pg + L1_GSmain_ext_SPEI4pg + L2_GSmain_ext_SPEI4pg + W_GSmain_ext_SPEI4pg +
                               W_L1_GSmain_ext_SPEI4pg + W_L2_GSmain_ext_SPEI4pg + W_SPEI4pg + W_L1_SPEI4pg + W_L2_SPEI4pg + elevation_cell + 
                               rough_cell + area_cell + as.factor(use_primary) + dis_river_cell + as.factor(shared) +  as.factor(border) + as.factor(any_mineral) +
                               ELF + W_elevation_cell + W_rough_cell + W_area_cell + W_ELF + as.factor(W_any_mineral) + as.factor(W_shared)  + 
                               W_dis_river_cell + as.factor(W_use_primary) + as.factor(country_largest_share):as.numeric(year), # country specific linear time trends
                             data = as.data.frame(geoconflict_main_if), 
                             index=c("cell","year"),
                             listw = listwmat,
                             model="pooling",
                             effect = "time",  
                             dynamic = TRUE,
                             spatial.error="none",
                             zero.policy = TRUE, 
                             lag=TRUE)

summary(mod3)
################



################
# 4th column 

# Preparation adding column for country x year fixed effects (country_year_fe)
geoconflict_main_fe <- geoconflict_main_if %>% 
  mutate(country_year_fe = paste0(country_largest_share, year))

mod4 <- spml(ANY_EVENT_ACLED ~ lag(ANY_EVENT_ACLED) + SPEI4pg + L1_SPEI4pg + L2_SPEI4pg + GSmain_ext_SPEI4pg + L1_GSmain_ext_SPEI4pg + L2_GSmain_ext_SPEI4pg + W_GSmain_ext_SPEI4pg +
                              W_L1_GSmain_ext_SPEI4pg + W_L2_GSmain_ext_SPEI4pg + W_SPEI4pg + W_L1_SPEI4pg + W_L2_SPEI4pg + elevation_cell + 
                              rough_cell + area_cell + as.factor(use_primary) + dis_river_cell + as.factor(shared) +  as.factor(border) + as.factor(any_mineral) +
                              ELF + W_elevation_cell + W_rough_cell + W_area_cell + W_ELF + as.factor(W_any_mineral) + as.factor(W_shared)  + 
                              W_dis_river_cell + as.factor(W_use_primary) + as.factor(country_year_fe), # country x year FE
                             data = as.data.frame(geoconflict_main_fe), 
                             index= c("cell","year"),
                             listw = listwmat,
                             model = "pooling",  # no fixed effects
                             effect = "individual", 
                             spatial.error="none",
                             zero.policy = TRUE, 
                             dynamic = TRUE,
                             lag = TRUE, 
                             Hess = TRUE,
                             local=list( parallel = T)) # makes it faster

summary(mod4)


################ Does not work computationally
# 5th column 

mod5 <- spml(ANY_EVENT_ACLED ~ lag(ANY_EVENT_ACLED) + SPEI4pg + L1_SPEI4pg + L2_SPEI4pg + GSmain_ext_SPEI4pg + L1_GSmain_ext_SPEI4pg + L2_GSmain_ext_SPEI4pg + W_GSmain_ext_SPEI4pg +
               W_L1_GSmain_ext_SPEI4pg + W_L2_GSmain_ext_SPEI4pg + W_SPEI4pg + W_L1_SPEI4pg + W_L2_SPEI4pg + elevation_cell + 
               rough_cell + area_cell + as.factor(use_primary) + dis_river_cell + as.factor(shared) +  as.factor(border) + as.factor(any_mineral) +
               ELF + W_elevation_cell + W_rough_cell + W_area_cell + W_ELF + as.factor(W_any_mineral) + as.factor(W_shared)  + 
               W_dis_river_cell + as.factor(W_use_primary) + country_year_fe, 
             data = as.data.frame(geoconflict_main_fe), 
             index= c("cell","year"),
             listw = listwmat,
             model = "within",  # CELL fixed effects
             effect = "individual", 
             spatial.error="none",
             zero.policy = TRUE, 
             dynamic = TRUE,
             lag = TRUE, 
             Hess = TRUE,
             local=list( parallel = T)) # makes it faster

summary(mod5)







