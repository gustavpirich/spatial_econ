library(splm)
data("usaww")
usalw <- mat2listw(usaww, style="W")

# `splm` provides multiple useful functions, like:    

#  `bsjktest`        - Baltagi et al. test for spatial error correlation      
#  `bsktest`         - Baltagi et al. test for spat.err. or random effects  
#  `sphtest`         - Hausmann test for spatial models
#  `rwtest`          - Pesaran test for spatial dependence
#  `slmtest`         - Anselin test for nested spatial dependence (error in lag or vice versa)

# Also other useful little helpers are available:  

#  `slag`            - Calculates spatial lag      
#  `vcov.splm`       - Returns variance-covariance matrix
#  `effects.splm`    - Extracts fixed effects

# Functions for estimations would be:  

#  `spml`            - ML estimation of spatial panels
#  `spgm`            - GMM estimation of spatial panels
#  `spreml`          - ML estimation w/ serially correlated errors  


# Notation:
# λ : Spatial autoregressive parameter
# σ² ("phi") : Random effects
# ρ : Spatial error component
# ψ (psi) : Serial correlation


# As a basic introduction, let's set up a base model. We will estimate a simple 
# production function, linking gross social product (`gsp`) to  private capital 
# stock (`pcap`), public capital (`pc`), employment (`emp`), and unemployment 
# rate (`unemp`).
fm3 <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp

# The weights matrix is still `48x48`, while our data has 816 observations (`n x T`)! 
# We assume the neighborhood structure to be constant over time, thus `usalw` has 
# no notion of time, it defines the cross-sectional relation of observations only. 
# The `splm` package fits this one-time weights matrix to our multi-period data by 
# using Kronecker multiplication (see slides from lecture).
dim(usaww)

## Testing in spatial panels

# Similar to non-panel models, we can run through some testing steps to determine
# the source of the spatial dependence (if there is one). Note however, again, 
# that the selection of the model should rather be guided by theory (what effects
# are you interested in, where do they emanate, ...?)


# `slmtest()` contains Anselin's tests for spatial dependence, the panel equivalent 
# of the `LMtest()` from earlier on. Note that without any prior modification, 
# `slmtest()` assumes that the data can be pooled, i.e. it fixes only the dimension
# of the W model to fit to the panel data. If you however suspect that there are 
# significant individual effects in your panel then you need to get rid of them 
# before, i.e., using a manual `within()` transformation.

# Spatial lag structure?
slmlag <- slmtest(fm3, data=Produc, listw=usalw, test="lml")
slmlag
# Spatial error structure?
slmerr <- slmtest(fm3, data=Produc, listw=usalw, test="lme")
slmerr

# Robust versions:
# Spatial lag allowing for spatial error?
slmle <- slmtest(fm3, data=Produc, listw=usalw, test="rlml")
slmle
# Spatial error allowing for spatial lag?
slmel <- slmtest(fm3, data=Produc, listw=usalw, test="rlme")
slmel


# In panel models, we usually include some interactive effects, in the form of 
# fixed or random effects. Akin to the Hausmann test in standard panels, we can
# use the spatial version of Mutl/Paffermayr, provided by the `sphtest()` function.
# Note that the alternative hypothesis is that _one model is inconsistent_ (the 
# RE model), so if this is the case, we should use the FE estimator.

# For error specification
spherr <- sphtest(x=fm3, data=Produc, listw=usalw, 
                  spatial.model="error", method="ML")
spherr

# For lag specification
sphlag <- sphtest(x=fm3, data=Produc, listw=usalw, 
                  spatial.model="lag", method="ML")
sphlag

# For both lag and error specification
spherrlag <- sphtest(x=fm3, data=Produc, listw=usalw, 
                     spatial.model="sarar", method="ML")
spherrlag

# Digging into the error model, we can use some tests to determine, whether 
# there is spatial correlation or some form of error heterogeneity due to the 
# random region effects (or both)

# "LMH" H0:  ρ=σ²=0? (RE and spatial correlation non-zero?)
lmjoint <- bsktest(x=fm3, data=Produc, listw=usalw, test="LMH")
lmjoint

# "LM1" H0:  σ²=0? assuming ρ=0 (Are there ind effects, assuming no spatial corr.?)
lm1 <- bsktest(x=fm3, data=Produc, listw=usalw, test="LM1")
lm1

# "LM2" H0:  ρ=0? assuming σ²=0 (Is there spatial corr., assuming no ind effects?)
lm2 <- bsktest(x=fm3, data=Produc, listw=usalw, test="LM2")
lm2

# "CLMlambda" H0: ρ=0? (Is there spatial correlation given the possibility of RE?)
clml <- bsktest(x=fm3, data=Produc, listw=usalw, test="CLMlambda")
clml

# "CLMmu" H0:  σ²=0? (Are there RE, given the possibility of spatial correlation?)
clmm <- bsktest(x=fm3, data=Produc, listw=usalw, test="CLMmu")
clmm

# In addition, we can check for serial correlation using `bsjktest()`, a refined
# version o `bsktest()`. It tests in addition for correlation along the time 
# dimension, which is often disregarded in spatial panel setups.

# "J" H0: ρ=σ²=ψ=0? (RE, spatial or serial correlation non-zero?) 
bsjktest(x=fm3, data=Produc, listw=usalw, test="J")

# "C.1" H0: ρ=0? assuming ψ≠0 (serial correlation) and σ²>0 (RE) 
bsjktest(x=fm3, data=Produc, listw=usalw, test="C.1")

# "C.2" H0: ψ=0? assuming ρ≠0 (spatial correlation) and σ²>0
bsjktest(x=fm3, data=Produc, listw=usalw, test="C.2")

# "C.3" H0: σ²=0? assuming ρ≠0 and ψ≠0
bsjktest(x=fm3, data=Produc, listw=usalw, test="C.3")


## Estimating spatial panel models

# Full model including random effects (`model=`), lagged dependent (`lag=TRUE`) 
# and spatial error (`spatial.error=`):  
sararre <- spml(fm3, Produc, listw = usalw, model = "random", lag = TRUE, 
                spatial.error = "b") # Baltagi errors
summary(sararre)

# But, following the Hausmann test, we opt for a fixed effects specification
sararfemod <- spml(fm3, Produc, listw=usalw, lag=TRUE, spatial.error="b", 
                   model="within", effect="individual") # unit foxed effects
summary(sararfemod)

# Alternatively we can specify time FE
sarartfemod <- spml(fm3, Produc, listw=usalw, lag=TRUE, spatial.error = "b", 
                    model="within", effect="time") 
summary(sarartfemod)

# or both
sararbfemod <- spml(fm3, Produc, listw=usalw, lag=TRUE, spatial.error = "b", 
                    model="within", effect="twoways") 
summary(sararbfemod)


## Checking impacts 

# For newer versions of splm, there is `impacts()` for panel models. 
# Make sure that your weights (`usalw`) are row-standardized (`style="W"`). 
imp <- splm::impacts(sararbfemod, listw=usalw, time=1986)
summary(imp, zstats=TRUE, short=TRUE)

# `effects()` extracts the individual FE, which are not visible in the `summary()`  
eff <- effects(sararbfemod)
eff

## Resorting to GMM when ML struggles
GM_sararfemod <- spgm(fm3, Produc, listw=usaww, model="within", lag=TRUE, 
                      spatial.error=FALSE)
summary(GM_sararfemod)
imp.gm <- splm::impacts(GM_sararfemod, listw=usalw, time=1986)
imp.gm






# Promsising --------------------------------------------------------------


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


dist_threshold <- 180  # distance in meters
nb <- dnearneigh(st_coordinates(data_sf), 0, dist_threshold)
w_bin_180 <- nb2listw(nb, style = "B", zero.policy = TRUE)

# Save or use the spatial weights matrix
saveRDS(w_bin_180, "W_bin_180.rds")
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
               igraph, generics, knitr, kableExtra, formatR,readxl, haven)


setwd("~/Desktop/GITHUB/spatial_econ/data/03_assignment/dataset")

raster_Africa <- read_sf("~/Desktop/GITHUB/spatial_econ/data/03_assignment/dataset/raster_Africa.shp")

geoconflict_main <- read_dta("~/Desktop/GITHUB/spatial_econ/data/03_assignment/dataset/geoconflict_main.dta")

intersect_coord <- read_dta("~/Desktop/GITHUB/spatial_econ/data/03_assignment/dataset/intersect_coord.dta")

library(fixest)

geoconflict_main$year <- as.factor(geoconflict_main$year)
geoconflict_main$cell <- as.factor(geoconflict_main$cell)

mod1 <- feols(ANY_EVENT_ACLED ~ SPEI4pg + L1_SPEI4pg + L2_SPEI4pg + GSmain_ext_SPEI4pg + L1_GSmain_ext_SPEI4pg + L2_GSmain_ext_SPEI4pg + elevation_cell + rough_cell + area_cell + as.factor(use_primary) + 
        dis_river_cell + as.factor(shared) +  as.factor(border) + as.factor(any_mineral) + (ELF) + i((country_largest_share), as.numeric(year)) | as.factor(year), 
      data = geoconflict_main, 
      panel.id=c('cell', 'year'))

etable(mod1)


mod2 <- feols(ANY_EVENT_ACLED ~ 1 + SPEI4pg + L1_SPEI4pg + L2_SPEI4pg + GSmain_ext_SPEI4pg + L1_GSmain_ext_SPEI4pg + L2_GSmain_ext_SPEI4pg + elevation_cell + rough_cell + area_cell + use_primary + 
                dis_river_cell + as.factor(shared) +  as.factor(border) + as.factor(any_mineral)+ ELF + i(country_largest_share, as.numeric(year)) + 
                W_elevation_cell + W_rough_cell + W_area_cell + W_ELF + W_any_mineral + 
                W_L1_GSmain_ext_SPEI4pg + W_L2_GSmain_ext_SPEI4pg + W_SPEI4pg + W_L1_SPEI4pg + W_L2_SPEI4pg | as.factor(year), 
              data = geoconflict_main, 
              panel.id=c('cell', 'year'))

etable(mod2)
