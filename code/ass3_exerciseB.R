
# Exercise B - Replication of Table 2 -------------------------------------



# loading packages
pacman::p_load(spatialreg, fixest, splm, stringi, stringr, stringdist, haven, sf, dplyr, fuzzyjoin, 
               comparator, digest, zoomerjoin, ggplot2, tidyr, ggthemes, viridis, 
               fixest, conleyreg, plm, stargazer, magrittr, tidyverse, tmap, spdep, SDPDmod,
               igraph, generics, knitr, kableExtra, formatR,readxl, haven, units, data.table)


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

# visualizing the spillovers and impacts ----------------------------------

# For newer versions of splm, there is `impacts()` for panel models. 
# Make sure that your weights (`usalw`) are row-standardized (`style="W"`). 
geoconflict_main_weights <- geoconflict_main_if %>%
  st_as_sf(coords = c("lat", "lon"), crs = "WGS84")
neighboors <- dnearneigh(st_centroid(geoconflict_main_weights), d1 = 0, d2 = 180, use_s2 = TRUE)

listwmat <- nb2listw(neighboors, style ="W", zero.policy = TRUE)

imp <- impacts(mod3, listw=listwmat, time=1998)
summary(imp, zstats=TRUE, short=TRUE)



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
               W_L1_GSmain_ext_SPEI4pg + W_L2_GSmain_ext_SPEI4pg + W_SPEI4pg + W_L1_SPEI4pg + W_L2_SPEI4pg  + country_year_fe, 
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






# translation of R code ---------------------------------------------------

# Set up
rm(list = ls())  # Clear workspace
library(data.table)

name <- "figure2a"
cdate <- format(Sys.Date(), "%Y.%m.%d")

# Open log file
log_file <- file.path("${logs}", paste0(name, "_", cdate, ".log"))
sink(log_file)

# Coefficients
Y_lag <- 0.121
W_Y <- 0.0229
GSmain_ext_SPEI4pg <- mod3$coefficients[["GSmain_ext_SPEI4pg"]]
L1_GSmain_ext_SPEI4pg <- -0.0367
L2_GSmain_ext_SPEI4pg <- -0.00925
W_GSmain_ext_SPEI4pg <- -0.0042
W_L1_GSmain_ext_SPEI4pg <- 0.00648
W_L2_GSmain_ext_SPEI4pg <- 0.000522

# Load data
data <- geoconflict_main
data <- data %>%
  filter(lat >= -6.5 & lat < 2.5 & lon > 19.5 & lon < 29.5) %>%
  st_as_sf(coords = c("lat", "lon"))


# Weighting matrix
W_bin_180 <- dnearneigh(data, 0,180)

# Set initial values
data[, c(
  paste0(c("ANY_EVENT_cell", "SPEI4pg", "L1_SPEI4pg", "L2_SPEI4pg", "GSmain_ext_SPEI4pg", 
           "L1_GSmain_ext_SPEI4pg", "L2_GSmain_ext_SPEI4pg", "W_SPEI4pg", "W_L1_SPEI4pg",
           "W_L2_SPEI4pg", "W_GSmain_ext_SPEI4pg", "W_L1_GSmain_ext_SPEI4pg", 
           "W_L2_GSmain_ext_SPEI4pg"), "_0") := 0
)]

# Generate shock
data[, shock := ifelse(lon == 24.5 & lat == -2.5, 1, 0)]
data[, GSmain_ext_SPEI4pg_2 := ifelse(shock == 1, -0.365, 0)]

# Loop through time steps
for (j in 2:21) {
  j_lag1 <- j - 1
  j_lag2 <- j - 2
  j_lead1 <- j + 1
  
  # Update own cell
  data[, `:=`(
    GSmain_ext_SPEI4pg_lead1 := 0,
    Y_lag := shift(ANY_EVENT_cell, n = j_lag1),
    GSmain_ext_SPEI4pg := get(paste0("GSmain_ext_SPEI4pg_", j)),
    L1_GSmain_ext_SPEI4pg := get(paste0("GSmain_ext_SPEI4pg_", j_lag1)),
    L2_GSmain_ext_SPEI4pg := get(paste0("GSmain_ext_SPEI4pg_", j_lag2))
  )]
  
  # Generate spatial lags
  var_to_lag <- c("GSmain_ext_SPEI4pg", "L1_GSmain_ext_SPEI4pg", "L2_GSmain_ext_SPEI4pg")
  lag_data <- data[, ..var_to_lag, keyby = .(cell)][, lapply(.SD, shift, n = 1), .SDcols = var_to_lag]
  
  # Merge spatial lags to main dataset
  data <- merge(data, lag_data, by = "cell", all.x = TRUE)
  
  # Autoregressive term
  data[, `:=`(
    term1 := Y_lag^2,
    term5 := GSmain_ext_SPEI4pg^2,
    term6 := L1_GSmain_ext_SPEI4pg^2,
    term7 := L2_GSmain_ext_SPEI4pg^2,
    term11 := W_GSmain_ext_SPEI4pg^2,
    term12 := W_L1_GSmain_ext_SPEI4pg^2,
    term13 := W_L2_GSmain_ext_SPEI4pg^2
  )]
  
  data[, Ytilde := rowSums(.SD), .SDcols = paste0("term", 1:13)]
  
  # Calculate Y2
  I <- diag(nrow(data))
  IminusW <- I - W_Y * as.matrix(W_bin_180)
  Y2 <- solve(IminusW) %*% as.matrix(data$Ytilde)
  data[, paste0("Y_pred_", j) := Y2]
}

# Keep necessary columns
data <- data[, c("cell", "lat", "lon", grep("Y_pred_", names(data), value = TRUE), "shock", "_ID")]

# Plot Figure 2a
shock_data <- data[shock == 1]
shock_data[, t := 1:.N]
shock_data[, W_Y_pred_ := W_Y_pred_ / 7.4]
shock_data <- shock_data[t <= 7]
shock_data[, year := t - 2]

plot_data <- shock_data[, .(mean_Y_pred_ = mean(Y_pred_), mean_W_Y_pred_ = mean(W_Y_pred_)), by = "year"]
plot_data <- melt(plot_data, id.vars = "year")
plot_data$year <- as.factor(plot_data$year)

library(ggplot2)
ggplot(plot_data, aes(x = year, y = value, color = variable)) +
  geom_line(size = 1) +
  scale_x_discrete(name = "Years since shock") +
  scale_y_continuous(labels = scales::comma_format()) +
  scale_color_manual(values = c("black", "black"), labels = c("cell with shock", "neighboring cell")) +
  theme_minimal() +
  theme(legend.position = "top") +
  labs(title = "Figure 2a")

# Plot Figure 2b
years <- c(2, 3, 5, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14)
pdf("${output}/Figure2b.pdf")
par(mfrow = c(4, 4))
for (i in years) {
  spmap(Y_pred_[[i]], ...)
}
dev.off()

# Close log
sink()
library(data.table)
library(spdep)
library(ggplot2)
library(dplyr)

# Coefficients
Y_lag <- 0.121
W_Y <- 0.0229
GSmain_ext_SPEI4pg <- 0.0207
L1_GSmain_ext_SPEI4pg <- -0.0367
L2_GSmain_ext_SPEI4pg <- -0.00925
W_GSmain_ext_SPEI4pg <- -0.0042
W_L1_GSmain_ext_SPEI4pg <- 0.00648
W_L2_GSmain_ext_SPEI4pg <- 0.000522

# Load data
data <- fread("geoconflict_main.dta")

# Filter cells
data <- data[lat >= -6.5 & lat < 2.5 & lon > 19.5 & lon < 29.5]

# Weighting matrix
coords <- data[, c("lon", "lat")]
W_bin_180 <- nb2listw(spdep::dnearneigh(coords, 0, 180), style = "B")

# Set up variables
data[, c("Y_lag", "GSmain_ext_SPEI4pg", "L1_GSmain_ext_SPEI4pg", "L2_GSmain_ext_SPEI4pg") := .(NA, NA, NA, NA)]

# Generate shock
data[, shock := ifelse(lon == 24.5 & lat == -2.5, 1, 0)]
data[shock == 1, GSmain_ext_SPEI4pg_2 := -0.365]

# Generate lagged variables
for (j in 2:21) {
  j_lag1 <- j - 1
  j_lag2 <- j - 2
  j_lead1 <- j + 1
  
  data[Y_lag := shift(ANY_EVENT_cell, n = j_lag1, fill = NA)]
  data[GSmain_ext_SPEI4pg := shift(GSmain_ext_SPEI4pg, n = j, fill = NA)]
  data[L1_GSmain_ext_SPEI4pg := shift(GSmain_ext_SPEI4pg, n = j_lag1, fill = NA)]
  data[L2_GSmain_ext_SPEI4pg := shift(GSmain_ext_SPEI4pg, n = j_lag2, fill = NA)]
  
  # Spatial lags
  var_to_lag <- c("GSmain_ext_SPEI4pg", "L1_GSmain_ext_SPEI4pg", "L2_GSmain_ext_SPEI4pg")
  for (var in var_to_lag) {
    lagged <- data[, get(var)]
    W_lagged <- lagged %*% W_bin_180
    names(W_lagged) <- paste0("W_", var)
    data <- cbind(data, W_lagged)
  }
  
  # Autoregressive term
  terms <- c("Y_lag^2", "GSmain_ext_SPEI4pg^2", "L1_GSmain_ext_SPEI4pg^2", "L2_GSmain_ext_SPEI4pg^2",
             "W_GSmain_ext_SPEI4pg^2", "W_L1_GSmain_ext_SPEI4pg^2", "W_L2_GSmain_ext_SPEI4pg^2")
  data[, Ytilde := rowSums(.SD, na.rm = TRUE), .SDcols = terms]
  IminusW_inv <- solve(diag(nrow(data)) - W_Y * W_bin_180)
  Y2 <- as.matrix(data$Ytilde) %*% IminusW_inv
  data[[paste0("Y_pred_", j)]] <- Y2
}

# Figure 2a
shock_data <- data[shock == 1]
shock_data <- melt(shock_data, id.vars = c("lat", "lon"))
shock_data[, W_Y_pred_ := W_Y_pred_ / 7.4]
shock_data <- shock_data[1:7, ]  # Keep data for first 7 years
shock_data[, year := 1:7]
plot_2a <- ggplot(shock_data, aes(x = year)) +
  geom_line(aes(y = Y_pred_, linetype = "cell with shock"), color = "black") +
  geom_line(aes(y = W_Y_pred_, linetype = "neighboring cell"), color = "black", linetype = "dashed") +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.005), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(0, 5, by = 1), expand = c(0, 0), name = "Years since shock") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(title = "Figure 2a") +
  scale_linetype_manual(name = "", values = c("solid", "dashed"), labels = c("cell with shock", "neighboring cell"))
ggsave("Figure2a.pdf", plot_2a)

# Figure 2b
year_periods <- c(2, 3, 5, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14)
plots <- lapply(year_periods, function(i) {
  period <- i - 2
  sp_map <- ggplot(data, aes(x = lon, y = lat)) +
    geom_polygon(aes(group = group), fill = "white", color = "black") +
    geom_point(data = data, aes(color = get(paste0("Y_pred_", i))), size = 1) +
    scale_color_gradientn(colors = colorRampPalette(c("white", "blue"))(9),
                          breaks = c(-0.008, 0.00005, 0.0001, 0.00015, 0.0003, 0.0005, 0.001, 0.003, 0.0075, 0.10),
                          labels = format(c(-0.008, 0.00005, 0.0001, 0.00015, 0.0003, 0.0005, 0.001, 0.003, 0.0075, 0.10), nsmall = 5),
                          name = paste("t =", period)) +
    theme_void() +
    theme(legend.position = "bottom") +
    labs(title = paste("Figure 2b - t =", period))
  return(sp_map)
})
plot_2b <- do.call(gridExtra::grid.arrange, plots)
ggsave("Figure2b.pdf", plot_2b)


mod3 <- readRDS("./models/mod3.rds")

mod4 <- readRDS("./models/mod4.rds")
# how can we visualize the spatial spillover effects
listwmat <- nb2listw(neighboors, style ="B", zero.policy = TRUE)

imp <- splm::impacts(mod4, listw=listwmat, time = "year")
summary(imp, zstats=TRUE, short=TRUE)
mod4

w_sparse <- as(as_dgRMatrix_listw(lw), "CsparseMatrix")

# Let's do the calculation of impacts by hand to see where they come from
lambda <- as.numeric(mod4$arcoef)
beta.hat <- coefficients(mod4)[-1]
N <- nrow(geoconflict_main)
I <- diag(N)
ones <- as.matrix(rep(1, N))
t_ones <- t(ones)
w_sparse <- as(as_dgRMatrix_listw(listwmat), "CsparseMatrix")

S_inv <- solve(I - lambda * w_sparse) # inverse of spatial filter


# Using the inversion, we can build the matrix of partial derivatives dy/dx
# Using the example of the estimated parameter for income estimate, we scale its 
# marginal effect beta by the spatial structure given by (I-lambda W)⁻¹ (S_inv)
S_W <-  S_inv %*% (I * beta.hat["INC"])

# The direct effect is the sum of the diagonal scaled by `N`
S_W_income_direct <- sum(diag(S_W)) / N
S_W_income_direct

# Total effects are the sum of all derivatives (`sum(S_W)`) or in matrix notation:  
S_W_income_tot <- (t_ones %*% (S_W) %*% ones) / N
S_W_income_tot

# The indirect effect is their difference (sum of all off-diagonal elements).  
S_W_income_indirect <- S_W_income_tot - S_W_income_direct
S_W_income_indirect



# Alternative with traces (better for large W)
W <- as(listwmat, "CsparseMatrix")
trMatc <- trW(W, type = "mult",
              m = 30) # number of powers
mod_1.sar.imp2 <- impacts(mod4, 
                          tr = trMatc, # trace instead of listw
                          R = 300, 
                          Q = 30) # number of power series used for approximation
summary(mod_1.sar.imp2, zstats = TRUE, short = TRUE)

# Number of years
T <- length(unique(geoconflict_main$year))

# impacts
sarre.mod.imp <- impacts(mod4,
                         listw = listwmat,
                         time = T)
summary(sarre.mod.imp)                         






# DERIVATION OF SHOCKS

#first number denotes period second number denotes neighboor
# shock in period 0 of standard deviation in growing season SEPI (without contemporanously affecting SEPI over the year)
risk_cell_0_0 <- -1*mod4$coefficients["GSmain_ext_SPEI4pg"]*sd(geoconflict_main_weights$GSmain_ext_SPEI4pg)

risk_cell_0_1 <- -1*mod4$coefficients["W_GSmain_ext_SPEI4pg"]*sd(geoconflict_main_weights$GSmain_ext_SPEI4pg) + risk_cell_0_0*mod4$arcoef 

risk_cell_0_2 <- -1*mod4$coefficients["W_GSmain_ext_SPEI4pg"]*mod4$coefficients["W_GSmain_ext_SPEI4pg"]*sd(geoconflict_main_weights$GSmain_ext_SPEI4pg) + mod4$arcoef*risk_cell_0_1 + mod4$arcoef*mod4$arcoef*risk_cell_0_0


risk_cell_1_0 <- risk_cell_0_0*mod4$coefficients["lag(ANY_EVENT_ACLED)"] + -1*mod4$coefficients["L1_GSmain_ext_SPEI4pg"]*sd(geoconflict_main_weights$GSmain_ext_SPEI4pg) 

risk_cell_1_1 <- risk_cell_0_1*mod4$coefficients["lag(ANY_EVENT_ACLED)"] + -1*mod4$coefficients["W_L1_GSmain_ext_SPEI4pg"]*sd(geoconflict_main_weights$GSmain_ext_SPEI4pg) + risk_cell_1_0*mod4$arcoef 

risk_cell_1_2 <- risk_cell_0_2*mod4$coefficients["lag(ANY_EVENT_ACLED)"] + -1*mod4$coefficients["W_L1_GSmain_ext_SPEI4pg"]*mod4$coefficients["W_L1_GSmain_ext_SPEI4pg"]*sd(geoconflict_main_weights$GSmain_ext_SPEI4pg) + mod4$arcoef*risk_cell_1_1 + mod4$arcoef*mod4$arcoef*risk_cell_1_0


risk_cell_2_0 <- risk_cell_1_0*mod4$coefficients["lag(ANY_EVENT_ACLED)"] + risk_cell_0_0 * mod4$coefficients["lag(ANY_EVENT_ACLED)"] *mod4$coefficients["lag(ANY_EVENT_ACLED)"]  + -1*mod4$coefficients["L2_GSmain_ext_SPEI4pg"]*sd(geoconflict_main_weights$GSmain_ext_SPEI4pg)

risk_cell_2_1 <- risk_cell_1_1*mod4$coefficients["lag(ANY_EVENT_ACLED)"] + risk_cell_0_1 * mod4$coefficients["lag(ANY_EVENT_ACLED)"] *mod4$coefficients["lag(ANY_EVENT_ACLED)"]  + -1*mod4$coefficients["W_L2_GSmain_ext_SPEI4pg"]*sd(geoconflict_main_weights$GSmain_ext_SPEI4pg) + risk_cell_2_0*mod4$arcoef 

risk_cell_2_2 <- -1*mod4$coefficients["W_L2_GSmain_ext_SPEI4pg"]*mod4$coefficients["W_L2_GSmain_ext_SPEI4pg"]*sd(geoconflict_main_weights$GSmain_ext_SPEI4pg) + mod4$arcoef*risk_cell_2_1 +  mod4$arcoef*mod4$arcoef*risk_cell_2_0 + risk_cell_2_1*mod4$arcoef +risk_cell_2_0*mod4$arcoef*mod4$arcoef


grid1 <- matrix(NA, nrow = 5, ncol = 5)
grid2 <- matrix(NA, nrow = 5, ncol = 5)
grid3 <- matrix(NA, nrow = 5, ncol = 5)

# Define the values for each cell
risk_cell_0_0 
risk_cell_0_1 
risk_cell_0_2 

# Fill the outer cells with risk_cell_0_2
grid1[, c(1, 5)] <- risk_cell_0_2
grid1[c(1, 5), ] <- risk_cell_0_2

# Fill the inner cells with risk_cell_0_1
grid1[-c(1, 5), -c(1, 5)] <- risk_cell_0_1

# Fill the central cell with risk_cell_0_0
grid1[3, 3] <- risk_cell_0_0

# Convert the data to a data frame
df <- expand.grid(x = 1:5, y = 1:5)
df$z <-  as.vector(grid1)   # Replace with your own values

# Plot the grid with ggplot2
s1 <- ggplot(df, aes(x = x, y = y, fill = z)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +  # Adjust colors as needed
  theme_void() +
  labs(title = "Shock Growing Season SEPI for central cell at time 0")



# Define the values for each cell
risk_cell_1_0 
risk_cell_1_1 
risk_cell_1_2 

# Fill the outer cells with risk_cell_0_2
grid2[, c(1, 5)] <- risk_cell_1_2
grid2[c(1, 5), ] <- risk_cell_1_2

# Fill the inner cells with risk_cell_0_1
grid2[-c(1, 5), -c(1, 5)] <- risk_cell_1_1

# Fill the central cell with risk_cell_0_0
grid2[3, 3] <- risk_cell_1_0

# Convert the data to a data frame
df <- expand.grid(x = 1:5, y = 1:5)
df$z <- as.vector(grid2)   # Replace with your own values

# Plot the grid with ggplot2
s2 <- ggplot(df, aes(x = x, y = y, fill = z)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +  # Adjust colors as needed
  theme_void() +
  labs(title = "Standard Deviation Growing Season SEPI Shock after 1 year")



# Define the values for each cell
risk_cell_2_0 
risk_cell_2_1 
risk_cell_2_2 

# Fill the outer cells with risk_cell_0_2
grid3[, c(1, 5)] <- risk_cell_2_2
grid3[c(1, 5), ] <- risk_cell_2_2

# Fill the inner cells with risk_cell_0_1
grid3[-c(1, 5), -c(1, 5)] <- risk_cell_2_1

# Fill the central cell with risk_cell_0_0
grid3[3, 3] <- risk_cell_2_0

# Convert the data to a data frame
df <- expand.grid(x = 1:5, y = 1:5)
df$z <- as.vector(grid3)   # Replace with your own values

# Plot the grid with ggplot2
s3 <- ggplot(df, aes(x = x, y = y, fill = z)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +  # Adjust colors as needed
  theme_void() +
  labs(title = "Shock Growing Season SEPI after two years")






library(ggplot2)
library(gridExtra)


# Combine all three grids into one list
all_grids <- list(grid1, grid2, grid3)

# Calculate the overall range of values across all grids
overall_range <- range(unlist(all_grids))

# Convert data for each grid to data frames
df1 <- expand.grid(x = 1:5, y = 1:5)
df1$Prob <- as.vector(grid1)

df2 <- expand.grid(x = 1:5, y = 1:5)
df2$Prob <- as.vector(grid2)

df3 <- expand.grid(x = 1:5, y = 1:5)
df3$Prob <- as.vector(grid3)

# Plot grids next to each other with a consistent color scale using viridis palette
grid.arrange(
  ggplot(df1, aes(x = x, y = y, fill = Prob)) +
    geom_tile() +
    scale_fill_viridis_c(option = "magma", limits = overall_range) +
    theme_void() +
    theme(legend.position = "none") +
    labs(title = "SEPI Growing Season Shock at t = 0"),
  
  ggplot(df2, aes(x = x, y = y, fill = Prob)) +
    geom_tile() +
    scale_fill_viridis_c(option = "magma", limits = overall_range) +
    theme_void() +
    theme(legend.position = "none") +
    labs(title = "t = 1"),
  
  ggplot(df3, aes(x = x, y = y, fill = Prob)) +
    geom_tile() +
    scale_fill_viridis_c(option = "magma", limits = overall_range) +
    theme_void() +
    labs(title = "t = 2", fill = "Pr(Conflict)"),
  nrow = 1
)






