#------------------------------------------>
#- Spatial Economics
#-      Tutorial 2: Spatial stuff
#-
#- Original Authors (i.e. many thanks to):
#-  Mathias Moser (matmoser@wu.ac.at)
#-  Franziska Disslbacher
#-    (Franziska.disslbacher@wu.ac.at)
#- 
#-  Adapted by:
#-  Lukas Vashold
#-    (lukas.vashold@wu.ac.at)
#-  Nikolas Kuschnig
#-    (nikolas.kuschnig@wu.ac.at)
#-  Martin Prinz
#------------------------------------------>

library(dplyr)

# create folder for the current unit/project
dir.create("spatial_unit2")
setwd("spatial_unit2")

# create subfolders for code and data respectively
dir.create("code")
dir.create("data")


# Useful packages when working with spatial data: -------------------------


# sf (simple features):

  #simple features to a formal standard that describes how objects in the real 
  #world can be represented in computers, with emphasis on the spatial geometry of these objects
  
  # Simple features is a widely supported data model 
  #that underlies data structures in many GIS applications 
  #including QGIS and PostGIS. A major advantage of this is 
  #that using the data model ensures your work is cross-transferable 
  #to other set-ups, for example importing from and exporting to 
  #spatial databases based on the idea that all geometry are based on points: 
  #x and y coordinates + z  coordinate (for altitude)  + m coordinate
  ##=vecotr data principle: The geographic vector data model is 
  #based on points located within a coordinate reference system (CRS).
  
  #Advantages: 
  #Fast reading and writing of data.
  #Enhanced plotting performance
  #sf objects can be treated as data frames in most operations.
  #sf functions can be combined using %>% operator and works well with the tidyverse collection of R packages.
  #sf function names are relatively consistent and intuitive (all begin with st_).
  
  #++commands: st_......() (st = spatio-temporal)
  #%>% sf functions can be combined with pipes
  
# maptools:

  #Set of tools for manipulating and reading geographic data.

# rgdal:

  #Provides bindings to the 'Geospatial' Data Abstraction Library ('GDAL'). This builds the basis
  # for many geographic operations (like st_transform from last unit).
  # NOTE: This is outdated and will be phased out and integrated directly into the sf package

# spdep:

  #A collection of functions to create spatial weights matrix objects from polygon 
  #'contiguities', from point patterns by distance and tessellations, for summarizing these 
  #objects, and for permitting their use in spatial data analysis; a collection of tests 
  #for spatial 'autocorrelation'; exact tests for global and local 'Moran's I'; and functions 
  #for estimating spatial simultaneous autoregressive' ('SAR') lag and error models, impact 
  #measures for lag models, weighted and 'unweighted' 'SAR' and 'CAR' spatial regression models.

# maptools:

  #Alternative set of tools for manipulating and reading geographic data, in particular 'ESRI Shapefiles'; 
  #C code used from 'shapelib'. We will focus on "sf" as our primary package for these tasks.

# maps:

  #Display of maps. Projection code and larger maps are in separate packages.

# rgeos:

  #Interface to Geometry Engine

# GISTools:

  # Some mapping and spatial data manipulation tools - in particular drawing choropleth maps 
  # with nice looking legends, and aggregation of point data to polygons.

# tmap + tmaptools:

  # Packages for creating thematic maps, alternative to ggplot2 as visualization tool



# Shapefiles --------------------------------------------------------------


# The shapefile contains data for most of the spatial applications, as points, lines, and polygons.

# create a new empty object called 'temp' in which to store a zip file
# containing boundary data
temp <- tempfile(fileext = ".zip")

download.file("http://ec.europa.eu/eurostat/cache/GISCO/distribution/v2/nuts/download/ref-nuts-2021-03m.shp.zip", temp)

outDir <- "./data"
unzip(temp, exdir=outDir)

library(sf)
list.files(path="./data", pattern = '*.shp')

# let us choose projection WGS 84 (EPSG 4326) which is visible the file name 
# between the year (2021) and the level of the data (NUTS 2):
unzip("./data/NUTS_RG_03M_2021_4326_LEVL_2.shp.zip", exdir=outDir)
shp <- st_read(dsn = "./data", layer ="NUTS_RG_03M_2021_4326_LEVL_2") #reads in the shapefile

# The shapefile looks just like a normal dataframe with an additional "geometry" column:
head(shp)

# Save crs for later
old_crs <- st_crs(shp)

# There are also other functions to read shapefiles, but st_read() covers most of the use cases
# You can check the difference by looking at class()
class(shp)

# If it is "sf" then you can use the functions from the sf package, which start with st_*
# If it is "sp" then you are using the legacy sp format, which we do not cover here

# Using the plot function we can create a plot of the shapefile
# st_geometry() is used to just plot the borders (not variables or such)
plot(st_geometry(shp))

# the map seems kind of "squeezed"

shp <- st_transform(shp, 
                   st_crs("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +no_defs"))

# we can also exclude all oversea territories
overseas <- c("FRY1", "FRY2", "FRY3", "FRY4", "FRY5", "FRZZ", 
              "PT20", "PT30", "PTZZ", 
              "ES70", "ESZZ", 
              "NO0B", "NOZZ")
shp <- shp[! shp$NUTS_ID %in% overseas, ]

plot(st_geometry(shp))

# different projections are in use: e.g.
  # WGS84 (the projection we used before)
      # is usually used by organizations which 
      # provide GIS data for the entire globe or many countries, it is also used by Google Earth
  # MGI Lambert (EPSG:31287)
      # shapefiles from Statistik Austria use this projection
  # UTM, Zone 10 (EPSG: 32610)
      # Zone 10 is used in the Pacific Northwest
# for more details, see https://www.nceas.ucsb.edu/~frazier/RSpatialGuides/OverviewCoordinateReferenceSystems.pdf

plot(st_geometry(st_transform(shp, 
                 st_crs("epsg:31287"))))

# The CRS can also mirror/rotate the map!
plot(st_geometry(st_transform(shp, 
                 st_crs("epsg:32610"))))



# Types of spatial data ---------------------------------------------------

class(shp)

# The sf class can represent various types of spatial data:
#   POINT: A point in space
#   LINESTRING: Points connected by straight lines
#   POLYGON: Sequence of points that create closed area
#   MULTI-/POINT/LINESTRING/...: Collections of multiple POINTS, LINESTRINGS, ...

# For example our shp is a MULTIPOLYGON (many POLYGONs, one for each region)
class(shp$geometry)

# POINTS Example
# Center-point of each polygon...
centroids <- st_centroid(shp)
class(centroids$geometry)     # The centroids are POINTs
plot(st_geometry(centroids), pch=".", cex = 2, col = "red")
plot(st_geometry(shp), add=TRUE, lwd=0.1)


# Let's select just Vienna (1 observation)
poly <- filter(shp, NUTS_NAME=="Wien")
class(poly$geometry)

# Extract the POINTs of the polygon between which the lines are drawn:
poly_points <- st_cast(poly, "POINT")
class(poly_points$geometry)

nrow(poly_points)   # To draw Vienna, this shapefile needs 30 points
plot(st_geometry(poly))
plot(st_geometry(poly_points), col="red", add=TRUE)



# Adding data to the shapefile --------------------------------------------

# install.packages("eurostat")
library(eurostat)

# Data table of contents
EurostatTOC <- get_eurostat_toc()
head(EurostatTOC)

# tgs00026 (Disposable income of private households)
df <- get_eurostat("tgs00026", time_format = "raw")
df$TIME_PERIOD <- eurotime2num(df$TIME_PERIOD)

df <- df[df$TIME_PERIOD == 2018, ]     # Difference to which()?

df$value_cat <- cut_to_classes(df$values, n = 5)

df <- df[!df$geo %in% overseas,]

# Combine regions with data
# For sf objects we use *_join, like in dplyr:
shp <- left_join(shp, df, by=c("NUTS_ID"="geo"))

# Calculate polygon centroids, and extract the coordinates as data.frame
coords <- st_coordinates(st_centroid(shp))

# Eurostat also provides shapefiles directly
shp_eurostat <- get_eurostat_geospatial(output_class="sf", resolution="20", nuts_level=2, year=2021)



# Creating spatial weights matrices ---------------------------------------

# One can create spatial weights matrices or with own functions
# however, the package spdep also provides some routines for this task
library(spdep)
# 3 steps: poly ---poly2nb()--> neighbours list ----nb2listw--> weights list


# Contiguity

## Queen
queen_nb <- poly2nb(shp, row.names=shp$NUTS_ID, queen=TRUE)
W.list.queen <- nb2listw(queen_nb, style = "W")
# there are regions with no neighbors (e.g. Iceland), can be problematic
W.list.queen <- nb2listw(queen_nb, style = "W", zero.policy=TRUE)
W.queen <- listw2mat(W.list.queen)

# we can also plot weights matrices
plot(st_geometry(shp))
plot(queen_nb, coords, add=TRUE, col="green", cex=0.5)

## Rook
rook_nb <- poly2nb(shp, row.names=shp$NUTS_ID, queen=FALSE)

## Bishop (difference between the other twos)
bishop_nb <- diffnb(rook_nb, queen_nb)


# KNN (k-nearest-neighbours)
coords <- st_coordinates(st_centroid(shp))
k.near <- knearneigh(coords, k=5)
k5 <- knn2nb(k.near)
W.list.k <- nb2listw(k5, style = "W")

plot(st_geometry(shp))
plot(k5, coords, add=TRUE, col="green", cex=0.5)


# Distance band (in meters - depends on CRS)
distw <- dnearneigh(coords, 0, 1000000, row.names=shp$NUTS_ID)
summary(distw)

# We see that Iceland does not get a neighbour with 1000m distance
# What is the minimum distance we need so everybody gets at least one
# neighbor? We search for the nearest neighbor of each observation
# and then select the distance of those that are farthest apart
# (which is Norway-Iceland)
k1 <- knearneigh(coords, k=1)
k1 <- knn2nb(k1)
link.max <- max(unlist(nbdists(k1, coords=coords)))
link.max

distw <- dnearneigh(coords, 0, link.max, row.names=shp$NUTS_ID)
W.list <- nb2listw(distw, style="W", zero.policy=FALSE)

plot(st_geometry(shp))
plot(distw, coords, add=TRUE, col="green", cex=0.5)

# (Raw) Distances
dnbdists <- nbdists(distw, coords)
## Inverse distance
gl <- lapply(dnbdists, function(x) 1/x)
W.list.gl <- nb2listw(distw, glist=gl, zero.policy=FALSE)



# Handling spatial data ---------------------------------------------------
 
# install.packages("rgeos")
library(rgeos)

# the rgeos (stands for “R interface to the Geometry Engine - Open Source library”) package
# is a set of tools for geometric operations the you can either call directly or more
# comforatbly through the sf package.

# NOTE: rgeos will be retired in October 2023, with sf taking over. Left here as a reference

# Examples:
  # gArea: calculate area of a shape (st_area)
  # gDistance: distance between items (st_distance)
  # gCentroid: Collapse polygons to their centroids (st_centroid)
  
# units:
# rgeos isn’t a spatial library – it just does geometry. 
# It takes x-y coordinates and applies geometric formula. 
# Thus the units will always be the units of the x-y coordinates, which come from your projection. So here, we 
# can check the projection to find our units:

st_crs(shp)
# here, in our case, LENGTHUNIT["metre"] says our units are in meters

# area for each region:
st_area(shp) # computes area in each polygon

# for all areas at once - first "union" them together into one region
st_union(shp) %>% st_area()

# Dissolving is the process of dissolving sub-boundaries
# For example aggregating regions to a country.
# For this we use group_by and summary:
shp_country <- shp %>% 
  group_by(CNTR_CODE) %>% 
  summarise(CNTR_CODE = first(CNTR_CODE))
plot(st_geometry(shp_country))

# union of two shapefiles
# there are other ways to create subsets of a shapefile, e.g.
# let's create two shapefiles of Austria and Switzerland
shp_AT <- filter(shp, CNTR_CODE=="AT")
plot(st_geometry(shp_AT), main="Austria")

shp_CH <- filter(shp, CNTR_CODE=="CH")
plot(st_geometry(shp_CH), main="Switzerland")

# now, we can unify those two shapefiles using bind_rows() as in dplyr
shp_AT_CH <- bind_rows(shp_AT, shp_CH)
plot(st_geometry(shp_AT_CH))



# Export shapefiles -------------------------------------------------------


# if you want to save your shapefile which is now amended with income data:
dir.create("./shapefile EU")

# Save a shapefile (.shp) - beware this is an old format and has some limitations
# The error you get is because it does not support special characters in the column names
st_write(shp, dsn = "./shapefile EU/shapefile.shp")

# Better: Use GeoJSON (also faster and just one file)
st_write(shp, dsn = "./shapefile EU/shapefile.geojson")



# Issues in space ---------------------------------------------------------

# MAUP
## Modifiable Areal Unit Problem


# use now only those units where we have observations of the variable of interest 
# excluding e.g. Albania, Turkey, Switzerland
shp <- shp[which(!is.na(shp$values)), ]

# have to create neighbor lists anew given that we excluded regions
queen_nb <- poly2nb(shp, row.names=shp$NUTS_ID, queen=TRUE)
W.list.queen <- nb2listw(queen_nb, style = "W", zero.policy=TRUE)

coords <- st_coordinates(st_centroid(shp))
k1 <- knearneigh(coords, k=1)
k1 <- knn2nb(k1)
link.max <- max(unlist(nbdists(k1, coords=coords)))
distw <- dnearneigh(coords, 0, link.max, row.names=shp$NUTS_ID)
W.list <- nb2listw(distw, style="W", zero.policy=FALSE)

k.near <- knearneigh(coords, k=5)
k5 <- knn2nb(k.near)
W.list.k <- nb2listw(k5, style = "W")

## NUTS2:
cor(shp$values, lag.listw(W.list.queen, shp$values, NAOK=TRUE), use="complete.obs")

## NUTS1:
# Create the new, shorter ID for NUTS1 regions
shp$id_nuts1 <- substr(shp$NUTS_ID, 1, 3)

# Before proceeding we "correct" the shapefile in case some boundaries overlap (often not needed)
shp <- st_buffer(shp, 0)
# "Union" the shapefile by this ID and aggregate the data to NUTS1 (mean)
shp_nuts1 <- shp %>% group_by(id_nuts1) %>% summarise(values=mean(values))

W.list.queen.n3 <- nb2listw(poly2nb(shp_nuts1, 
                                    row.names = shp_nuts1$id_nuts1, queen=TRUE), 
                            style = "W", zero.policy=TRUE)

# lower correlation when we aggregate from NUTS2 to NUTS1
cor(shp_nuts1$values, lag.listw(W.list.queen.n3, shp_nuts1$values, NAOK=TRUE), use="complete.obs")


# Spatial Autocorrelation

### GLOBAL Moran's I ###

# Moran's test for spatial autocorrelation using a spatial weights matrix in 
# weights list form. The assumptions underlying the test are sensitive to the 
# form of the graph of neighbour relationships and other factors. 
# Results may be checked against those of moran.mc permutations.

# I = (n sum_i sum_j w_ij (x_i - xbar) (x_j - xbar)) / (S0 sum_i (x_i - xbar)^2)

# Computing the  Moran's I manually
WS <- listw2mat(W.list) # transform our weights list into an W matrix
u <- shp$values # choose the variable of interest
u.mean <- mean(u) # calculate the mean
MI.u.num <- (u - u.mean)%*%WS%*%(u - u.mean) # calculate the numerator of Moran's I
MI.u.den <- (u - u.mean)%*%(u - u.mean) # calculate the denominator of Moran's I
MI.u <- MI.u.num/MI.u.den # divide the first by the latter to get Moran's I
MI.u

#creating first spatial lag
lag.values <- lag.listw(W.list, shp$values)

#creating Morans'I plot manually
plot(shp$values, lag.values, xlab="Disposable income", ylab="Lag of disp. income")
abline(v=mean(shp$values), lty=2)
abline(h=mean(lag.values, na.rm = T), lty=2)
abline(a=mean(lag.values, na.rm = T)-MI.u*mean(shp$values), b=MI.u)

# Automated Moran's I statistic and plot:
# significance here is evaluated using 
# linear regression based logic and assumptions
moran.test(shp$values, listw = W.list, alternative = "greater", 
           randomisation = FALSE)
# if the parameter randomisation is set to TRUE,
# the variance of I calculated under the assumption
# of randomisation, if FALSE normality
moran.test(shp$values, listw = W.list, alternative = "greater")
moran.plot(shp$values, listw = W.list)

###equal to: 
lm.morantest(lm(shp$values ~ 1), listw = W.list)

# instead, we can also use Monte Carlo simulation,
# which is the preferred method. 
# The way it works that the values are randomly assigned
# to the polygons, and the Moran’s I is 
# computed. This is repeated several times to establish 
# a distribution of expected values. 
# The observed value of Moran’s I is then 
# compared with the simulated distribution to see how likely 
# it is that the observed values could be considered a random draw.
moran.mc(shp$values, listw = W.list, alternative = "greater", nsim = 100)


### LOCAL Moran's I ###

# Calculating Local Moran's I for each spatial unit and evaluating the statistical 
# significance for each I_i The local spatial statistic Moran's I is calculated 
# for each zone based on the spatial weights object used. 
# The values returned include a Z-value, and may be used as a diagnostic tool. 
# The statistic is:

# I_i = \frac{(x_i-\bar{x})}{{∑_{k=1}^{n}(x_k-\bar{x})^2}/(n-1)}{∑_{j=1}^{n}w_{ij}(x_j-\bar{x})}

localI <- localmoran(shp$values,
                     listw = W.list, 
                     alternative = "greater")
# low value indicate no clustering

# NOTE: if we have N = 100 regions and choice p = 0.05 and even if HO was true
# probability of obtaining a false result = 5 %
# We desire a false positive rate of 0 under H0 = TRUE
# p-values apply to I_i --> N independent tests, each having a false positive rate of 0.05
# therefore, the probability of not getting a false postive result is (1-p) and for n obs: (1_p)^n

# what is the probability of getting one or more false postive when looking at all countries denoted by p*
# P* = 1-(1-p)^n
# --> P*  = n*p for small p
# p = 1-(1-p*)^(1/n) or --> p* / n for small p
p.adjust(localI[,5], method = "bonferroni")

localI.data <- data.frame("NUTS_ID" = rownames(localI), 
                          "cluster" = attr(localI, "quadr")[,1], 
                          "p.value" = p.adjust(localI[,5], method = "bonferroni"))

# Plotting maps -----------------------------------------------------------


# some references: 
# https://bookdown.org/nicohahn/making_maps_with_r5/docs/introduction.html
# https://ggplot2-book.org/maps


### base R 

# Plotting k-nearest neighbors
plot(st_geometry(shp), border="grey20")
title("k-nearest Neighb.")
plot(W.list.k$neighbours, 
     coords, add=TRUE, col="red", cex=0.5)

# Easy plotting of ugly looking maps
brks <- quantile(shp$values,
                 probs=seq(0,1,0.25))
cols <- c("yellow", "green", "blue", "red")
plot(st_geometry(shp), border="brown", 
     axes=TRUE, 
     col=cols[findInterval(shp$values, 
                           brks, all.inside = TRUE)], 
     main="Disposable income 2018")
# add legend as before as usual
legend(x = -3555219, y= 10878743, 
       legend = c("Bottom 25%", "25th - Median", 
                  "Median - 75th", "75th - Top"), 
       fill = c("yellow",
                "brown",
                "blue", "red"))


### ggplot2 for bit more advanced map plotting

# may take a look at the ggplot-cheat-sheet

# install.packages("ggplot2")
# install.packages("ggthemes")
# install.packages("viridis")
library(ggplot2)
library(ggthemes) # for map theme used in ggplot (theme_map())
library(viridis) # viridis color scale

p <- ggplot(data = shp) +
  geom_sf(aes(fill=values)) +
  theme_map() +
  labs(x=NULL, y=NULL,
       title="Disposable Income", 
       subtitle="NUTS2", 
       caption = "Source: Eurostat")

p

q <- p + scale_fill_viridis(option="magma",
                            direction = -1)
q

q <- p + 
  theme(legend.position = "bottom") +
  scale_fill_viridis(option = "magma",
                     direction = -1, 
                     name = "Disposable income",
                     guide = guide_colorbar(direction="horizontal", 
                                            barheight = unit(2, units= "mm"),
                                            barwidth = unit(30, units= "mm"),
                                            draw.ulim = F, title.position="top"))

q


# plotting local Moran's I 

shp.moran <- left_join(shp, localI.data)

# plotting all local Moran's I
p_localI <- ggplot(data = shp.moran) +
  geom_sf(aes(fill = cluster)) +
  theme_map() +
  labs(x=NULL, y=NULL,
       title="Local Moran's I for Disposable Income", 
       subtitle="All cluster", 
       caption = "Source: Based on Eurostat")
p_localI

shp.moran <- shp.moran %>% 
  mutate(cluster_sig = replace(cluster, p.value > 0.05, NA))

p_localI_sig <- ggplot(data = shp.moran) +
  geom_sf(aes(fill = cluster_sig)) +
  theme_map() +
  labs(x=NULL, y=NULL,
       title = "Local Moran's I for Disposable Income", 
       subtitle = "Only significant cluster", 
       caption = "Source: Based on Eurostat")
p_localI_sig

