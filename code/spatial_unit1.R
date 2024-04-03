#------------------------------------->
#- Spatial Economics
#-      Unit 1: Prep Script
#-
#- Original Authors (i.e. many thanks to):
#-  Mathias Moser (matmoser@wu.ac.at)
#-  Franziska Disslbacher
#-    (franziska.disslbacher@wu.ac.at)
#-  Anna Stelzer
#-
#-  Adapted by:
#-  Lukas Vashold
#-    (lukas.vashold@wu.ac.at)
#-  Martin Prinz
#------------------------------------->

#############################################################################################
##################################                         ##################################
##################################   BASE R INTRODUCTION   ##################################
##################################                         ##################################
#############################################################################################



# Working Directory and Workspace -----------------------------------------


#As with the standard R GUI, RStudio employs the notion of a global default working directory. 
#Normally this is the user home directory (typically referenced using ~ in R).

#The current working directory is displayed by RStudio within the title region of the Console.
#Or you can ask R for the current working directory location using the command:

getwd()

# we can create a new folder for the current project we are working on as follows
dir.create("spatial_unit1")


#There are a number of ways to change the current working directory.
#However, the most commonly used one is using the "setwd" function:

setwd("spatial_unit1")

#Or you may use the "Session | Set Working Directory | To Source File Location". This will also 
#change directory location to the folder where your codefile is saved.

#In order to provide a nice folder structure it is reasonable to employ the project with 
#subfolders for your data, images, ...
#To navigate between subfolders, relative or absolute paths can be used.
#It is best practice to have the user running the script begin in a consistent directory 
#on their machine and then use relative file paths from that directory to access files.

# let us create two subfolders in our working directory for R code and data
dir.create("code")
dir.create("data")
# when creating new folders, use lowercase letters in the folder names:
# Windows ignores upper/lowercase but not Mac, which may lead to broken code,
# especially when using the same code on different computers


### Download files from Canvas ###

# download files 'bankstellenverzeichnis.xlsx', 'cornwell.dta', 'GOV_DEBT.csv' and
# 'acs2006.csv' and place it in the folder you created

# you may also want to save this file within the 'code' folder

#Relative paths start  with a "." at your working directory. To navigate into a folder use
#"/Foldername", to get out of a folder use another "."
# "." is the current folder (or working directory)
# ".." is one level above the current folder
#Example:

setwd("./code")

# but for now: let's stay in our unit1 folder as working directory
# ".." changes back in the folder above the current folder
setwd("..")


# let's create a simple scalar object
a <- 5

# save the workspace to the file .RData in the cwd 
save.image()
# we can also specify a certain name for the .RData file
save.image(file="dummy.RData")

# load a workspace into the current session
# if you don't specify the path, the cwd is assumed 
load("dummy.RData")

#q() # quit R. You will be prompted to save the workspace.


#Some other standard commands for managing the workspace are:

# list.files() lists all files that are currently in the working directory
list.files()

# ls() lists all files that are currently in the environment
ls()

# save specific objects (not the whole workspace) to a file
# if you don't specify the path, the wd is assumed 
save(a,file="myfile.RData")

# load a workspace into the current session
# if you don't specify the path, the wd is assumed 
load("myfile.RData")

# while save() and load() always save and restore both an object's representation as well as
# its name - this may be inconvenient when loading an object with the same name as one already 
# in the enviroment; instead, we can use saveRDS() and loadRDS() which only saves the 
# representation of an object, not the name
saveRDS(a, file="a_object.RData")
b <- readRDS("a_object.RData")
# note the difference: when using load("a_object.RData"), we cannot assign a name for the object

# rm() removes elements from the environment
rm(a)
# removes all objects currently in the working directory
rm(list=ls())




# Packages ----------------------------------------------------------------

#Packages are collections of R functions, data, and compiled code in a well-defined format. 
#The directory where packages are stored is called the library. R comes with a standard set 
#of packages. Others are available for download and installation. Once installed, they have 
#to be loaded into the session to be used.

.libPaths() # get library location
library()   # see all packages installed
search()    # see packages currently loaded

#Adding Packages:
#You can expand the types of analyses you do be adding other packages. A complete list of 
#contributed packages is available from CRAN.

#Follow these steps:

#Download and install a package (you only need to do this once).
#To use the package, invoke the library(package) command to load it into the current session. 
#(You need to do this once in each session, unless you customize your environment to 
#automatically load it each time.)


# On MS Windows:
# use the console to install packages, e.g.
# install the package dplyr
# (install only if necessary)
# install.packages("dplyr")
# and load it in the current session
library(dplyr)


#List of commonly used packages:

#To manipulate data:
## dplyr - Essential shortcuts for subsetting, summarizing, rearranging, and joining 
#together data sets. dplyr is our go to package for fast data manipulation.
## tidyr - Tools for changing the layout of your data sets. Use the gather and spread 
#functions to convert your data into the tidy format, the layout R likes best.

#To visualize data:
## ggplot2 - R's famous package for making beautiful graphics. ggplot2 lets you use the 
#grammar of graphics to build layered, customizable plots.

#To model data:
## car - car's Anova function is popular for making type II and type III Anova tables.
## mgcv - Generalized Additive Models
## lme4/nlme - Linear and Non-linear mixed effects models
## multcomp - Tools for multiple comparison testing

#To report results:
## xtable - The xtable function takes an R object (like a data frame) and returns the 
#latex or HTML code you need to paste a pretty version of the object into your documents. 
#Copy and paste, or pair up with R Markdown.

## stargazer - Like xtable, stargazer takes R objects (mainly lists of regression results)
#and returns the latex or HTML code you need to paste a pretty version of the object into 
#your documents.

#For Spatial data:
## sf - Tools for loading and using spatial data including shapefiles
## sp, maptools - (Old) tools for loading and using spatial data including shapefiles.
## maps - Easy to use map polygons for plots.
## ggmap - Download street maps straight from Google maps and use them as a background in 
#your ggplots.

#For Time Series:
## zoo - Provides the most popular format for saving time series objects in R.
## xts - Very flexible tools for manipulating time series data sets.


# Rtools: A collection of tools necessary for building R packages in Windows
# Available for download at http://cran.r-project.org/bin/windows/Rtools/
# Select the .exe download link from the table that corresponds to your version of R
# Note: If you're not sure what version of R you have, open or restart R and it's the 
# first thing that comes up in the console
# Once the download completes, open the .exe file to begin the installation
# Unless you really know what you are doing, you should just go with the default 
# selections at each step of the installation
# Once the Rtools installation completes, (close and re-)open RStudio
# Install the devtools R package if you have not previously done so
# If you aren't sure, enter find.package("devtools") in the console
# To install devtools, use 
# install.packages("devtools")
library(devtools)
# use
find_rtools()
# This should return TRUE in the console if your Rtools installation worked properly

# packages are usually maintained by their creators, this means they might need updating from time
# to time

# list all packages where an update is available
old.packages()

# update all available packages
# update.packages()

# update, without prompts for permission/clarification
# update.packages(ask = FALSE)

# update only a specific package use install.packages()
#install.packages("plotly")



# Help Functions ----------------------------------------------------------

# R includes extensive facilities for accessing documentation and searching for help.
#The "help()" function and "?"-help operator in R provide access to the documentation pages 
#for R functions, data sets, and other objects, both for packages in the standard R 
#distribution and for contributed packages. To access documentation for the standard lm 
#(linear model) function, for example, enter the command: 

help(lm) #or 
help("lm")#, or 
?lm

#You may also use the help() function to access information about a package in your library — 
#for example: 

help(package="MASS")

#which displays an index of available help pages for the package along with other information.

#Help pages for functions usually include a section with executable examples illustrating how 
#the functions work. You can execute these examples in the current R session via the example() 
#command: 

example(lm)

#Many packages include vignettes, which are discursive documents meant to illustrate and 
#explain facilities in the package. You can discover vignettes by accessing the help page 
#for a package, or via the browseVignettes() function. The command: 

browseVignettes()

#opens a list of vignettes from all of your installed packages in your browser.

#The "help()" function and "?" operator are useful only if you already know the name of the 
#function that you wish to use. There are also facilities in the standard R distribution for 
#discovering functions and other objects.

#The "help.search()" function scans the documentation for packages installed in your library. 
#The (first) argument to help.search() is a character string or regular expression. For example: 

help.search("^glm") 

#searches for help pages, vignettes, and code demos.

#It is also possible to use a general search site like Google, by qualifying the search with “R” 
#or the name of an R package (or both). It can be particularly helpful to paste an error message 
#into a search engine to find out whether others have solved a problem that you encountered.

rm(list=ls())






# Basic Commands ----------------------------------------------------------


### R as calculator

# NUMERIC FUNCTIONS

3 + 3         #addition; code works regardless of the spaces, but it is recommended to use them
3 - 5         #subtraction
3 * 5         #multiplication
3 / 5         #division

3 ^ 5         #exponentiation
sqrt(81)      #square root
243 ^ (1/5)   #a-th root
sin(pi/2)     #sine
cos(0)        #cosine
tan(0)        #tangent
log(1)        #natural logarithm
exp(1)        #e^x

ceiling(1.2)  #next higher integer
floor(1.2)    #next lower integer
abs(-1)       #absolute value
round(2.45,1) #rounds given number of digits


# MATRIX ALGEBRA

set.seed(1234) # setting a seed here for reproducibility
m <- round(matrix(sample(1:20, 9), nrow=3),2)
dim(m)
nrow(m); ncol(m)

# All normal operations apply to matrices as well (element-by-element)
m^2
3 * m
sqrt(m)

# Matrix-specific operations

# Matrix multiplication and transpose
m %*% matrix(1:3, ncol=1)
m %*% t(matrix(1:3, nrow=1))

t(m) %*% m

# Inverse ("solving"), determinant, eigenvalues/-vectors
det(m)
solve(m)
eigen(m)

# The trace
sum(diag(m))

# Special matrices
diag(3)
matrix(1, nrow=3, ncol=3)
kronecker(diag(3), m)


# CHARACTER FUNCTIONS

toupper("Test")                     #returns string of uppercase letters
tolower("Test")                     #returns string of lowercase letters
substr("Test", start=2, stop=3)     #returns string from start to stop
strsplit("Test", split="")          #splits string given the split character
paste("T", "e", "s", "t", sep="")   #combines strings with given seperator


# STATISTICAL PROBABILITY FUNCTIONS

rnorm(1,0,1)        #generates 1 random number from the standard normal distribution
dnorm(0)            #normal density function (for standard normal on default)
pnorm(0)            #cumulative normal probability (area ander PDF left from defined value)
qnorm(0.8)          #normal quantile (value at the p percentile of normal distribution)

#r-, d-, p-, q- always evoke the above described commands for a given distributen 
#Other distributions may vary in their parametrization.

#commonly used distributions:

#Normal:    -norm
#Uniform:   -unif
#Beta:      -beta
#gamma:     -gamma
#Binomial:  -binom
#Poisson:   -pois
#Weibull:   -weibull


# OTHER STATISTICAL FUNCTIONS

x <- seq(1:10)

mean(x)	        #arithmetic mean
sd(x)	          #standard deviation
var(x)          #variance
median(x)	      #median
quantile(x)	    #quantiles (quartiles on default)
range(x)	      #range
sum(x)	        #sum
min(x)	        #minimum
max(x)	        #maximum


# OTHER USEFUL FUNCTIONS

# seq(from , to, by)	          #generate a sequence
x <- seq(1:10)
# rep(x, ntimes)	              #repeat x n times
y <- rep(x, 2)
# cut(x, n)	                  #divide continuous variable in factor with n levels 
table(cut(y, 4))



# Data Types --------------------------------------------------------------


# R has a wide variety of data types including scalars, vectors (numerical, character, logical), 
# matrices, data frames, and lists.


### Vectors ###

a <- c(1,2,5.3,6,-2,4) # numeric vector
b <- c("one","two","three") # character vector
c <- c(TRUE,TRUE,TRUE,FALSE,TRUE,FALSE) #logical vector

#Refer to elements of a vector using subscripts.

a[c(2,4)] # 2nd and 4th elements of vector


### Matrices ###

# All columns in a matrix must have the same mode(numeric, character, etc.) and the same length. 
# The general format is:

# mymatrix <- matrix(vector, nrow=r, ncol=c, byrow=FALSE, 
#                    dimnames=list(char_vector_rownames, char_vector_colnames))

# byrow=TRUE indicates that the matrix should be filled by rows. byrow=FALSE indicates that the 
# matrix should be filled by columns (the default). dimnames provides optional labels for the 
# columns and rows.


# generates 5 x 4 numeric matrix 
y<-matrix(1:20, nrow=5,ncol=4)

# another example
cells <- c(1,26,24,68)
rnames <- c("R1", "R2")
cnames <- c("C1", "C2") 
mymatrix <- matrix(cells, nrow=2, ncol=2, byrow=TRUE,
                   dimnames=list(rnames, cnames))

#Identify rows, columns or elements using subscripts.

y[,4] # 4th column of matrix
y[3,] # 3rd row of matrix 
y[2:4,1:3] # rows 2,3,4 of columns 1,2,3


### Arrays ###

#Arrays are similar to matrices but can have more than two dimensions. 
#See help(array) for details.


### Data Frames ###

#A data frame is more general than a matrix, in that different columns can have different modes
#(numeric, character, factor, etc.). This is similar to SAS and SPSS datasets.

d <- c(1,2,3,4)
e <- c("red", "white", "red", NA)
f <- c(TRUE,TRUE,TRUE,FALSE)

myframe <- data.frame(d,e,f)
names(myframe) <- c("ID","Color","Passed") # variable names

#There are a variety of ways to identify the elements of a data frame .

myframe[,1:3] # columns 1, 2, 3 of data frame
myframe[,c("ID","Color")] # columns ID and Age from data frame
myframe$ID # variable x1 in the data frame


### Lists ###

#An ordered collection of objects (components). A list allows you to gather a variety of 
#(possibly unrelated) objects under one name.
# objects in a list can be of different length (e.g. a list of data frames with varying 
# dimensions)

# example of a list with 4 components - 
# a string, a numeric vector, a matrix, and a scaler 
mylist <- list(name="Fred", mynumbers=a, mymatrix=y, age=5.3)

# example of a list containing two lists 
# v <- c(list1,list2)

#Identify elements of a list using the [[]] convention.

mylist[[2]] # 2nd component of the list
mylist[["mynumbers"]] # component named mynumbers in list

# lists names of all elements in a list
names(mylist)

# we can also access multiple items of a list using single brackets 
mylist[c(1:2)]


### Factors ###

#Tell R that a variable is nominal by making it a factor. The factor stores the nominal 
#values as a vector of integers in the range [ 1... k ] (where k is the number of unique 
#values in the nominal variable), and an internal vector of character strings (the original 
#values) mapped to these integers.

# variable gender with 20 "male" entries and 30 "female" entries 
# (the rep command replicates a given argument a given number of times)
# entries are saved as characters
gender <- c(rep("male",20), rep("female", 30)) 
# transform characters into factors
gender <- factor(gender) 

# stores gender as 20 1s and 30 2s and associates
# 1=female, 2=male internally (alphabetically)
# R now treats gender as a nominal variable 
summary(gender)

#R will treat factors as nominal variables in statistical proceedures and graphical analyses. 


### Useful Functions ###

length(a) # number of elements or components
str(mylist)    # structure of an object 
class(y)  # class or type of an object
names(myframe)  # names

c(1,2,3)       # combine objects into a vector
cbind(a, b) # combine objects as columns
rbind(a, b) # combine objects as rows 

mymatrix    # prints the object

which(a>=4) # returns indices in an object which satisfy a given condition



# Importing Data ----------------------------------------------------------

# download GOV_DEBT.csv, cornwell.dta and bankstellenverzeichnis.xlsx from Canvas and
# save it in the /data folder

# Importing data into R is fairly simple. For Stata and Systat, use the foreign package. For 
# SPSS and SAS I would recommend the Hmisc package for ease and functionality.


### From A Comma Delimited Text File (CSV) ###

# CSV (comma-separated Values) can use different delimiters:
# comma (,), semicolon (;), tabs (  ),...
# if the exact data formate is not known, it might be useful to check which delimiter was used
readLines("./data/GOV_DEBT.csv", 3) # shows the first 3 lines of a file
# first row contains variable names, comma is separator 
# note the / instead of \ on mswindows systems 

mydata.csv <- read.csv("./data/GOV_DEBT.csv", sep=",", quote="\"", dec=".", header=TRUE)

#The 2 commands down below can also be used to import CSV-Files.
?read.csv
# csv2 is commonly used in the german speaking world, it is built similarly to csv but is separated by semi-
# colon instead of a comma and uses "," instead of "." to indicate decimals
?read.csv2


### From Excel ###

#One of the best ways to read an Excel file is to export it to a comma delimited file and import 
#it using the method above. Alternatively you can use the xlsx package to access Excel files. 
#The first row should contain variable/column names.

# read in the first worksheet from the workbook myexcel.xlsx
# first row contains variable names
# install.packages("readxl")
library(readxl)
mydata.excel <- read_excel("./data/bankstellenverzeichnis.xlsx")


### From Stata ###

# input Stata file
# install.packages("foreign")
library(foreign)
mydata.stata <- read.dta("./data/cornwell.dta")

# the foreign package can also be used to load SPSS files (read.spss())


### From an URL ###

# if we found a data source online, we can also download it directly into R
# for example, for a csv file, we can use read.csv and include an URL instead of a file path
mydata.url <- read.csv("https://archive.ics.uci.edu/ml/machine-learning-databases/communities/communities.data",
                       header=FALSE)
# using download.file, we can download any file and save it locally - it is, however not yet loaded intot the R
# environment
download.file("http://data.statistik.gv.at/web/meta.jsp?dataset=OGDEXT_NUTS_1_STATISTIK_AUSTRIA_NUTS1_20160101.zip", 
              destfile="./data/Stat_Austria_NUTS2016.zip")


# there are also many data bases which offer packages which help downloading data, e.g.

# World Bank Development Indicators
# install.packages("WDI")
library(WDI)

# Eurostat database
# install.packages("eurostat")
library(eurostat)

# ILO statistics
# install.packages("Rilostat")
library(Rilostat)

# IMF statistics
# install.packages("IMFData")
library(IMFData)


# Exporting Data ----------------------------------------------------------

# There are numerous methods for exporting R objects into other formats . For Stata, you will 
# need to load the foreign packages. For Excel, you will need the xlsx package.


### To an R Data File ###

saveRDS(mydata.excel, file="../data/Bankstellenverzeichnis.rds")


### To A Tab Delimited Text File ###

write.table(mydata.stata, "../data/cornwell.csv", sep=",")




# Functions ---------------------------------------------------------------


# One of the great strengths of R is the user's ability to add functions. In fact, many of the 
# functions in R are actually functions of functions. The structure of a function is given below.

#myfunction <- function(arg1, arg2, ... ){
#statements
#return(object)
#}

# For example:

#setting up a funtion with two arguments:

myfu <- function(arg1,arg2) {  
  A <- paste(arg1,arg2,sep=" ")
  return(A)
}

#running the self written function
myfu(arg1 = "hello", arg2 = "world")



# Handling Data -----------------------------------------------------------


# use the World Bank Development Indicator download option
help("WDI-package")
?WDI

data <- WDI(country="all", indicator=c("SP.POP.TOTL","NY.GDP.MKTP.KD", "SI.POV.GINI"),
            start=2000, end=2016, extra=TRUE)

data <- rename(data, "pop"="SP.POP.TOTL", "gdp"="NY.GDP.MKTP.KD", "gini"="SI.POV.GINI")

# analyse the basic structure of the data
ncol(data)                      # number of columns/variables
nrow(data)                      # number of observations
summary(data) # summary statistics, including no. NAs
head(data) # shows the first rows of some data
tail(data) # shows the last rows of some data
by(data, data$region, summary) # group data by a certain variable (here: region) and apply a function
levels(as.factor(data$country)) # name levels of a factor
table(as.factor(data$region)) # counts the number of observations with per factor level
unique(data$iso2c) # unique values in a vector
duplicated(data$year) # logical indicator whether a values in a vector occurs more than once


### Creating new variables ###

#Use the assignment operator <- to create new variables. 
#Example: we want a varibale gdp per capita

data$gdppc <- data$gdp / data$pop

#Example: Add an ID column:

data$ID <- seq(1,nrow(data))

#providing a new name, will generate a new collumn in the dataframe, while providing an 
#existing name will lead to overwriting.

#Generally, the vector must have the same lenght in order to be uniteable.


### Sorting ###

#To sort a data frame in R, use the order( ) function. By default, sorting is ASCENDING. 
#Prepend the sorting variable by a minus sign to indicate DESCENDING order.

#order itself, just provides a permutation which rearranges its argument into ascending 
#or descending order

order(data$gdppc)       #ascending
order(-data$gdppc)      #descending


# sort data by gdppc (GDP per capita)
newdata <- data[order(data$gdppc),]

#sort by gdppc (ascending) and gini (descending)
newdata <- data[order(data$gdppc, -data$gini),]


### Merging ###

# ADDING COLUMNS

#To merge two data frames (datasets) horizontally, use the merge function. In most cases, 
#you join two data frames by one or more common key variables (i.e., an inner join).

#In order to show this the dataset is first split into 2 sets:
data11 <- data[,c(1:7,ncol(data))]
data12 <- data[,c(8:ncol(data))]

# merge two data frames by ID
datajoin <- merge(data11,data12,by="ID")

# merge two data frames by ID and Country
#total <- merge(data.frame.A,data.frame.B,by=c("ID","Country"))

# if we are sure that the order of the observations are identical in two data frames, we can
# also use cbind() to join columns

# ADDING ROWS
#To join two data frames (datasets) vertically, use the rbind function. The two data frames 
#must have the same variables, but they do not have to be in the same order.

#In order to show this the dataset is first split into 2 sets:
data21 <- data[c(1:2000),]
data22 <- data[c(2001:4488),]

databind <- rbind(data21, data22)

#If a data frame A has variables that data frame B does not, then either:

# Delete the extra variables in data frame A or
# Create the additional variables in data frame B and set them to NA (missing)
# before joining them with rbind( ).

# whenever joining data frames, make sure that variables and observations are in the same 
# order before joining!


### Aggregating ###

#It is relatively easy to collapse data in R using one or more BY variables and 
#a defined function.

# aggregate data frame data by country, returning means for numeric variables

aggdata <-aggregate(data, by=list(data$country), 
                    FUN=mean, na.rm=TRUE)
head(aggdata)
# careful! taking means like this does not always make sense (e.g. for years)
# furthermore, we can only aggregate numerical variables, otherwise R returns NAs

#When using the aggregate() function, the by variables must be in a list (even if there 
#is only one). The function can be built-in or user provided.


### Subsetting ###

#R has powerful indexing features for accessing object elements. These features can be used 
#to select and exclude variables and observations. The following code snippets demonstrate 
#ways to keep or delete variables and observations.

# SELECTING (KEEPING) VARIABLES

# select variables ID, country and gdppc
myvars <- c("ID", "country", "gdppc")
newdata <- data[myvars]

# EXCLUDING (DROPPING) VARIABLES

# exclude variables ID, country and gdppc
myvars <- names(data) %in% c("ID", "country", "gdppc") 
newdata <- data[!myvars]

# SELECTING OBSERVATIONS

# first 5 observations
newdata <- data[1:5,]

# based on variable values
newdata <- data[ which(data$year==2010 & data$gdppc >= 15000), ]

# SUBSET FUNCTION

#The subset( ) function is the easiest way to select variables and observations:

# using subset function to subset for the year 2010
newdata <- subset(data, year==2010)

# use the subset function to subset the data frame for the year 2008 and countries with lower than average GDP!
newdata <- subset(data, year == 2008 | gdp < mean(data$gdp), 
                  select=c(country, gini))
# we can remove NAs by using na.exclude():
newdata <- na.exclude(newdata)



# if Conditions -----------------------------------------------------------

# if:

# An if statement in R consists of three elements:
# The keyword if
# A single logical value between parentheses (or an expression that leads to a single logical value)
# A block of code between braces that has to be executed when the logical value is TRUE

# if (test_expression) {
#   statement
# }

x <- 5
if(x > 0){
  print("Positive number")
}

# or

if ( any(is.na(data$gdp)) ) {
  data <- subset(data, !is.na(data$gdp))
}
# 'data' shrunk from 4488 to 4122 observations

# if - else:

# if (test_expression) {
#   statement1
# } else {
#   statement2
# }

x <- -5
if(x > 0){
  print("Non-negative number")
} else {
  print("Negative number")
}

# if - else statements can be nested in various ways, e.g.

# if ( test_expression1) {
#   statement1
# } else if ( test_expression2) {
#   statement2
# } else if ( test_expression3) {
#   statement3
# } else {
#   statement4
# }

# or

x <- 0
if (x < 0) {
  print("Negative number")
} else if (x > 0) {
  print("Positive number")
} else {
  print("Zero")
}


# ifelse():

# There is a vector equivalent form of the if…else statement in R, the ifelse() function.

data$income.group.pc <- as.factor(ifelse(data$gdppc >= mean(data$gdppc), "high income pc", "low income pc"))

# note: when using ifelse(), make sure that your data does not include any NAs - for any NA, the logical
# statement can neither be confirmed nor denied which results in mulfunctioning of the code
# as an alternative, you could use which():
data$income.group.pc <- NA
data$income.group.pc[which(data$gdppc>=mean(data$gdppc))] <- as.character("high income pc")
data$income.group.pc[which(data$gdppc<mean(data$gdppc))] <- as.character("low income pc")



# Repeating Calculations --------------------------------------------------


#In R you have multiple options when repeating calculations: vectorized operations, loops, 
#and apply functions.


### Looping ###

#The most commonly used loop is the "for" loop. It is used to apply the same function calls 
#to a collection of objects. HIn R, for loops take an interator variable and assign it successive 
#values from a sequence or vector. For loops are most commonly used for iterating over the 
#elements of an object (list, vector, etc.)

for(i in 1:10) {
  print(i)
}

#This loop takes the "i" variable and in each iteration of the loop gives it values 1, 2, 3, …,
#10, executes the code within the curly braces, and then the loop exits.

#Another example:

x <- c("a", "b", "c", "d")

for(j in seq_along(x)){
  print(x[j])
}

#This loop takes the "j" varibale and in each iteration of the loop (based on the length of x),
#it executes the code within the curly braces (print the j-th element of x), and then exits 
#the loop.

#Other typs of loop:

#while-loop:  #While loops begin by testing a condition. If it is true, then they execute the 
#loop body.Once the loop body is executed, the condition is tested again, and 
#so forth, until the condition is false, after which the loop exits.

#repeat-loop: #repeat initiates an infinite loop right from the start. These are not commonly 
#used in statistical or data analysis applications but they do have their uses. 
#The only way to exit a repeat loop is to call break.

# whenever possible:
# try to avoid using loops and opt for vectorized function applications, e.g. using apply:


### Apply Functions ###

#R has some functions which implement looping in a compact form to make your life easier.

#lapply(): Loop over a list and evaluate a function on each element
#sapply(): Same as lapply but try to simplify the result
#apply(): Apply a function over the margins of an array
#tapply(): Apply a function over subsets of a vector
#mapply(): Multivariate version of lapply


#--- LAPPLY [lapply(X, FUN, ...)]

#This function takes three arguments: (1) a list X; (2) a function (or the name of a function) 
#FUN; (3) other arguments via its ... argument.

#Here’s an example of applying the mean() function to all elements of a list. If the original 
#list has names, the the names will be preserved in the output.

x <- list(a = 1:5, b = rnorm(10))

lapply(x, mean)


#--- SAPPLY [sapply(X, FUN, ...)]

#The sapply() function behaves similarly to lapply(); the only real difference is in the return 
#value. sapply() will try to simplify the result of lapply() if possible.

#Notice that lapply() returns a list (as usual), but that each element of the list has length 1.
#Here’s the result of calling sapply() on the same list:

sapply(x, mean)

# the default value for the arguemtn "simplify" in sapply() is TRUE,
# coerces the result to a vector, matric or higher dimensional array if possible
# if we use sapply(..., SIMPLIFY=FALSE), R returns a list
sapply(x, mean, simplify = FALSE)


#--- APPLY [sapply(X, MARGIN, FUN, ...)]

#The apply() function is used to evaluate a function (often an anonymous one) over the margins 
#of an array. It is most often used to apply a function to the rows or columns of a matrix 
#(which is just a 2-dimensional array).


#Here I create a 20x10 matrix of Normal random numbers. I then compute the mean of each column.

x <- matrix(rnorm(200), 20, 10)
apply(x, 2, mean)  # Take the mean of each column


#I can also compute the sum of each row.

apply(x, 1, sum)   ## Take the mean of each row



# Plotting Data -----------------------------------------------------------


# first, tidy up the data before plotting it
plotdata <- data[ which(data$region!="Aggregates" & data$year==2014), ] # exclude aggregates, only 2014
plotdata <- na.exclude(plotdata)

# In R, graphs are typically created interactively.
#Example:

# Creating a Graph
plot(plotdata$gdppc, plotdata$gini) 
abline(lm(plotdata$gini~plotdata$gdppc))
title("Regression of mean Gini on mean GDP p.c.")

#The plot( ) function opens a graph window and plots weight vs. miles per gallon. 
#The next line of code adds a regression line to this graph. The final line adds a title.

# check the plot function for all it's arguments, there are many ways to personalize a plot
help(plot)


### Saving Graphs ###

#You can save the graph via code using one of the following functions:

#pdf("mygraph.pdf")	          pdf file
#win.metafile("mygraph.wmf")	windows metafile
#png("mygraph.png")	          png file
#jpeg("mygraph.jpg")	        jpeg file
#bmp("mygraph.bmp")	          bmp file
#postscript("mygraph.ps")	    postscript file

pdf("Plot.pdf")
plot(plotdata$gdppc, plotdata$gini) 
abline(lm(plotdata$gini~plotdata$gdppc))
title("Regression of mean Gini on mean GDP p.c.")
dev.off()

# in order to save a plot, we need to use the structure above: open an empty pdf file, write
# the plot inside and close it again



### Histogram & Density Plot ###


#--- HISTOGRAM

#You can create histograms with the function hist(x) where x is a numeric vector of values 
#to be plotted. The option freq=FALSE plots probability densities instead of frequencies. 
#The option breaks= controls the number of bins.

# Simple Histogram
hist(plotdata$gdppc)

# Colored Histogram with Different Number of Bins
hist(plotdata$gini, breaks=10, col="red")

# Add a Normal Curve 
x <- plotdata$gini 
h<-hist(x, breaks=10, col="red", xlab="Gini", 
        main="Histogram with Normal Curve") 
xfit<-seq(min(x),max(x),length=40) 
yfit<-dnorm(xfit,mean=mean(x),sd=sd(x)) 
yfit <- yfit*diff(h$mids[1:2])*length(x) 
lines(xfit, yfit, col="blue", lwd=2)


#--- DENSITY PLOT

#Kernal density plots are usually a much more effective way to view the distribution of 
#a variable. Create the plot using plot(density(x)) where x is a numeric vector.

# Kernel Density Plot
d <- density(plotdata$pop) # returns the density data 
plot(d) # plots the results

# Filled Density Plot
plot(d, main="Kernel Density of Population")
polygon(d, col="red", border="blue")



### Dot Plots ###

#Create dotplots with the dotchart(x, labels=) function, where x is a numeric vector and 
#labels is a vector of labels for each point. You can add a groups= option to designate a 
#factor specifying how the elements of x are grouped. If so, the option gcolor= controls 
#the color of the groups label. cex controls the size of the labels.

# Simple Dotplot
dotchart(plotdata$gini,labels=plotdata$country,cex=.7,
         main="Gini coefficients in the world", 
         xlab="Gini")




### Bar Plots ###

#Create barplots with the barplot(height) function, where height is a vector or matrix. 
#If height is a vector, the values determine the heights of the bars in the plot. If height 
#is a matrix and the option beside=FALSE then each bar of the plot corresponds to a column of 
#height, with the values in the column giving the heights of stacked “sub-bars”. If height is 
#a matrix and beside=TRUE, then the values in each column are juxtaposed rather than stacked. 
#Include option names.arg=(character vector) to label the bars. The option horiz=TRUE to 
#create a horizontal barplot.

# Simple Bar Plot 
counts <- table(plotdata$region)
barplot(counts, main="Regional Distribution", xlab="Number of countries in a region")

# Simple Horizontal Bar Plot with Added Labels 
counts <- table(plotdata$region)
barplot(counts, main="Regional Distribution", horiz=TRUE,
        names.arg=levels(data$region))

# if we want to use colours in our plots, the package RColorBrewer offers some nicer choices than 
# R standard
# install.packages("RColorBrewer")
library("RColorBrewer")
display.brewer.all()

# Stacked Bar Plot with Colors and Legend
counts <- table(plotdata$income, plotdata$region)
row.names(counts) <- levels(data$income)
barplot(counts, main="Distribution by region and income",
        xlab="Region", 
        col=brewer.pal(length(levels(data$income)),"Set3"), # selects n=no. income levels colours out of palette Set3
        legend = rownames(counts))

# Grouped Bar Plot
barplot(counts, main="Distribution by region and income",
        xlab="Region", 
        col=brewer.pal(length(levels(data$income)),"Set3"),
        legend = rownames(counts), beside=TRUE)



### Line Charts ###

#Line charts are created with the function lines(x, y, type=) where x and y are numeric 
#vectors of (x,y) points to connect. type= can take the following values:


# p	    #points
# l	    #lines
# o	    #overplotted points and lines
# b, c	#points (empty if "c") joined by lines
# s, S	#stair steps
# h	    #histogram-like vertical lines
# n	    #does not produce any points or lines

#The lines( ) function adds information to a graph. It can not produce a graph on its own. 
#Usually it follows a plot(x, y) command that produces a graph.

#By default, plot( ) plots the (x,y) points. Use the type="n" option in the plot( ) command, 
#to create the graph with axes, titles, etc., but without plotting the points.

#Example
#In the following code each of the type= options is applied to the same dataset. The plot( ) 
#command sets up the graph, but does not plot the points.

x <- plotdata$gdppc; y <- x # specify data
par(pch=22, col="red") # plotting symbol and color 
par(mfrow=c(2,4)) # all plots on one page 
opts = c("p","l","o","b","c","s","S","h") 
for(i in 1:length(opts)){ 
  heading = paste("type=",opts[i]) 
  plot(x, y, type="n", main=heading) 
  lines(x, y, type=opts[i]) 
}

#Next, we demonstrate each of the type= options when plot( ) sets up the graph and does plot 
#the points.

x <- plotdata$gdppc; y <- x # specify data
par(pch=22, col="blue") # plotting symbol and color 
par(mfrow=c(2,4)) # all plots on one page 
opts = c("p","l","o","b","c","s","S","h") 
for(i in 1:length(opts)){ 
  heading = paste("type=",opts[i]) 
  plot(x, y, main=heading) 
  lines(x, y, type=opts[i]) 
}

#As you can see, the type="c" option only looks different from the type="b" option 
#if the plotting of points is suppressed in the plot( ) command.

dev.off() #to reset the plot options



### Boxplots ###

#Boxplots can be created for individual variables or for variables by group. The format 
#is boxplot(x, data=), where x is a formula and data= denotes the data frame providing the 
#data. An example of a formula is y~group where a separate boxplot for numeric variable y is 
#generated for each value of group. Add varwidth=TRUE to make boxplot widths proportional to 
#the square root of the samples sizes. Add horizontal=TRUE to reverse the axis orientation.


# Boxplot of GDP by region
boxplot(gdp~region,data=plotdata, main="World Income Data", 
        xlab="GDP per capita", ylab="Gini")



### Scatterplots ###

#There are many ways to create a scatterplot in R. The basic function is plot(x, y), where x 
#and y are numeric vectors denoting the (x,y) points to plot.

# Simple Scatterplot
plot(data$gdppc, data$gini, main="World Income Data", 
     xlab="GDP per capita", ylab="Gini", pch=19)

# Add fit lines
abline(lm(gini~gdppc,data=plotdata), col="red") # regression line (y~x) 

#--- Scatterplot Matrices


# Basic Scatterplot Matrix
pairs(~gdppc+factor(region)+gini,data=plotdata, 
      main="Simple Scatterplot Matrix")



# DPLYR -------------------------------------------------------------------


#To install dplyr
#install.packages("dplyr")

#To load dplyr
library(dplyr)

# download and save the folder "ACS 2006" in the "./data" folder
# ACS: American Community Survey, 2006
# the file below contains household data, 691333 observations and 188 variables
ACS06 <-  read.csv("./data/acs2006.csv")


#dplyr verbs	    Description
#select()	        select columns
#filter()	        filter rows
#arrange()	      re-order or arrange rows
#mutate()	        create new columns
#summarise()	    summarise values
#group_by()	      allows for group operations in the 'split-apply-combine' concept
#Select a set of columns: 
ACS06 <- select(ACS06,
                RT, SERIALNO, PUMA, ST, ADJUST, WGTP, NP, ACR, BDS,
                CONP, ELEP, GASP, FINSP, MHP, FULP, MRGP, RMS, SMP, WATP, FES, FINCP,
                RNTP, GRNTP, GRPIP, HHL, HHT, HINCP, LNGI, SMOCP, TAXP, MV, PLM)
head(ACS06)

# when using base R, we can use str() to get information about the data structure
str(ACS06)
# the dplyr option is called glimpse(), which has better formatting, and adapts to your screen width
glimpse(ACS06)

#To select all the columns except a specific column, use the '-' (subtraction) operator:
head(select(ACS06, -RT))

#To select a range of columns by name, use the ':' (colon) operator:
head(select(ACS06, PUMA:HINCP))

#To select all columns that start with the character string 'sl', use the function starts_with():
head(select(ACS06, starts_with("f")))

#ends_with() = Select columns that end with a character string
#contains() = Select columns that contain a character string
#matches() = Select columns that match a regular expression
#one_of() = Select columns names that are from a group of names

# Filtering

# calculate mean HH income
mean.income <- mean(ACS06$HINCP, na.rm=TRUE)

# filter only by one condition
filter(ACS06, HINCP >= mean.income)
# in base R, we would write the following to get the same result
# ACS06[ACS06$HINCP >= mean.income, ]

# filter by multiple conditions
filter(ACS06,  HINCP < mean.income, RMS <= 4)

#Filter the rows for PUMA 700 or 100
filter(ACS06, PUMA %in% c(700, 100))

# recode(): replace numeric values based on their position, and character values by their name
# let's use it to recode the variable PLM (complete plumbing facilities)
ACS06$PLM <- recode(ACS06$PLM, '1'="yes", '2'="no")
# note that we now created a character column, if we need factors, we use recode_factor()
# recode_factor() works like recode but creates factors with levels, we can use it to recode the variable LNGI
# (Linguistic isolation), where 1 indicates "Not linguistically isolated" and 2 
# "Linguistically isolated"
ACS06$LNGI <- recode_factor(ACS06$LNGI, '1' = "not isolated", '2' = "isolated")

# rename(): renames variable names
# this code renames the variable ST as state and WGTP as weight
head(rename(ACS06, state=ST, weight=WGTP))

# setNames(): sets the names on an object and returns the object
head(setNames(select(ACS06, PUMA:BDS), c("Puma", "state", "adjust.dollar", "weight",
                                         "persons", "lotsize", "bedrooms")))

# slice(): select rows by position
slice(ACS06, c(1,3,6))
head(ACS06)

# join(): combines two data frames by joining columns, there are several methods, check help(join)
# for all of them; the most useful usually are left_join() and right_join()
# both need a unique identifier which allows to match data from one to the other data frame, in 
# our case: SERIALNO which identifies each household
ACS06.x <- select(ACS06, colnames(ACS06)[1:15])
# let's only select a subset of the observations of the remaining variables
ACS06.y1 <- select(ACS06, colnames(ACS06)[c(2, 16:32)])
ACS06.y2 <- filter(ACS06.y1, HINCP >0)

# left_join(x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"), ...):
# return all rows from x, and all columns from x and y. Rows in x with no match in y will have NA 
# values in the new columns. If there are multiple matches between x and y, all combinations of the 
# matches are returned.
nrow(left_join(ACS06.x, ACS06.y2, by="SERIALNO")) # keeps all observations from ACS06.x

# right_join(x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"), ...):
# return all rows from y, and all columns from x and y. Rows in y with no match in x will have NA 
# values in the new columns. If there are multiple matches between x and y, all combinations of the 
# matches are returned.
nrow(right_join(ACS06.x, ACS06.y2, by="SERIALNO")) # keeps only observations from ACS06.y


# bind_rows()/bind_cols(): work similar to rbind() and cbind(), but have some advantages, e.g.
# bind_rows allow to combine data frames with different dimensions, inserting NAs when data is missing
# also, it can be provided with a list of data frames as arguments, which can be very useful
columns <- as.list(ACS06)

# let's put them back together and assign the original column names
head(bind_cols(columns))

head(bind_cols(ACS06.x, ACS06.y1[,-1]))

# Replacing NAs
sum(is.na(ACS06$HINCP)) # there are 90906 observations with missing HINCP
sum(is.na(coalesce(ACS06$HINCP, 0L))) # no more NAs

# case_when works as a vectorised if-function, which can also combine several if statements
# the following creates a new variable which includes factor "high" for household income greater
# than the mean HH income and "low" for HH income at most as high as the  mean HH income
ACS06$Inc <- as.factor(case_when(ACS06$HINCP > mean.income ~ "high",
                                 ACS06$HINCP <= mean.income ~ "low"))

# Pipe operator: %>% 

#Here's an example you have seen:
head(select(ACS06, PUMA, HINCP))

#Now in this case, we will pipe the ACS data frame to the function that will select two columns 
#(PUMA and HINCP) and then pipe the new data frame to the function head() which will return 
#the head of the new data frame.
ACS06 %>% select(PUMA, HINCP) %>% head

#You might want to write functions in separated rows to further ease the legibility of your code:
ACS06 %>% 
  select(PUMA, HINCP) %>% 
  head

# You will soon see how useful the pipe operator is when we start to combine many functions.

# Sorting

#To arrange (or re-order) rows by a particular column such as household incomee, list the name 
#of the column you want to arrange the rows by:
ACS06 %>% 
  arrange(HINCP) %>% 
  head

#Now, we will select three columns from ACS06, arrange the rows by the household income and then 
#arrange the rows by the number of rooms. Finally show the head of the final data frame:
ACS06 %>% 
  select(SERIALNO, HINCP, RMS) %>%
  arrange(HINCP, RMS) %>% 
  head

#Same as above, except here we filter the rows for households with more than 6 rooms instead 
#of showing the head of the final data frame:
ACS06 %>% 
  select(SERIALNO, HINCP, RMS) %>%
  arrange(HINCP, RMS) %>% 
  filter(RMS > 6) %>% 
  head

#Something slightly more complicated: same as above, except arrange the rows in the HINCP 
#column in a descending order. For this, use the function desc():
ACS06 %>% 
  select(SERIALNO, HINCP, RMS) %>%
  arrange(desc(HINCP), RMS) %>% 
  filter(RMS > 6) %>% 
  head


# Editing variables

#The mutate() function will add new columns to the data frame. Create a new column called 
#RMSPP (rooms per person)
ACS06 %>% 
  mutate(RMSPP = RMS / NP) %>%
  head

#You can generate many new columns using mutate (separated by commas). Here we add a second column 
#called ELEPINC which is the share of monthly electricity costs in household income:
ACS06 %>% 
  mutate(RMSPP = RMS / NP,
         ELEPINC = ELEP / HINCP) %>%
  head

#if you only want to keep some old and the new variable, use transmute()
ACS06 %>% 
  transmute(SERIALNO, 
            RMSPP = RMS / NP,
            ELEPINC = ELEP / HINCP) %>%
  head

# Aggregating

#The summarise() function will create summary statistics for a given column in the data frame 
#such as finding the mean. For example, to compute the average number of people in a household, apply 
#the mean() function to the column sleep_total and call the summary value avg_NP:
ACS06 %>% 
  summarise(avg_NP = mean(NP))

#There are many other summary statistics you could consider such sd(), min(), max(), median(), 
#sum(), n() (returns the length of vector), first() (returns first value in vector), last() 
#(returns last value in vector) and n_distinct() (number of distinct values in vector).

#EXAMPLE:
ACS06 %>% 
  summarise(avg_NP = mean(NP), 
            min_NP = min(NP),
            max_NP = max(NP),
            total = n())
ACS06 %>% 
  select(SERIALNO, HINCP, RMS, NP) %>%
  filter(NP!=0) %>%
  mutate(RMSPP = RMS / NP) %>%
  arrange(HINCP, RMS) %>% 
  filter(RMS > 6) %>%
  summarise(avg_NP = mean(NP), 
            avg_RMSPP = mean(RMSPP),
            total=n())

# summarise_each allows you to apply the same summary function to multiple columns at once

# Grouped Operations

#The group_by() verb is an important function in dplyr. As mentioned before it's related to 
#the concept of 'split-apply-combine'. We literally want to split the data frame by some variable 
#(e.g. HHT, household type), apply a function to the individual data frames and then combine the output.

#Let's do that: split the msleep data frame by the household type, then ask for the same summary 
#statistics as above. We expect a set of summary statistics for each household type.

ACS06 %>% 
  group_by(HHT) %>%
  summarise(avg_HINCP = mean(HINCP), 
            min_NP = min(NP), 
            max_NP= max(NP),
            total = n())

#All datasets and summary statistics generated above can easily be saved by using the "<-" 
#operator on the whole piping command:

ACS06_stats <- ACS06 %>% 
  group_by(HHT) %>%
  summarise(avg_HINCP = mean(HINCP), 
            min_NP = min(NP), 
            max_NP= max(NP),
            total = n())

# Aggregate on multiple columns

#There are scoped variants of summarise(), mutate() and transmute(), they apply operations on a 
#selection of variables.
#summarise_all(), mutate_all() and transmute_all() apply the functions to all (non-grouping) columns.
#summarise_at(), mutate_at() and transmute_at() allow you to select columns using the same name-based 
#select_helpers just like with select().
#summarise_if(), mutate_if() and transmute_if() operate on columns for which a predicate returns TRUE.

# _all: applies function to all variables
ACS06 %>%
  group_by(HHT) %>%
  summarise_all(mean)

# _at: applies function to specified variables
ACS06 %>% 
  group_by(HHT) %>%
  summarise_at(c("HINCP", "NP"), mean, na.rm = TRUE)

ACS06 %>% 
  mutate_at(vars(matches("HINCP")), log)

# _if: applies function if a condition is met
ACS06 %>% 
  summarise_if(is.numeric, mean, na.rm = TRUE)

ACS06%>%
  mutate_if(is.factor, as.character) %>%
  glimpse


# Estimation --------------------------------------------------------------


# Estimating OLS
library(MASS)
data(Boston)
names(Boston)

## Simple Regression
# Seek to predict
# -> medv (median house value)
# by
# -> lstat (percent of households with low socioeconomic status)

# scatterplot
plot(medv~lstat, data=Boston)

# fit linear model by least squares
lm.fit <- lm(medv~lstat, data=Boston)

# basic information about the model
lm.fit

# add regression line to scatterplot
abline(lm.fit, col='red')

# more detailed output: std. errors, p-values, R2..
summary(lm.fit)

# what other information is stored in lm.fit?
names(lm.fit)

# parameter estimates
lm.fit$coefficients
coef(lm.fit)

# confidence interval for the coefficient estimates
confint(lm.fit)

# prediction of medv for a given value of lstat
predict(lm.fit, data.frame(lstat=(c(5,10,15))), interval="confidence")
predict(lm.fit, data.frame(lstat=(c(5,10,15))), interval="prediction")

# replicate the above manually
x <- Boston$lstat
y <- Boston$medv
n <- length(y)

# calculate the estimators
b1 <- sum((x - mean(x)) * (y - mean(y))) / sum((x - mean(x))^2)
b0 <- mean(y) - b1 * mean(x)

# generate fitted values, residuals
y.fit <- b0 + b1 * x
e <- y - y.fit
s2 <- var(e)

# calculate the standard errors of b0, b1
se.b0 <- sqrt(s2*(1/n + mean(x)^2/sum((x-mean(x))^2)))
se.b1 <- sqrt(s2/sum((x-mean(x))^2))

# calculate t-statistic, check significance
t.stat <- c((b0-0)/se.b0, (b1-0)/se.b1)
pt(abs(t.stat), n-2, lower.tail=FALSE)

# calculate the confidence intervals using the t-distribution
th <- qt(0.975, n-2)
conf.b0 <- c(b0-th*se.b0, b0+th*se.b0)
conf.b1 <- c(b1-th*se.b1, b1+th*se.b1)

# predict
y.hat <- b0 + b1*c(5, 10, 15)

# model accuracy
rss <- sum((y-y.fit)^2)
tss <- sum((y-mean(y))^2)
r2 <- 1-(rss/tss)

# first spatial view: simple dummy variable

## "chas" is a dummy variable where 1 indicates the regiona round the Charles river
plot(medv~chas, Boston)
lm2 <- lm(medv ~ lstat + as.factor(chas), data=Boston)


# a simple gravity model:
# "dis" is distance to an employment centre
plot(medv~dis, Boston)
lm3 <- lm(medv ~ lstat + dis, data=Boston, weights = Boston$zn/sum(Boston$zn))

# Multiple OLS by hand
X <- Boston[,c("crim", "indus", "age")]
class(X)

X <- as.matrix(X)
class(X)
X <- cbind(intercept=1, X)

solve(t(X)%*%X) %*% t(X)%*%y

# weighted OLS
solve(t(X)%*%diag(Boston$zn)%*%X) %*% (t(X)%*%diag(Boston$zn)%*%y)


### Estimating Maximum Likelihood ###

logL <- function(par, XX=X, yy=y) {
  sigma <- par[length(par)]
  beta <- par[-length(par)]
  -sum(-0.5*log(sigma) - 0.5*log(2*pi) - (1/(2*sigma)) * (yy-XX%*%matrix(beta))^2)
}

optim(rep(0.1, 5), logL, control=list(maxit=100000, abstol=10^-6))

# Weighted ML
W <- Boston$zn/sum(Boston$zn)

logLw <- function(par, XX=X, yy=y, ww=W) {
  sigma <- par[length(par)]
  beta <- par[-length(par)]
  -sum( ww*(-0.5*log(sigma) - 0.5*log(2*pi) - (1/(2*sigma)) * (yy-XX%*%matrix(beta))^2))
}

optim(rep(0.1, 5), logLw, method="BFGS", control=list(maxit=100000, abstol=10^-6))$par



logL2 <- function(par, XX=X, yy=y) {
  sigma <- par[length(par)]
  beta <- par[-length(par)]
  -sum(log(dnorm(yy, XX%*%matrix(beta), sigma)))
}
optim(c(rep(1, 4), 10), logL2, control=list(maxit=100000, abstol=10^-6))


