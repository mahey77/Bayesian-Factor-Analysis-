#install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
library(here)
library(rio)
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()) # need this bit to run on multiple cores

##############################################################

stockwatson <- import(here("data","stockwatson.rda"))
View(stockwatson)
# Remove 2008, and the first two quarters w missing data
stockwatson <-  stockwatson[-c(1:2, 197:200),]
# Just pick out first 10 variables:
stockwatson10 <- stockwatson[,1:10]
# How many PCs?
screeplot(prcomp(stockwatson10), npcs=10, type="l", cex.axis = 1.5)
# Looks like k = 3 would be okay

scr <- "dynamic factor model.stan"
scr2 <- "dynamic factor model prediction.stan"

#####################################################################
# Sampler 1
data  <-  list(n = 194,
               p = 10,
               k = 3,
               y = stockwatson10,
               a_idio = 2.2,
               b_idio = sqrt(0.1/2.2),
               C_Lambda = 1)
model_3factors <- stan(file = scr, data = data, iter = 50000, chains=2, seed=784236)
print(model_3factors, pars=c("diag_Sigma", "PLT_Lambda_diag", "PLT_Lambda_lowertri", "Gamma"))

install.packages("plyr")
source("Graphical_Diagnostics.R")
output = as.array(model_3factors, pars=c("diag_Sigma", "PLT_Lambda_diag", "PLT_Lambda_lowertri", "Gamma"))
diagnostics(output, rows=2, lag.max=50)

#############################################################################
data_predictive <-  list(n = 190, # not 194
               p = 10,
               k = 3,
               y = stockwatson10[1:(190),], # just pick out first n - n_star rows
               a_idio = 2.2,
               b_idio = sqrt(0.1/2.2),
               C_Lambda = 1,
               n_star = 4)

#############################################################
# Predicted Value Comparison
model_3factors_predictive <- stan(file = scr2, data = data_predictive, iter = 50000, chains=2, seed=784236)             
y_star_new <-  extract(model_3factors_predictive)$y_star
y_star_new <-  y_star_new[(49000:50000),,]
dim(y_star_new)
y_star_new <-  colMeans(y_star_new, dims=1)
y_star_new <-  as.array(y_star_new[,1:3])
y_star_new
stockwatson_final <- stockwatson10[,1:3]  
stockwatson_final <- stockwatson_final[191:194,]
colnames(stockwatson_final) <-  c("RealGDP", "CPI", "InterestRate")
rownames(stockwatson_final) <-  c("Q1", "Q2","Q3","Q4" )
stockwatson_final <- as.data.frame(stockwatson_final)
library(gt)
library(dplyr)
stockwatson_final|> gt() |>
  tab_stubhead(label = "landmass")
y_star_new <- as.data.frame(y_star_new)
colnames(y_star_new) <-  c("RealGDP", "CPI", "InterestRate")
library(gt)
library(dplyr)
y_star_new|> gt() |>
  tab_stubhead(label = "landmass")
