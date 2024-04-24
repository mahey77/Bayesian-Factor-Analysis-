

######################
library(psych)
library(dplyr)


stockwatson <- import(here("data","stockwatson.rda"))
stockwatson <-  stockwatson[-c(1:2, 197:200),]
# Just pick out first 10 variables:
stockwatson10 <- stockwatson[,1:10]

prcomp(stockwatson10)
sd <- apply(stockwatson10, 2, sd)
sd
S <- cor(stockwatson10)
eigdec <- eigen(S)
eigdec
eig <- eigdec$values
eig
eigvec <-  eigdec$values
sum(eig)
# sum(eig)  = 40 as expected
######################################################
pc_loading <- eigdec$vectors
rownames(pc_loading) <- colnames(StockandWatson)
pc_loading

###################################################
# Variances in percentage
eig <- eigdec$values
variance <- eig*100/sum(eig)
# Cumulative variances
cumvar <- cumsum(variance)
eig2 <- data.frame(Eigenvalue = eig, Variance = variance,
                   Cumulative_Variance = cumvar)
eig2
colnames(eig2) <-  c("Eigenvalue","Variance","Cumulative Variance")
eig2 <- data.frame(eig2)
install.packages("gt")
library(gt)
library(dplyr)
eig2 |> gt() |>
  tab_stubhead(label = "landmass")


###############################################
barplot(eig2[, 2], names.arg=1:nrow(eig2),
        xlab = "Number of Principal Components",
        ylim = c(0,max(eig2[,2])*1.2) ,
        ylab = "Percentage of Overall Variance")
        
# Add connected line segments to the plot
lines(x = 1:nrow(eig2), eig2[, 2],
      type="b", pch=18, col = "black")


