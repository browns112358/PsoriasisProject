library(GEOquery)
library(glmnet)
library(lasso2)
library(illuminaHumanv4.db)
library(matrixStats)

rm(list=ls())
cat("\014")

#following the R script generated by GEO2R on the GSE47598 page
#fetch the data matrix from GEO.  This contains a bunch of stuff including the expression for each gene/individual.
gse <- getGEO("GSE47598", GSEMatrix = TRUE) 

#I believe this is formatting to make the gse list convertible later with the exprs function
if (length(gse) > 1) idx <- grep("GPL10558", attr(gse, "names")) else idx <- 1
gse <- gse[[idx]]

#show(gse)
#This is the classification vector (0 means patient, 1 is control)
y =     c(0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,1,0,1,1,1)

#Convert the GSE list into a matrix (and traspose so rows are individuals (sample) and columns are genes (features).
#eset, which is the matrix we are interested in in huge so viewing it all takes forever in R.  This just
# let me verify that the data matches the GSM files
eset <- t(exprs(gse))
show(eset[,1:5])

rm(gse)

#now do lasso on data set
grid =10^seq(2,-8, length=100)

lasso.geo = glmnet(as.matrix(eset), y, alpha=1, family="binomial", lambda=grid, standardize = TRUE)

plot(lasso.geo, xvar = 'lambda')
show(summary(lasso.geo))

which(coef(lasso.geo)[,100] != 0)
ind100 = which(coef(lasso.geo)[,100] != 0)
coefs100.lasso = coef(lasso.geo)[ind100[-1],100]
plot(coefs100.lasso)

which(coef(lasso.geo)[,50] != 0)
ind50 = which(coef(lasso.geo)[,50] != 0)
coefs50.lasso = coef(lasso.geo)[ind50[-1],50]
plot(coefs50.lasso)


ridge.geo = glmnet(as.matrix(eset), y, alpha=0, family="binomial", lambda=grid, standardize = TRUE)
plot(ridge.geo, xvar = 'lambda')