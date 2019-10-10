library(secure)
library(MASS)
require(ggplot2)
require(reshape2)
library(splines)
source("R/rvsmanopt.R")
set.seed(1)

data(CellCycle)
# Imputation By Random Forests
# library(imputeMissings)
# Y <- data.frame(CellCycle$Y)
# imp.rf.Y <- imputeMissings::impute(Y, method="randomForest")
# imp.Y <- as.matrix(imp.rf.Y)
# save(imp.Y, file="application/data/imputed_CellCycle_Y.RData")

# load dataset ------------------------------------------------------------
load("application/data/imputed_CellCycle_Y.RData")
X <- CellCycle$X
Y <- imp.Y
n <- nrow(Y)
p <- ncol(X)
q <-  ncol(Y)

# Analysis of Experimentary Comfirmed TFs  --------------------------------

# fit for CellCycle Dataset
# SeCURE
max.est.rank <- min(n, p, q)
control <- secure.control(elnetAlpha = 1, spU=1.0,spV=1.0)
fit.secure.al <- secure.path(Y, X,nrank = max.est.rank, nlambda = 100, orthXU = T, orthV = T, control = control)

control <- secure.control(spU=1.0,spV=1.0)
fit.secure.ae <- secure.path(Y, X,nrank = max.est.rank, nlambda = 100, orthXU = T, orthV = T, control = control)

# RVSManOpt
parameters <- rvs.manopt.parameters(alpha=0.90, U.gamma=2, V.gamma=2, D.gamma=0.5)
fit.rvs.manopt <- rvs.manopt.path(X, Y, nrank = max.est.rank, parameters)

cat("estimated rank: ", dim(fit.rvs.manopt$U)[2])

# Comparison with confirmed genes
genes.confirmed <- c("ACE2","SWI4","SWI5","SWI6","MBP1","STB1","FKH1","FKH2",
                     "NDD1","MCM1","ABF1","BAS1","CBF1","GCN4","GCR1","GCR2",
                     "LEU3","MET31","REB1","SKN7","STE12")
genes.identify.secure.al <- colnames(X)[apply(fit.secure.al$U,1,sum)!=0]
genes.identify.secure.ae <- colnames(X)[apply(fit.secure.ae$U,1,sum)!=0]
genes.identify.rvs.manopt <- colnames(X)[apply(fit.rvs.manopt$U,1,sum)!=0]

genes.identify.secure.al.intersect <- intersect(genes.identify.secure.al, genes.confirmed)
genes.identify.secure.ae.intersect <- intersect(genes.identify.secure.ae, genes.confirmed)
genes.identify.rvs.manopt.intersect <- intersect(genes.identify.rvs.manopt, genes.confirmed)

# the number of confirmed genes selected.
num.confirmed.genes.identify.secure.al <- length(genes.identify.secure.al.intersect) 
num.confirmed.genes.identify.secure.ae <- length(genes.identify.secure.ae.intersect) 
num.confirmed.genes.identify.rvs.manopt <- length(genes.identify.rvs.manopt.intersect) 

# the number of overall genes selected.
num.all.genes.identify.secure.al <- sum(apply(fit.secure.al$U,1,sum)!=0)  
num.all.genes.identify.secure.ae <- sum(apply(fit.secure.ae$U,1,sum)!=0)  
num.all.genes.identify.rvs.manopt <- sum(apply(fit.rvs.manopt$U,1,sum)!=0)  

# Selection rate
rate.secure.al <- num.confirmed.genes.identify.secure.al/num.all.genes.identify.secure.al
rate.secure.ae <- num.confirmed.genes.identify.secure.ae/num.all.genes.identify.secure.ae
rate.rvs.manopt <- num.confirmed.genes.identify.rvs.manopt/num.all.genes.identify.rvs.manopt

cat("secure(AL) rate: ", rate.secure.al, "num.all.genes.identify.secure.al: ", num.all.genes.identify.secure.al, "num.confirmed.genes.identify.secure.al", num.confirmed.genes.identify.secure.al, "\n")
cat("secure(AE) rate: ", rate.secure.ae, "num.all.genes.identify.secure.ae: ", num.all.genes.identify.secure.ae, "num.confirmed.genes.identify.secure.ae", num.confirmed.genes.identify.secure.ae, "\n")
cat("rvs.manopt rate: ", rate.rvs.manopt, "num.all.genes.identify.rvs.manopt: ", num.all.genes.identify.rvs.manopt, "num.confirmed.genes.identify.rvs.manopt", num.confirmed.genes.identify.rvs.manopt,"\n")


# Plot of all confirmed TF effect ---------------------------------------------
# for(i in 1:length(genes.confirmed)){
#   filename <- paste0("results/yeast/", genes.confirmed[i], ".pdf")
#   pdf(filename, width = 10, height = 10)
#   plot(seq(0,119,7),fit.rvs.manopt$C.est[which(colnames(X) == genes.confirmed[i]),],type="l",
#        xlab="",ylab="",cex.axis = 2)
#   abline(h=0,col="grey")
#   dev.off()
# }

