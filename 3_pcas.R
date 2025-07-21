
rm(list=ls())
setwd("~/University of Oregon Dropbox/Lauren Ponisio/")
setwd("parasite_networks")
library(corrr)
library(ggcorrplot)
library(FactoMineR)
library(factoextra)

load(file="../parasite_networks/data/sp_mets.RData")

net.cols <- c("zweighted.betweenness",
              "zweighted.closeness",
              "zd",
              "zdegree",
              "zHBOverlap")

fin.cols <- !is.na(sp.network.metrics[, net.cols])
keep.rows <- apply(fin.cols, 1, function(x) all(x))
sp.network.metrics <- sp.network.metrics[keep.rows,]

sp.network.metrics[, net.cols] <-
  apply(sp.network.metrics[, net.cols], 2, scale)

net.pca <- princomp(sp.network.metrics[, net.cols])
summary(net.pca)

fviz_eig(net.pca, addlabels = TRUE)

fviz_pca_var(net.pca, col.var = "black")

fviz_cos2(net.pca, choice = "var", axes = 1:3)


# Variance explained by each PC
pca.var <- net.pca$sdev^2  

# Proportion of variance explained by each PC
pca.pvar <- pca.var/sum(pca.var) 

# Cumulative variance explained plot
plot(cumsum(pca.pvar),
     xlab = "Principal component",
     ylab = "Cumulative Proportion of variance explained",
     ylim = c(0,1), type = 'b')
grid()

# Add a horizontal line
abline(h = 0.95, col = "blue")
