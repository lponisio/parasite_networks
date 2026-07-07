rm(list=ls())
source("lab_paths.R")
setwd(local.path)
setwd("parasite_networks")

library(corrr)
library(ggcorrplot)
library(FactoMineR)
library(factoextra)
library(ggplot2)

load(file="../parasite_networks/data/sp_mets.RData")

# Folder for PCA PDFs
pca.fig.dir <- file.path("figures", "pcas")
dir.create(pca.fig.dir, recursive = TRUE, showWarnings = FALSE)

net.cols <- c("zweighted.betweenness",
              "zweighted.closeness",
              "zd",
              "zdegree",
              "zHBOverlap",
              "zproportional.generality")

fin.cols <- !is.na(sp.network.metrics[, net.cols])
keep.rows <- apply(fin.cols, 1, function(x) all(x))
sp.network.metrics <- sp.network.metrics[keep.rows,]

sp.network.metrics[, net.cols] <-
  apply(sp.network.metrics[, net.cols], 2, scale)

net.pca <- princomp(sp.network.metrics[, net.cols])
summary(net.pca)

# Eigenvalue / scree plot
p.eig <- fviz_eig(net.pca, addlabels = TRUE)

ggsave(
  filename = file.path(pca.fig.dir, "pca_eigenvalues.pdf"),
  plot = p.eig,
  width = 7,
  height = 5,
  units = "in"
)

# PCA variable contribution plot
p.var <- fviz_pca_var(net.pca, col.var = "black")

ggsave(
  filename = file.path(pca.fig.dir, "pca_variables.pdf"),
  plot = p.var,
  width = 7,
  height = 5,
  units = "in"
)

# Cos2 plot
p.cos2 <- fviz_cos2(net.pca, choice = "var", axes = 1:3)

ggsave(
  filename = file.path(pca.fig.dir, "pca_variable_cos2_axes_1_3.pdf"),
  plot = p.cos2,
  width = 7,
  height = 5,
  units = "in"
)

# Variance explained by each PC
pca.var <- net.pca$sdev^2

# Proportion of variance explained by each PC
pca.pvar <- pca.var / sum(pca.var)

# Cumulative variance explained plot
pdf(
  file = file.path(pca.fig.dir, "pca_cumulative_variance_explained.pdf"),
  width = 7,
  height = 5
)

plot(
  cumsum(pca.pvar),
  xlab = "Principal component",
  ylab = "Cumulative proportion of variance explained",
  ylim = c(0, 1),
  type = "b"
)
grid()
abline(h = 0.95, col = "blue")

dev.off()
