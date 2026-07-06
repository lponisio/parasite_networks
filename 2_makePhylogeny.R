## Generate phylogenies and plot them
rm(list=ls())
source("lab_paths.R")
setwd(local.path)
setwd("parasite_networks")

library(ape)

## load tree from:
## Henriquez Piskulich, Patricia Andrea; Hugall, Andrew F.; Stuart-Fox, Devi (2023).
## A supermatrix phylogeny of the world’s bees (Hymenoptera: Anthophila) [Dataset].
## Dryad. https://doi.org/10.5061/dryad.80gb5mkw1

load(file="data/sp_mets.RData")

## Folder for phylogeny plots
phylo.fig.dir <- file.path("figures", "phylogenies")
dir.create(phylo.fig.dir, recursive = TRUE, showWarnings = FALSE)

## Helper to save a phylogeny plot as a PDF
plotSaveTree <- function(tree,
                         filename,
                         main = NULL,
                         type = "phylogram",
                         show.tip.label = TRUE,
                         cex = 0.7,
                         width = 8,
                         height = 10) {
  pdf(file = filename, width = width, height = height)
  par(mar = c(2, 2, 3, 2))
  ape::plot.phylo(
    tree,
    type = type,
    show.tip.label = show.tip.label,
    cex = cex,
    no.margin = FALSE,
    main = main
  )
  dev.off()
}

phylo <- ape::read.tree("data/BEE_mat7_fulltree.nwk")

## Clean up unwanted portion of labels
pattern <- "(_n\\d+m\\d+_[A-Za-z0-9]+)?$"
phylo$tip.label <- gsub(pattern, "", phylo$tip.label)

## Replace underscore with space
phylo$tip.label <- gsub("_", " ", phylo$tip.label)

## ## Plot the full bee phylogeny.
## ## Tip labels are suppressed because the full tree is very large.
## plotSaveTree(
##   tree = phylo,
##   filename = file.path(phylo.fig.dir, "full_bee_phylogeny.pdf"),
##   main = "Full bee phylogeny",
##   type = "fan",
##   show.tip.label = FALSE,
##   width = 12,
##   height = 12
## )

makeSaveTree <- function(spec.dats, phylo, genus){
  genus_name <- paste(genus, collapse = "_")

  spec.dats <- spec.dats[
    spec.dats$Genus %in% genus,]

  species_to_keep <- na.omit(unique(spec.dats$GenusSpecies))
  print(paste(genus_name, "species in dataset:", length(species_to_keep)))

  phylo_tips <- phylo$tip.label
  print(paste("tips in full phylogeny:", length(phylo_tips)))

  species_in_tree <- species_to_keep[species_to_keep %in% phylo_tips]
  missing_species <- species_to_keep[!species_to_keep %in% phylo_tips]

  ## Only keep tips that match our species
  phylo <- ape::keep.tip(phylo, species_in_tree)
  print(paste(genus_name, "tips retained:", length(phylo$tip.label)))
  print(paste("missing", paste(missing_species, collapse = ", ")))

  ## Calculate covariance matrix
  phylo_matrix <- ape::vcv.phylo(phylo)

  ## Save pruned tree and covariance matrix
  save(
    phylo,
    phylo_matrix,
    file = sprintf("data/%s_phylogeny.Rdata", genus_name)
  )

  ## Plot pruned phylogeny
  plotSaveTree(
    tree = phylo,
    filename = file.path(phylo.fig.dir,
                         sprintf("%s_phylogeny.pdf", genus_name)),
    main = sprintf("%s phylogeny", genus_name),
    type = "phylogram",
    show.tip.label = TRUE,
    cex = 0.8,
    width = 8,
    height = 10
  )

  return(phylo)
}

## Synonym, vancouverensis is more accurate, but a new name so not in
## the tree
sp.network.metrics$GenusSpecies[sp.network.metrics$GenusSpecies ==
                                "Bombus vancouverensis"] <- "Bombus bifarius"

bombus.phylo <- makeSaveTree(sp.network.metrics, phylo, "Bombus")

melissodes.phylo <- makeSaveTree(sp.network.metrics, phylo, "Melissodes")
