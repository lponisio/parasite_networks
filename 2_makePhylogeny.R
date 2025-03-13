## Generate phylogeny of Bombus
rm(list=ls())
source("lab_paths.R")
setwd(local.path)
setwd("parasite_networks")

## load tree from :
## Henriquez Piskulich, Patricia Andrea; Hugall, Andrew F.; Stuart-Fox, Devi (2023). 
## A supermatrix phylogeny of the world’s bees (Hymenoptera: Anthophila) [Dataset]. 
## Dryad. https://doi.org/10.5061/dryad.80gb5mkw1

load(file="data/sp_mets.RData")

phylo <- ape::read.tree("data/BEE_mat7_fulltree.nwk")

##clean up unwanted portion of labels
pattern <- "(_n\\d+m\\d+_[A-Za-z0-9]+)?$"
phylo$tip.label <- gsub(pattern, "", phylo$tip.label)

## replace underscore with space
phylo$tip.label <- gsub("_", " ", phylo$tip.label)

makeSaveTree <- function(spec.dats, phylo, genus){
  spec.dats <- spec.dats[
    spec.dats$Genus %in% genus,]

  species_to_keep <- na.omit(unique(spec.dats$GenusSpecies))
  print(length(species_to_keep))
  species_to_keep

  phylo_tips <- phylo$tip.label
  print(length(phylo_tips))

  ## only keep tips that match our species
  phylo <- ape::keep.tip(phylo, species_to_keep[species_to_keep %in%
                                                phylo_tips])
  print(length(phylo$tip.label))

  print(paste("missing", species_to_keep[!species_to_keep %in%
                                         phylo$tip.label]))
  ## Calculate covariance matrix
  phylo_matrix <- ape::vcv.phylo(phylo)
  save(phylo, phylo_matrix, file = sprintf("data/%s_phylogeny.Rdata",
                                           genus))
  return(phylo)
}


bombus.phylo <- makeSaveTree(sp.network.metrics, phylo, "Bombus")

melissodes.phylo <- makeSaveTree(sp.network.metrics, phylo, "Melissodes")
