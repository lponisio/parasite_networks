## Generate phylogeny of Bombus
rm(list=ls())
source("lab_paths.R")
setwd(local.path)
setwd("cascades-meadows/analysis/Phylogeny")


## load tree from :
## Henriquez Piskulich, Patricia Andrea; Hugall, Andrew F.; Stuart-Fox, Devi (2023). 
## A supermatrix phylogeny of the world’s bees (Hymenoptera: Anthophila) [Dataset]. 
## Dryad. https://doi.org/10.5061/dryad.80gb5mkw1
load('../../data/spec_net.Rdata')

phylo <- ape::read.tree("../../../cascades-meadows_saved/data/BEE_mat7_fulltree.nwk")

##clean up unwanted portion of labels
pattern <- "(_n\\d+m\\d+_[A-Za-z0-9]+)?$"
phylo$tip.label <- gsub(pattern, "", phylo$tip.label)

## replace underscore with space
phylo$tip.label <- gsub("_", " ", phylo$tip.label)

## i think we need to drop the tips that aren't in our dataset
## changing sp labels to match phylogeny
species_to_keep <- na.omit(unique(spec.net$GenusSpecies))
print(length(species_to_keep))


#the 6 species getting dropped are still present at this step

phylo_tips <- phylo$tip.label
print(length(phylo_tips))

#only keep tips that match our species
phylo <- ape::keep.tip(phylo, species_to_keep[species_to_keep %in% phylo_tips])
print(length(phylo$tip.label))

## Calculate covariance matrix
phylo_matrix <- ape::vcv.phylo(phylo)
save(phylo, phylo_matrix, file = "../../data/community_phylogeny.Rdata")
###############################################################################

## Do a bombus only phylogeny
species_to_keep <- na.omit(unique(spec.net$GenusSpecies[spec.net$Genus == "Bombus"]))
print(length(species_to_keep))


#the 6 species getting dropped are still present at this step

phylo_tips <- phylo$tip.label
print(length(phylo_tips))

#only keep tips that match our species
phylo <- ape::keep.tip(phylo, species_to_keep[species_to_keep %in% phylo_tips])
print(length(phylo$tip.label))

## Calculate covariance matrix
phylo_matrix <- ape::vcv.phylo(phylo)

save(phylo, phylo_matrix, file = "../../data/bombus_phylogeny.Rdata")
