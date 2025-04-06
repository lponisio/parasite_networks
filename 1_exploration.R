library(ggplot2)
library(ggpubr)
library(gplots)
library(bipartite)
library(RColorBrewer)
library(tidyverse)

rm(list=ls())
setwd("~/University of Oregon Dropbox/Lauren Ponisio/")
setwd("parasite_networks")

source("src/misc.R")

load(file="../parasite_networks/data/sp_genus_site_mets.RData")

## ********************************************************
##  Heat maps
## ********************************************************

load(file="../parasite_networks/data/sp_genus_site_mets.RData")

## focus on only genera consistently surveyed across projects
genera <- c("Bombus", "Apis", "Melissodes", "Andrena", "Svastra")

par.mets <- par.mets[par.mets$Genus %in% genera,]

parasites <- c("AscosphaeraSpp", "ApicystisSpp",
               "CrithidiaExpoeki",
               "CrithidiaBombi", "CrithidiaSpp",
               "NosemaCeranae", "NosemaBombi")

parasite.cols <- c( "SpCrithidiaPresence",
                   paste0("Sp", parasites))

par.mets <- par.mets[par.mets$SpScreened > 5,]

par.mets <- par.mets[, c("GenusSpecies", "Genus",
                                   "ProjectSubProject",
                                   "SampleRound", "Year",
                                   "SpScreened",
                                   parasite.cols)]

sum.screened <- par.mets  %>%
  group_by(ProjectSubProject, GenusSpecies, Genus) %>%
  summarise(
    SpCrithidiaPresence= sum(SpCrithidiaPresence),
    SpApicystisSpp = sum(SpApicystisSpp),
    SpAscosphaeraSpp= sum(SpAscosphaeraSpp),
    SpCrithidiaExpoeki = sum(SpCrithidiaExpoeki),
    SpCrithidiaBombi = sum(SpCrithidiaBombi),
    SpCrithidiaSpp = sum(SpCrithidiaSpp),
    SpNosemaCeranae = sum(SpNosemaCeranae),
    SpNosemaBombi = sum(SpNosemaBombi),
    SpScreened = sum(SpScreened)
  )

sum.screened[, parasite.cols ] <- sum.screened[, parasite.cols
                                               ]/sum.screened$SpScreened

sp.sum <- par.mets  %>%
    group_by(GenusSpecies) %>%
    summarise(
        TotalScreened = sum(SpScreened)
    )

## add N to species/genus
sum.screened$GenusSpecies <- paste0(sum.screened$GenusSpecies, " (",
                                    sp.sum$TotalScreened[match(
                                        sum.screened$GenusSpecies,
                                        sp.sum$GenusSpecies)],
                                    ")")


## sum by site
site.sum <- par.mets  %>%
    group_by(ProjectSubProject) %>%
    summarise(
        TotalScreened = sum(SpScreened)
    )


## add N to sites
sum.screened$ProjectSubProject <- paste0(sum.screened$ProjectSubProject, " (",
                               site.sum$TotalScreened[match(
                                   sum.screened$ProjectSubProject,
                                   site.sum$ProjectSubProject)],
                               ")")

sum.screened$SpScreened <- NULL

## make community matrices 
makeParMat <- function(parasite, screened, sp.col="GenusSpecies"){
    colnames(screened)[colnames(screened) == parasite] <-
        "parasite"
    colnames(screened)[colnames(screened) == sp.col] <- "SpCat"
    par.mat <- screened %>%
        select(SpCat, ProjectSubProject, parasite) %>%
        pivot_wider(names_from = SpCat, values_from =
                                            parasite)
    par.mat <- as.data.frame(par.mat)
    rownames(par.mat) <- par.mat$ProjectSubProject
    par.mat$ProjectSubProject <- NULL
    par.mat <- as.matrix(par.mat)
    return(par.mat)
}

par.mats.species <- lapply(parasite.cols, makeParMat, sum.screened)
names(par.mats.species) <- parasite.cols

## heat maps of # of infected individuals by species
plotParasiteMap <- function(){
    colfunc <- colorRampPalette(c("grey", "red"))
    if(parasite == "SpAscosphaeraSpp"){
        par(oma=c(15,4,3,10), mar=c(1,2,2,1),
            mgp=c(1.5,0.5,0))
    } else{
        par(oma=c(10,4,3,8), mar=c(1,2,2,1),
            mgp=c(1.5,0.5,0))
    }
    heatmap.2(bipartite::empty(par.mats.species[[parasite]]),
              trace="none",
              col=colfunc,
              breaks=seq(0, 1, 0.1))
    mtext(parasite, 3, line=2.5)
}
for(parasite in parasite.cols){
    pdf.f(plotParasiteMap,
          file=sprintf("figures/heatmaps/%s.pdf", parasite),
          width=10, height=7)
}



## ********************************************************
##  Network metrics vs. parasite prevalence
## ********************************************************

load(file="../parasite_networks/data/network_mets.RData")
load(file="../parasite_networks/data/sp_mets.RData")


plotNetsCrithidia <- function(genus, cor.dats){
  this.genus <- cor.dats[cor.dats$Genus == genus,]
  ## 1
  nodf <- ggplot(this.genus, aes(x = zweighted.NODF,
                                 y = PropGenusCrithidiaPresence,
                                 color = ProjectSubProject)) +
    geom_point(size=2) +
    labs(y="Crithidia spp. prevalence") 
  h2 <- ggplot(this.genus, aes(x = zH2 ,
                               y = PropGenusCrithidiaPresence,
                               color = ProjectSubProject)) +
    geom_point(size=2) +
    labs(y="Crithidia spp. prevalence") 

  ## 2
  niche.overlap.hl <- ggplot(this.genus,
                             aes(x = zniche.overlap.HL,
                                 y = PropGenusCrithidiaPresence,
                                 color = ProjectSubProject)) +
    geom_point(size=2) +
    labs(y="Crithidia spp. prevalence") 
  niche.overlap.ll <- ggplot(this.genus,
                             aes(x = zniche.overlap.LL,
                                 y = PropGenusCrithidiaPresence,
                                 color = ProjectSubProject)) +
    geom_point(size=2) +
    labs(y="Crithidia spp. prevalence") 

  ## 3
  complementarity.hl <- ggplot(this.genus,
                               aes(x = zfunctional.complementarity.HL,
                                   y = PropGenusCrithidiaPresence,
                                   color = ProjectSubProject)) +
    geom_point(size=2) + 
    labs(y="Crithidia spp. prevalence") 
  complementarity.ll <- ggplot(this.genus,
                               aes(x = zfunctional.complementarity.LL,
                                   y = PropGenusCrithidiaPresence,
                                   color = ProjectSubProject)) +
    geom_point(size=2) +
    labs(y="Crithidia spp. prevalence") 

  ## 4
  site.div <- ggplot(this.genus,
                     aes(x = SiteBeeDiversity,
                         y = PropGenusCrithidiaPresence,
                         color = ProjectSubProject)) +
    geom_point(size=2) +
    labs(y="Crithidia spp. prevalence") 
  site.abund <- ggplot(this.genus,
                       aes(x = SiteRelativeBeeAbundance,
                           y = PropGenusCrithidiaPresence,
                           color = ProjectSubProject)) +
    geom_point(size=2) +
    labs(y="Crithidia spp. prevalence")

  ## 5
 
  genus.abund <- ggplot(this.genus,
                        aes(x = GenusRelativeAbundance,
                            y = PropGenusCrithidiaPresence,
                            color = ProjectSubProject)) +
    geom_point(size=2) +
    labs(y="Crithidia spp. prevalence") 

  all.net.plots <- ggarrange(nodf, h2,
                             niche.overlap.hl, niche.overlap.ll,
                             complementarity.hl, complementarity.ll,
                             site.div, site.abund, genus.abund,
                             nrow=5,ncol=2,
                             labels =
                               c("A", "B", "C", "D", "E", "F", "G"),
                             common.legend = TRUE, legend="bottom")

  ggsave(all.net.plots,
         file=sprintf("figures/%s_Crithidia.pdf",genus),
         height=13, width=8)

}

plotNetsCrithidia("Bombus", network.metrics)
plotNetsCrithidia("Apis", network.metrics)
plotNetsCrithidia("Melissodes", network.metrics)




plotNetsApicystisSpp <- function(genus, cor.dats){
  this.genus <- cor.dats[cor.dats$Genus == genus,]

  ## 1
  nodf <- ggplot(this.genus, aes(x = zweighted.NODF,
                                 y = PropGenusApicystisSpp,
                                 color = ProjectSubProject)) +
    geom_point(size=2) +
    labs(y="Apicystis spp. prevalence") 
  h2 <- ggplot(this.genus, aes(x = zH2 ,
                               y = PropGenusApicystisSpp,
                               color = ProjectSubProject)) +
    geom_point(size=2) +
    labs(y="Apicystis spp. prevalence") 

  ## 2
  niche.overlap.hl <- ggplot(this.genus,
                             aes(x = zniche.overlap.HL,
                                 y = PropGenusApicystisSpp,
                                 color = ProjectSubProject)) +
    geom_point(size=2) +
    labs(y="Apicystis spp. prevalence") 
  niche.overlap.ll <- ggplot(this.genus,
                             aes(x = zniche.overlap.LL,
                                 y = PropGenusApicystisSpp,
                                 color = ProjectSubProject)) +
    geom_point(size=2) +
    labs(y="Apicystis spp. prevalence")
  
  ## 3
  complementarity.hl <- ggplot(this.genus,
                               aes(x = zfunctional.complementarity.HL,
                                   y = PropGenusApicystisSpp,
                                   color = ProjectSubProject)) +
    geom_point(size=2) + 
    labs(y="Apicystis spp. prevalence") 
  complementarity.ll <- ggplot(this.genus,
                               aes(x = zfunctional.complementarity.LL,
                                   y = PropGenusApicystisSpp,
                                   color = ProjectSubProject)) +
    geom_point(size=2) +
    labs(y="Apicystis spp. prevalence") 

  ## 4
  site.div <- ggplot(this.genus,
                     aes(x = SiteBeeDiversity,
                         y = PropGenusApicystisSpp,
                         color = ProjectSubProject)) +
    geom_point(size=2) +
    labs(y="Apicystis spp. prevalence") 
  site.abund <- ggplot(this.genus,
                       aes(x = SiteRelativeBeeAbundance,
                           y = PropGenusApicystisSpp,
                           color = ProjectSubProject)) +
    geom_point(size=2) +
    labs(y="Apicystis spp. prevalence")

  ## 5
  genus.abund <- ggplot(this.genus,
                        aes(x = GenusRelativeAbundance,
                            y = PropGenusApicystisSpp,
                            color = ProjectSubProject)) +
    geom_point(size=2) +
    labs(y="Apicystis spp. prevalence") 

  all.net.plots <- ggarrange(nodf, h2,
                             niche.overlap.hl, niche.overlap.ll,
                             complementarity.hl, complementarity.ll,
                             site.div, site.abund, genus.abund,
                             nrow=5,ncol=2,
                             labels =
                               c("A", "B", "C", "D", "E", "F", "G"),
                             common.legend = TRUE, legend="bottom")


  ggsave(all.net.plots, file=sprintf("figures/%s_ApicystisSpp.pdf",genus),
         height=11, width=8)

}

plotNetsApicystisSpp("Bombus", network.metrics)
plotNetsApicystisSpp("Apis", network.metrics)
plotNetsApicystisSpp("Melissodes", network.metrics)



plotSpNetsCrithidia <- function(genus, cor.dats){
  this.genus <- cor.dats[cor.dats$Genus == genus,]

  ## 1
  between <- ggplot(this.genus, aes(x = zweighted.betweenness,
                                    y = PropSpCrithidiaPresence,
                                    color = ProjectSubProject)) +
    geom_point(size=2) +
    labs(y="Crithidia spp. prevalence") 
  closeness <- ggplot(this.genus, aes(x = zweighted.closeness,
                                      y = PropSpCrithidiaPresence,
                                      color = ProjectSubProject)) +
    geom_point(size=2) +
    labs(y="Crithidia spp. prevalence") 

  ## 2
  nestedrank <- ggplot(this.genus,
                       aes(x = znestedrank,
                           y = PropSpCrithidiaPresence,
                           color = ProjectSubProject)) +
    geom_point(size=2) +
    labs(y="Crithidia spp. prevalence") 
  d <- ggplot(this.genus,
              aes(x = zd,
                  y = PropSpCrithidiaPresence,
                  color = ProjectSubProject)) +
    geom_point(size=2) +
    labs(y="Crithidia spp. prevalence") 

  ## 3
  zdegree <- ggplot(this.genus,
                    aes(x = proportional.generality,
                        y = PropSpCrithidiaPresence,
                        color = ProjectSubProject)) +
    geom_point(size=2) + 
    labs(y="Crithidia spp. prevalence") 
  degree <- ggplot(this.genus,
                   aes(x = degree,
                       y = PropSpCrithidiaPresence,
                       color = ProjectSubProject)) +
    geom_point(size=2) +
    labs(y="Crithidia spp. prevalence")

    ## 4
  site.div <- ggplot(this.genus,
                     aes(x = SiteBeeDiversity,
                         y = PropSpCrithidiaPresence,
                         color = ProjectSubProject)) +
    geom_point(size=2) +
    labs(y="Crithidia spp. prevalence") 
  site.abund <- ggplot(this.genus,
                       aes(x = SiteRelativeBeeAbundance,
                           y = PropSpCrithidiaPresence,
                           color = ProjectSubProject)) +
    geom_point(size=2) +
    labs(y="Crithidia spp. prevalence")

  ## 5
  genus.abund <- ggplot(this.genus,
                        aes(x = GenusRelativeAbundance,
                            y = PropSpCrithidiaPresence,
                            color = ProjectSubProject)) +
    geom_point(size=2) +
    labs(y="Crithidia spp. prevalence") 
    

  all.net.plots <- ggarrange(between, closeness,
                             nestedrank, d,
                             zdegree, degree,
                             site.div, site.abund,
                             genus.abund,
                             nrow=5,ncol=2,
                             labels = c("A", "B", "C", "D",
                                        "E", "F", "G"),
                             common.legend = TRUE, legend="bottom")

  ggsave(all.net.plots, file=sprintf("figures/Sp_%s_Crithidia.pdf",genus),
         height=11, width=8)

}

plotSpNetsCrithidia("Bombus", sp.network.metrics)
plotSpNetsCrithidia("Apis", sp.network.metrics)
plotSpNetsCrithidia("Melissodes", sp.network.metrics)



plotSpNetsApicystisSpp <- function(genus, cor.dats){
  this.genus <- cor.dats[cor.dats$Genus == genus,]

  ## 1

  between <- ggplot(this.genus, aes(x = zweighted.betweenness,
                                    y = PropSpApicystisSpp,
                                    color = ProjectSubProject)) +
    geom_point(size=2) +
    labs(y="Apicystis spp. prevalence") 

  
  closeness <- ggplot(this.genus, aes(x = zweighted.closeness,
                                      y = PropSpApicystisSpp,
                                      color = ProjectSubProject)) +
    geom_point(size=2) +
    labs(y="Apicystis spp. prevalence") 

  ## 2
  nestedrank <- ggplot(this.genus,
                       aes(x = znestedrank,
                           y = PropSpApicystisSpp,
                           color = ProjectSubProject)) +
    geom_point(size=2) +
    labs(y="Apicystis spp. prevalence") 


  d <- ggplot(this.genus,
              aes(x = zd,
                  y = PropSpApicystisSpp,
                  color = ProjectSubProject)) +
    geom_point(size=2) +
    labs(y="Apicystis spp. prevalence") 


  ## 3
  zdegree <- ggplot(this.genus,
                    aes(x = proportional.generality,
                        y = PropSpApicystisSpp,
                        color = ProjectSubProject)) +
    geom_point(size=2) + 
    labs(y="Apicystis spp. prevalence") 


  degree <- ggplot(this.genus,
                   aes(x = degree,
                       y = PropSpApicystisSpp,
                       color = ProjectSubProject)) +
    geom_point(size=2) +
    labs(y="Apicystis spp. prevalence")

  ## 4
  site.div <- ggplot(this.genus,
                     aes(x = SiteBeeDiversity,
                         y = PropSpApicystisSpp,
                         color = ProjectSubProject)) +
    geom_point(size=2) +
    labs(y="Apicystis spp. prevalence") 
  site.abund <- ggplot(this.genus,
                       aes(x = SiteRelativeBeeAbundance,
                           y = PropSpApicystisSpp,
                           color = ProjectSubProject)) +
    geom_point(size=2) +
    labs(y="Apicystis spp. prevalence")

  ## 5
  genus.abund <- ggplot(this.genus,
                        aes(x = GenusRelativeAbundance,
                            y = PropSpApicystisSpp,
                            color = ProjectSubProject)) +
    geom_point(size=2) +
    labs(y="Apicystis spp. prevalence") 
  

  all.net.plots <- ggarrange(between, closeness,
                             nestedrank, d,
                             zdegree, degree,
                             site.div, site.abund,
                             genus.abund,
                             nrow=5,ncol=2,
                             labels = c("A", "B", "C", "D",
                                        "E", "F", "G"),
                             common.legend = TRUE, legend="bottom")
  
  ggsave(all.net.plots,
         file=sprintf("figures/Sp_%s_ApicystisSpp.pdf",genus),
         height=11, width=8)

}

plotSpNetsApicystisSpp("Bombus", sp.network.metrics)
plotSpNetsApicystisSpp("Apis", sp.network.metrics)
plotSpNetsApicystisSpp("Melissodes", sp.network.metrics)
