library(ggplot2)
library(ggpubr)
rm(list=ls())
setwd("~/University of Oregon Dropbox/Lauren Ponisio/")
setwd("parasite_networks")


load(file="../parasite_networks/data/network_mets.RData")
load(file="../parasite_networks/data/sp_mets.RData")



plotNetsCrithidia <- function(genus, cor.dats){
  this.genus <- cor.dats[cor.dats$Genus == genus,]
  ## 1
  nodf <- ggplot(this.genus, aes(x = zweighted.NODF,
                                 y = PropGenusCrithidiaPresence,
                                 color = Project)) +
    geom_point(size=2) +
    labs(y="Crithidia spp. prevalence") 
  h2 <- ggplot(this.genus, aes(x = zH2 ,
                               y = PropGenusCrithidiaPresence,
                               color = Project)) +
    geom_point(size=2) +
    labs(y="Crithidia spp. prevalence") 

  ## 2
  niche.overlap.hl <- ggplot(this.genus,
                             aes(x = zniche.overlap.HL,
                                 y = PropGenusCrithidiaPresence,
                                 color = Project)) +
    geom_point(size=2) +
    labs(y="Crithidia spp. prevalence") 
  niche.overlap.ll <- ggplot(this.genus,
                             aes(x = zniche.overlap.LL,
                                 y = PropGenusCrithidiaPresence,
                                 color = Project)) +
    geom_point(size=2) +
    labs(y="Crithidia spp. prevalence") 

  ## 3
  complementarity.hl <- ggplot(this.genus,
                               aes(x = zfunctional.complementarity.HL,
                                   y = PropGenusCrithidiaPresence,
                                   color = Project)) +
    geom_point(size=2) + 
    labs(y="Crithidia spp. prevalence") 
  complementarity.ll <- ggplot(this.genus,
                               aes(x = zfunctional.complementarity.LL,
                                   y = PropGenusCrithidiaPresence,
                                   color = Project)) +
    geom_point(size=2) +
    labs(y="Crithidia spp. prevalence") 

  ## 4
  site.div <- ggplot(this.genus,
                     aes(x = SiteBeeDiversity,
                         y = PropGenusCrithidiaPresence,
                         color = Project)) +
    geom_point(size=2) +
    labs(y="Crithidia spp. prevalence") 
  site.abund <- ggplot(this.genus,
                       aes(x = SiteRelativeBeeAbundance,
                           y = PropGenusCrithidiaPresence,
                           color = Project)) +
    geom_point(size=2) +
    labs(y="Crithidia spp. prevalence")

  ## 5
 
  genus.abund <- ggplot(this.genus,
                        aes(x = GenusRelativeAbundance,
                            y = PropGenusCrithidiaPresence,
                            color = Project)) +
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
                                 color = Project)) +
    geom_point(size=2) +
    labs(y="Apicystis spp. prevalence") 
  h2 <- ggplot(this.genus, aes(x = zH2 ,
                               y = PropGenusApicystisSpp,
                               color = Project)) +
    geom_point(size=2) +
    labs(y="Apicystis spp. prevalence") 

  ## 2
  niche.overlap.hl <- ggplot(this.genus,
                             aes(x = zniche.overlap.HL,
                                 y = PropGenusApicystisSpp,
                                 color = Project)) +
    geom_point(size=2) +
    labs(y="Apicystis spp. prevalence") 
  niche.overlap.ll <- ggplot(this.genus,
                             aes(x = zniche.overlap.LL,
                                 y = PropGenusApicystisSpp,
                                 color = Project)) +
    geom_point(size=2) +
    labs(y="Apicystis spp. prevalence")
  
  ## 3
  complementarity.hl <- ggplot(this.genus,
                               aes(x = zfunctional.complementarity.HL,
                                   y = PropGenusApicystisSpp,
                                   color = Project)) +
    geom_point(size=2) + 
    labs(y="Apicystis spp. prevalence") 
  complementarity.ll <- ggplot(this.genus,
                               aes(x = zfunctional.complementarity.LL,
                                   y = PropGenusApicystisSpp,
                                   color = Project)) +
    geom_point(size=2) +
    labs(y="Apicystis spp. prevalence") 

  ## 4
  site.div <- ggplot(this.genus,
                     aes(x = SiteBeeDiversity,
                         y = PropGenusApicystisSpp,
                         color = Project)) +
    geom_point(size=2) +
    labs(y="Apicystis spp. prevalence") 
  site.abund <- ggplot(this.genus,
                       aes(x = SiteRelativeBeeAbundance,
                           y = PropGenusApicystisSpp,
                           color = Project)) +
    geom_point(size=2) +
    labs(y="Apicystis spp. prevalence")

  ## 5
  genus.abund <- ggplot(this.genus,
                        aes(x = GenusRelativeAbundance,
                            y = PropGenusApicystisSpp,
                            color = Project)) +
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
                                    color = Project)) +
    geom_point(size=2) +
    labs(y="Crithidia spp. prevalence") 
  closeness <- ggplot(this.genus, aes(x = zweighted.closeness,
                                      y = PropSpCrithidiaPresence,
                                      color = Project)) +
    geom_point(size=2) +
    labs(y="Crithidia spp. prevalence") 

  ## 2
  nestedrank <- ggplot(this.genus,
                       aes(x = znestedrank,
                           y = PropSpCrithidiaPresence,
                           color = Project)) +
    geom_point(size=2) +
    labs(y="Crithidia spp. prevalence") 
  d <- ggplot(this.genus,
              aes(x = zd,
                  y = PropSpCrithidiaPresence,
                  color = Project)) +
    geom_point(size=2) +
    labs(y="Crithidia spp. prevalence") 

  ## 3
  zdegree <- ggplot(this.genus,
                    aes(x = proportional.generality,
                        y = PropSpCrithidiaPresence,
                        color = Project)) +
    geom_point(size=2) + 
    labs(y="Crithidia spp. prevalence") 
  degree <- ggplot(this.genus,
                   aes(x = degree,
                       y = PropSpCrithidiaPresence,
                       color = Project)) +
    geom_point(size=2) +
    labs(y="Crithidia spp. prevalence")

    ## 4
  site.div <- ggplot(this.genus,
                     aes(x = SiteBeeDiversity,
                         y = PropSpCrithidiaPresence,
                         color = Project)) +
    geom_point(size=2) +
    labs(y="Crithidia spp. prevalence") 
  site.abund <- ggplot(this.genus,
                       aes(x = SiteRelativeBeeAbundance,
                           y = PropSpCrithidiaPresence,
                           color = Project)) +
    geom_point(size=2) +
    labs(y="Crithidia spp. prevalence")

  ## 5
  genus.abund <- ggplot(this.genus,
                        aes(x = GenusRelativeAbundance,
                            y = PropSpCrithidiaPresence,
                            color = Project)) +
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
                                    color = Project)) +
    geom_point(size=2) +
    labs(y="Apicystis spp. prevalence") 

  
  closeness <- ggplot(this.genus, aes(x = zweighted.closeness,
                                      y = PropSpApicystisSpp,
                                      color = Project)) +
    geom_point(size=2) +
    labs(y="Apicystis spp. prevalence") 

  ## 2
  nestedrank <- ggplot(this.genus,
                       aes(x = znestedrank,
                           y = PropSpApicystisSpp,
                           color = Project)) +
    geom_point(size=2) +
    labs(y="Apicystis spp. prevalence") 


  d <- ggplot(this.genus,
              aes(x = zd,
                  y = PropSpApicystisSpp,
                  color = Project)) +
    geom_point(size=2) +
    labs(y="Apicystis spp. prevalence") 


  ## 3
  zdegree <- ggplot(this.genus,
                    aes(x = proportional.generality,
                        y = PropSpApicystisSpp,
                        color = Project)) +
    geom_point(size=2) + 
    labs(y="Apicystis spp. prevalence") 


  degree <- ggplot(this.genus,
                   aes(x = degree,
                       y = PropSpApicystisSpp,
                       color = Project)) +
    geom_point(size=2) +
    labs(y="Apicystis spp. prevalence")

  ## 4
  site.div <- ggplot(this.genus,
                     aes(x = SiteBeeDiversity,
                         y = PropSpApicystisSpp,
                         color = Project)) +
    geom_point(size=2) +
    labs(y="Apicystis spp. prevalence") 
  site.abund <- ggplot(this.genus,
                       aes(x = SiteRelativeBeeAbundance,
                           y = PropSpApicystisSpp,
                           color = Project)) +
    geom_point(size=2) +
    labs(y="Apicystis spp. prevalence")

  ## 5
  genus.abund <- ggplot(this.genus,
                        aes(x = GenusRelativeAbundance,
                            y = PropSpApicystisSpp,
                            color = Project)) +
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
