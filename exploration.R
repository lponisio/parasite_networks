library(ggplot2)
library(ggpubr)
rm(list=ls())
setwd("~/University of Oregon Dropbox/Lauren Ponisio/")
setwd("parasite_networks")


load(file="../parasite_networks/data/network_mets.RData")


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
  cluster.hl <- ggplot(this.genus,
                       aes(x = zweighted.cluster.coefficient.HL,
                           y = PropGenusCrithidiaPresence,
                           color = Project)) +
    geom_point(size=2) +
    labs(y="Crithidia spp. prevalence") 


  cluster.ll <- ggplot(this.genus,
                       aes(x = zweighted.cluster.coefficient.LL,
                           y = PropGenusCrithidiaPresence,
                           color = Project)) +
    geom_point(size=2) +
    labs(y="Crithidia spp. prevalence") 

  all.net.plots <- ggarrange(nodf, h2,
                             niche.overlap.hl, niche.overlap.ll,
                             complementarity.hl, complementarity.ll,
                             cluster.hl, cluster.ll,
                             nrow=4,ncol=2,
                             labels = c("A", "B", "C", "D", "E", "F"),
                             common.legend = TRUE, legend="bottom")

  ggsave(all.net.plots, file=sprintf("figures/%s_Crithidia.pdf",genus),
         height=11, width=8)

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
    labs(y="ApicystisSpp spp. prevalence") 

  
  h2 <- ggplot(this.genus, aes(x = zH2 ,
                               y = PropGenusApicystisSpp,
                               color = Project)) +
    geom_point(size=2) +
    labs(y="ApicystisSpp spp. prevalence") 

  ## 2
  niche.overlap.hl <- ggplot(this.genus,
                             aes(x = zniche.overlap.HL,
                                 y = PropGenusApicystisSpp,
                                 color = Project)) +
    geom_point(size=2) +
    labs(y="ApicystisSpp spp. prevalence") 


  niche.overlap.ll <- ggplot(this.genus,
                             aes(x = zniche.overlap.LL,
                                 y = PropGenusApicystisSpp,
                                 color = Project)) +
    geom_point(size=2) +
    labs(y="ApicystisSpp spp. prevalence") 


  ## 3
  complementarity.hl <- ggplot(this.genus,
                               aes(x = zfunctional.complementarity.HL,
                                   y = PropGenusApicystisSpp,
                                   color = Project)) +
    geom_point(size=2) + 
    labs(y="ApicystisSpp spp. prevalence") 


  complementarity.ll <- ggplot(this.genus,
                               aes(x = zfunctional.complementarity.LL,
                                   y = PropGenusApicystisSpp,
                                   color = Project)) +
    geom_point(size=2) +
    labs(y="ApicystisSpp spp. prevalence") 



  ## 4
  cluster.hl <- ggplot(this.genus,
                       aes(x = zweighted.cluster.coefficient.HL,
                           y = PropGenusApicystisSpp,
                           color = Project)) +
    geom_point(size=2) +
    labs(y="ApicystisSpp spp. prevalence") 


  cluster.ll <- ggplot(this.genus,
                       aes(x = zweighted.cluster.coefficient.LL,
                           y = PropGenusApicystisSpp,
                           color = Project)) +
    geom_point(size=2) +
    labs(y="ApicystisSpp spp. prevalence") 

  all.net.plots <- ggarrange(nodf, h2,
                             niche.overlap.hl, niche.overlap.ll,
                             complementarity.hl, complementarity.ll,
                             cluster.hl, cluster.ll,
                             nrow=4,ncol=2,
                             labels = c("A", "B", "C", "D", "E", "F"),
                             common.legend = TRUE, legend="bottom")

  ggsave(all.net.plots, file=sprintf("figures/%s_ApicystisSpp.pdf",genus),
         height=11, width=8)

}

plotNetsApicystisSpp("Bombus", network.metrics)
plotNetsApicystisSpp("Apis", network.metrics)
plotNetsApicystisSpp("Melissodes", network.metrics)
