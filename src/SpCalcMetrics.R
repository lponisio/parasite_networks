library(igraph, quietly = TRUE)
library(bipartite, quietly = TRUE)



## function to simulate 1 null, and calculate statistics on it
SpCalcNullStat <- function(dat.web,
                         null.fun,...) {
    sim.web <- null.fun(dat.web)
    colnames(sim.web) <- colnames(dat.web)
    rownames(sim.web) <- rownames(dat.web)
    out.met <- spLevelWrap(sim.web,...)
    return(out.met)
}

calcZvals <- function(sp.level, true.stats, null.stats){
  ## calculate zvalues
  true.sp <- true.stats[[sp.level]]
  sp <- rownames(true.sp)
  rownames(true.sp) <- NULL
  
  null <- lapply(null.stats, function(x){
    x[[sp.level]]
  })

  all.mets <- c(list(true.sp), null)

  array <- array(unlist(all.mets), dim = c(nrow(all.mets[[1]]),
                                           ncol(all.mets[[1]]),
                                           length(all.mets)))

  array.mean <- apply(array, c(1, 2), mean, na.rm=TRUE)
  array.sd <- apply(array, c(1, 2), sd, na.rm=TRUE)

  zvalues <- (true.sp-array.mean)/array.sd
  
  colnames(zvalues) <- paste0("z", colnames(true.sp))

  out.mets <- cbind(true.sp, zvalues)
  out.mets$GenusSpecies <- sp
  out.mets$SpeciesType <- sp.level

  return(out.mets)
}


honeybeeOverlap <- function(web){
  ## web with pollinators as columns, plants as rows
  if("Apis mellifera" %in% colnames(web)){
    int.dist <- vegdist(t(web), method = "chao")
    int.sim <- 1 - as.matrix(int.dist)
    overlap <- int.sim["Apis mellifera", ]
  } else{
    overlap <- rep(0, length(colnames(web)))
    names(overlap) <- colnames(web)
                 
  }
  return(overlap)
}


spLevelWrap <- function(sim.web, ...){
  sp.lev <- specieslevel(sim.web,...)
  overlap <- honeybeeOverlap(sim.web)
  sp.lev$'higher level'$HBOverlap <- overlap
  sp.lev$'lower level'$HBOverlap <- NA
  return(sp.lev)
}


##  function that computes summary statistics on simulated null matrices
##  (nulls simulated from web N times)
SpCalcNetworkMetrics <- function(dat.web, N,
                               index, ...) {

    ## check that matrix is proper format (no empty row/col and no NAs)
    ## drop empty rows and columns
    dat.web <- as.matrix(bipartite::empty(dat.web))
    ## check to make sure emptied matrix is large enough
    ## to calculate statistics on
    if(all(dim(dat.web) >= 2)) {
        ## calculate null metrics
        null.stat <- replicate(N,
                               SpCalcNullStat(dat.web,
                                            null.fun= vaznull.fast,
                                            index=index, level="both", ...),
                               simplify=FALSE)
        ## calculate metrics from data
        true.stat <- spLevelWrap(dat.web,
                                index=index, level="both", ...)

      ## hb.overlap <- honeybeeOverlap(dat.web)

      mets.hl <- calcZvals("higher level", true.stat, null.stat)      
      mets.ll <- calcZvals("lower level", true.stat, null.stat)
      
      out <- rbind(mets.hl, mets.ll)
      ## out$HBOverlap <- hb.overlap[match(out$GenusSpecies, names(hb.overlap))]
      
      return(out)


    }
    return(rep(NA, length(index)*nrow(dat.web) + length(index)*ncol(dat.web)))
}

SpPrepDat <- function(cor.stats, spec.dat,
                    cols.to.keep,
                    net.type,
                    species.level= "higher level"){
  dats <- do.call(rbind, cor.stats)
  out <- data.frame(dats)
  out$Site <- sapply(strsplit(rownames(out), "\\."),
                     function(x) x[1])
  out$Year <-  as.factor(sapply(strsplit(rownames(out), "\\."),
                                function(x) x[2]))

  if(net.type == "YearSR"){
    out$SampleRound <-  as.factor(sapply(strsplit(rownames(out), "\\."),
                                         function(x) x[3]))
  }

  site.dats <- unique(spec.dat[, c(cols.to.keep)])
  site.dats$Site  <- as.character(site.dats$Site)
  out$Year  <- as.character(out$Year)
  site.dats$Year  <- as.character(site.dats$Year)

  site.dats <- as.data.frame(site.dats)

  out <- out[out$SpeciesType == species.level,]
  out$Genus <- sapply(strsplit(out$GenusSpecies, "[ ]"), function(x) x[1])
  
  colnames(site.dats)[colnames(site.dats) %in% colnames(out)]
  print(dim(out))
  net.mets <- merge(out, site.dats, all.xy=TRUE)
  
  print(dim(net.mets))
  rownames(net.mets) <- NULL
  return(net.mets)
}

