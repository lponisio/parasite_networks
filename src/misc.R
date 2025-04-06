
pdf.f <- function(f, file, ...) {
    cat(sprintf('Writing %s\n', file))
    pdf(file, ...)
    on.exit(dev.off())
    f()
}

## add transparency to named colors
add.alpha <- function(col, alpha=0.2){
    apply(sapply(col, col2rgb)/255, 2,
          function(x)
              rgb(x[1], x[2], x[3],
                  alpha=alpha))
}


standardize <- function(x)
(x-mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)

unstandardize <- function(x, orig){
    (x*sd(orig, na.rm=TRUE)) + mean(orig, na.rm=TRUE)
}




standardize.axis <- function(x, orig)
(x-mean(orig, na.rm=TRUE))/sd(orig, na.rm=TRUE)


plot.res <- function(mod, mod.name){
    ## function tp plot diagnostic figures for mcmc
    pdf(sprintf("figures/diagnostics/%s_Diag.pdf", mod.name),
        height=11, width=8.5)
    plot(mod,  N = 4, ask = FALSE)
    dev.off()
}

## function to clean up white-space in a column of data (replaces all
## instances of white-space with " " and empty cells with ""
fix.white.space <- function(d) {
  d <- as.character(d)
  remove.first <- function(s) substr(s, 2, nchar(s))
  d <- gsub("      ", " ", d, fixed=TRUE)
  d <- gsub("     ", " ", d, fixed=TRUE)
  d <- gsub("    ", " ", d, fixed=TRUE)
  d <- gsub("   ", " ", d, fixed=TRUE)
  d <- gsub("  ", " ", d, fixed=TRUE)

  tmp <- strsplit(as.character(d), " ")
  d <- sapply(tmp, function(x) paste(x, collapse=" "))

  first <- substr(d, 1, 1)
  d[first==" "] <- remove.first(d[first==" "])
  d
}


## makes a summary table of loo results
makeLooTable <- function(parasite, ## parasite name (length =1)
                         genus, ## bee genus (length=1)
                         abundance.order, ## what abundance var in
                         ## model? (variable length depending on the
                         ## number of var sets)
                         loo.results ## list of model loo results
                         ){
    ## loo.results must be a list with the results in the same order
    ## as abundance.order

    print(paste("abundance.order matches list length", length(abundance.order)==length(loo.results)))
    if(length(abundance.order) == 2){
        out <- data.frame(parasite=parasite,
                          genus=genus,
                          abundVar= abundance.order,
                          Estimate=round( c(loo.results[[1]]$estimates["looic", "Estimate"],
                                            loo.results[[2]]$estimates["looic", "Estimate"]), 2),
                          SE=round(c(loo.results[[1]]$estimates["looic", "SE"],
                                     loo.results[[2]]$estimates["looic", "SE"]),
                                   2))
    } else if(length(abundance.order) == 3){
        out <- data.frame(parasite=parasite,
                          genus=genus,
                          abundVar= abundance.order,
                          Estimate=round( c(loo.results[[1]]$estimates["looic", "Estimate"],
                                            loo.results[[2]]$estimates["looic", "Estimate"],
                                            loo.results[[3]]$estimates["looic", "Estimate"]), 2),
                          SE=round(c(loo.results[[1]]$estimates["looic", "SE"],
                                     loo.results[[2]]$estimates["looic", "SE"],
                                     loo.results[[3]]$estimates["looic", "SE"]), 2))
    }
    return(out)
}

makeGenusSubset <- function(spec.data, genus){
    out <- spec.data
    out$WeightsPar[out$Genus != genus] <- 0
    out$WeightsSp[out$Genus != genus] <- 0
    return(out)
}
