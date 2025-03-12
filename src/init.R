
save.dir0 <- "saved/"
if(!dir.exists(save.dir0)) {
  dir.create(save.dir0, showWarnings = FALSE)
}

save.dir <- "saved/tables"
if(!dir.exists(save.dir)) {
  dir.create(save.dir, showWarnings = FALSE)
}

fig.dir <- "figures"
if(!dir.exists(fig.dir)) {
  dir.create(fig.dir, showWarnings = FALSE)
}


fig.dir2 <- "figures/diagnostics"
if(!dir.exists(fig.dir2)) {
  dir.create(fig.dir2, showWarnings = FALSE)
}



plot.res <- function(mod, mod.name){
    ## function tp plot diagnostic figures for mcmc
    pdf(sprintf("figures/diagnostics/%s_Diag.pdf", mod.name),
        height=11, width=8.5)
    plot(mod,  N = 4, ask = FALSE)
    dev.off()
}
