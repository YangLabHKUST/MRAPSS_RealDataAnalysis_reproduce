run_RAPS_func <- function(dat, Threshold, exposure="Trait 1", outcome="Trait 2"){

  res = NULL

  for(Thresh in Threshold){

    indx = which(dat$pval.exp <= as.numeric(Thresh))
    nsnp = length(indx)
    MRres = setNames(data.frame(matrix(nrow=1,ncol=8)),
                     c("exposure", "outcome","Threshold",
                       "method", "nsnp", "b", "se","pval"))

    if(nsnp<=3) next
    cat("Threshold:",Thresh,"\t","nSNPs:", nsnp, "\n")

    #Raps
    res.raps <- try(TwoSampleMR::mr_raps(dat$b.exp[indx], dat$b.out[indx],
                                         dat$se.exp[indx], dat$se.out[indx],
                                         parameters = list(over.dispersion = TRUE,
                                                           loss.function = "huber")))

    if(inherits(res.raps,"try-error")){
      next
    }



    MRres$outcome = outcome
    MRres$exposure = exposure
    MRres$Threshold = Thresh
    MRres$method <- "RAPS"
    MRres$nsnp = length(indx)
    MRres$b <- res.raps$b
    MRres$se <- res.raps$se
    MRres$pval <- res.raps$pval
    res <- rbind(res,MRres)

  }
  return(res)
}


#########################################################################
library(readr)
library(dplyr)
library(ggplot2)

setwd("/home/share/xhhu/MR-APSS/reproduce")
pairs = read.table("/home/share/xhhu/MR-APSS/reproduce/tested_pairs",header = T)
colnames(pairs) = c("exposure", "outcome")

for(i in 1:nrow(pairs)){

  if(i > nrow(pairs)) break
  cat("pair", i, "\n")
  time=proc.time()
  trait1.name = as.character(pairs[i,"exposure"])
  trait2.name = as.character(pairs[i, "outcome"])

  dat= try(read.table(paste0("./MRdata/", trait1.name,"~",trait2.name), header = T))
  if(inherits(dat , 'try-error')) next

  if(nrow(dat) < 4) next
  MRres = run_RAPS_func(dat,
                              Threshold=c(5e-05,5e-06,5e-07,5e-08),
                              exposure=trait1.name,outcome=trait2.name)

  write.table(MRres, "./res/RAPS.MRres",
              quote=F, col.names = F, append = T,row.names = F)

}


