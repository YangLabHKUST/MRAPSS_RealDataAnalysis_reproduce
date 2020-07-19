#######################################

run_APSS_func <- function(clumped=NULL,
                          clumped_rev=NULL,
                          trait1.name=NULL,
                          trait2.name=NULL,
                          Sigma_err = diag(2),
                          Omega = matrix(0, 2,2),
                          Threshold=5e-08,
                          tol = 1e-08,
                          cor.selecBias = T){

  time <- proc.time()

  cat("Begining MR tests ... \n")
  res=NULL

  for(i in 1:length(Threshold)){

    cat("Selection Threshold =",  Threshold[i], "\n")
    if(cor.selecBias == T) {
      Thresh=Threshold[i]
    }else{
      Thresh=1
    }

    MRres = NULL
    MRres_rev = NULL

    if(!is.null(clumped)){

      test = subset(clumped, pval.exp < Threshold[i])

      MRres = try(MRAPSS::MRAPSS(test,
                                 exposure=trait1.name,
                                 outcome=trait2.name,
                                 Sigma_err= Sigma_err,
                                 Omega= Omega,
                                 Threshold=Thresh))

      if(inherits(MRres, 'try-error')) {
        MRres=NULL
      }
    }

    if(!is.null(clumped_rev)){

      test_rev = subset(clumped_rev, pval.exp <= Threshold[i])

      MRres_rev = try(MRAPSS::MRAPSS(test_rev,
                                     exposure=trait2.name,
                                     outcome=trait1.name,
                                     Sigma_err= matrix(rev(Sigma_err),2,2),
                                     Omega = matrix(rev(Omega),2,2),
                                     Threshold=Thresh))


      if(inherits(MRres_rev, 'try-error')){
        MRres_rev=NULL
      }


    }

    res0 = data.frame(exposure=c( MRres$exposure,
                                  MRres_rev$exposure),
                      outcome=c(MRres$outcome,
                                MRres_rev$outcome),
                      nSNP = c(nrow(MRres$MRdat),
                               nrow(MRres_rev$MRdat)),
                      pi0 = c(MRres$pi0,
                              MRres_rev$pi0),
                      sigma.sq= c(MRres$sigma.sq,
                                  MRres_rev$sigma.sq),
                      tau.sq= c(MRres$tau.sq,
                                MRres_rev$tau.sq),
                      beta = c(MRres$beta,
                               MRres_rev$beta),
                      se = c(MRres$beta.se,
                             MRres_rev$beta.se),
                      pval=c(MRres$pvalue,
                             MRres_rev$pvalue),
                      method=c(MRres$method,
                               MRres_rev$method))

    if(nrow(res0)!=0){
      res0$Threshold= Threshold[i]
    }
    res = rbind(res, res0)
  }

  time_all <- proc.time()-time

  print(res)
  cat("Time elapsed", time_all[3],"\n")

  return(res)

}


library(readr)
library(dplyr)
library(MRAPSS)
library(mvtnorm)

setwd("/home/share/xhhu/MR-APSS/reproduce")
pairs = read.table("/home/share/xhhu/MR-APSS/reproduce/tested_pairs",header = T)
colnames(pairs) = c("exposure", "outcome")
threshs = c(5e-05, 5e-06, 5e-07, 5e-08)


for(i in  a:(a+49)){

  cat('Batch', a , "~", a+49, "pair", i, "\n")

  if(i > nrow(pairs)) break

  time=proc.time()

  trait1.name = as.character(pairs[i,"exposure"])
  trait2.name = as.character(pairs[i, "outcome"])

  cat(trait1.name,"~",trait2.name,"\n")

  clumped = try(read.table(paste0("./MRdata/", trait1.name,"~",trait2.name), header = T))

  Sigma_err = try(as.matrix(read.table(paste0("./Omega_est/",trait1.name, "~", trait2.name,"_Sigma_LD"),
                                       header = F)))

  Sigma_err = matrix(Sigma_err[nrow(Sigma_err), ], 2,2)

  Omega = try(as.matrix(read.table(paste0("./Omega_est/", trait1.name, "~", trait2.name,"_Omega"),
                                   header = F)))
  Omega = matrix(Omega[nrow(Omega), ], 2, 2)

  if(inherits(clumped , 'try-error')) clumped = NULL

  if(inherits(Sigma_err, 'try-error')) next
  if(inherits(Omega, 'try-error')) next

  # APSS
  res = run_APSS_func(clumped=clumped,
                      trait1.name=trait1.name,
                      trait2.name=trait2.name,
                      Sigma_err=Sigma_err,
                      Omega=Omega,
                      Threshold=threshs,
                      cor.selecBias = T)


  write.table(res, "./res/MRAPSS.MRres",
              quote=F, col.names = F, append = T,row.names = F)

  # Omega = 0

  res_omega0 = run_APSS_func(clumped = clumped,
                             trait1.name=trait1.name,
                             trait2.name=trait2.name,
                             Sigma_err=Sigma_err,
                             Threshold=threshs,
                             cor.selecBias = T)

  write.table(res_omega0, "./res/MRAPSS_omega0.MRres",
              quote=F, col.names = F, append = T,row.names = F)

  # c = 0
  res_c0 = run_APSS_func(clumped=clumped,
                         trait1.name=trait1.name,
                         trait2.name=trait2.name,
                         Sigma_err=Sigma_err,
                         Omega=Omega,
                         Threshold=threshs,
                         cor.selecBias = F)

  write.table(res_c0, "./res/MRAPSS_c0.MRres",
              quote=F, col.names = F, append = T,row.names = F)

  # rho = 0
  res_rho0 = run_APSS_func(clumped=clumped,
                           trait1.name=trait1.name,
                           trait2.name=trait2.name,
                           Omega=Omega,
                           Threshold=threshs,
                           cor.selecBias = T)


  write.table(res_rho0, "./res/MRAPSS_rho0.MRres",
              quote=F, col.names = F, append = T,row.names = F)

  print(proc.time()-time)
}


