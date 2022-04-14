run_APSS_func <- function(clumped=NULL,
                          exposure=NULL,
                          outcome=NULL,
                          C = diag(2),
                          Omega = matrix(0, 2, 2),
                          IV.Threshold = 5e-08, # IV selection threshold
                          Threshold = 5e-08,    # threshold for correcting for selection bias
                          Cor.SelectionBias=T){
  
 
  res=NULL
  
  for(i in 1:length(Threshold)){
    
    start = proc.time()
    if(!is.null(clumped)){
      
      test = subset(clumped, pval.exp <= IV.Threshold[i])
      if(nrow(test) < 4 ) next
      test$Threshold = Threshold[i]
      cat("IV selection threshold:", IV.Threshold[i] ,"\n")
      
      
      MRres = try(MRAPSS::MRAPSS(test,
                                 exposure=exposure,
                                 outcome=outcome,
                                 C= C,
                                 Omega= Omega,
                                 Cor.SelectionBias = Cor.SelectionBias))
 
      print(proc.time()-start)
      
      if(inherits(MRres, 'try-error')) {
        MRres=NULL
      }
    }
    
    
    res0 = data.frame(exposure = MRres$exposure,
                      outcome = MRres$outcome,
                      nSNP = nrow(MRres$MRdat),
                      pi0 = MRres$pi0,
                      nvalid = nrow(MRres$MRdat)*MRres$pi0,
                      sigma.sq= MRres$sigma.sq,
                      tau.sq= MRres$tau.sq,
                      beta = MRres$beta,
                      se = MRres$beta.se,
                      pval= MRres$pvalue,
                      method = MRres$method,
                      Threshold = Threshold[i],
                      IVstrength = MRres$IVsignal.sum
    )
    
    if(nrow(res0)!=0){
      res0$IV.Threshold = IV.Threshold[i]
    }
    
    res = rbind(res, res0)
    
  }
  
  return(res)
  
}



# All the results will be saved to the file "NC_MRAPSS.MRres"
start = proc.time()
ts1 = c("AD", "ASD", "Daytime_Sleepiness", "Height_UKB",  "Intelligence", "RA",      
        "T2D", "Alcohol", "BMI", "Depression", "IBD", "MDD", "SCZ", "Angina", 
        "CAD", "HBP", "Income", "NEB", "Smoking", "Urate", "Anorexia", 
        "CD", "Height_GIANT", "Insomnia", "Neuroticism", "SWB")

ts2 = c("Hair_Light_Brown", "Hair_Dark_Brown",  "Hair_Black",    "Hair_Blonde",  "Tanning")

IV.Threshold = 5e-05

for( exposure in ts1){
  
  for( outcome in ts2){
    
    cat("Pair: ", exposure,"~", outcome,"\n")
    # read in GWAS summary data for IVs
    clumped = try(read.table(paste0("./MRdat/", exposure,"~",outcome), header = T))
    
    # read in background parameters Omega and C
    C = try(as.matrix(read.table(paste0("./pairs_bg_paras/",exposure, "~", outcome,"_C"),
                                 header = F)))
    
    C = matrix(C[nrow(C), ], 2,2)
    
    Omega = try(as.matrix(read.table(paste0("./pairs_bg_paras/", exposure, "~", outcome,"_Omega"),
                                     header = F)))
      
    Omega = matrix(Omega[nrow(Omega), ], 2, 2)
    
    if(inherits(clumped , 'try-error')) clumped = NULL
    if(inherits(C, 'try-error')) next
    if(inherits(Omega, 'try-error')) next
    if(nrow(clumped) < 4 ) next
    
    # The p-value threshold for selection bias correction
    Threshold =  ifelse(IV.Threshold==5e-05, unique(clumped$Threshold), IV.Threshold)
    
    # MR-APSS
    cat("Run MR-APSS ... \n")
    res = run_APSS_func(clumped = clumped,
                        exposure = exposure,
                        outcome = outcome,
                        C = C,
                        Omega=Omega,
                        IV.Threshold = IV.Threshold,
                        Threshold = Threshold,
                        Cor.SelectionBias = T)
    
    # MR-APSS(Omega=0)
    cat("Run MR-APSS (Omega = 0), i.e. MR-APSS not accounting for polygenicity and pleiotropy. ... \n")
    res_Omega0 = run_APSS_func(clumped = clumped,
                               exposure = exposure,
                               outcome = outcome,
                               C = C,
                               Omega=matrix(0, 2, 2),
                               IV.Threshold = IV.Threshold,
                               Threshold = Threshold,
                               Cor.SelectionBias = T)
    res_Omega0$method = "MR-APSS(Omega=0)"
     
     # MR-APSS(C=I)
     cat("Run MR-APSS (C = I), i.e. MR-APSS not accounting for sample structure ... \n")
     res_CI = run_APSS_func(clumped = clumped,
                                exposure = exposure,
                                outcome = outcome,
                                C = diag(2),
                                Omega = Omega,
                                IV.Threshold = IV.Threshold,
                                Threshold = Threshold,
                                Cor.SelectionBias = T)
     res_CI$method = "MR-APSS(C=I)"
     
     # MR-APSS not accounting for selection bias
     cat("Run MR-APSS not accounting for selection bias ... \n")
     res_selec0 = run_APSS_func(clumped = clumped,
                            exposure = exposure,
                            outcome = outcome,
                            C = C,
                            Omega = Omega,
                            IV.Threshold = IV.Threshold,
                            Threshold = Threshold,
                            Cor.SelectionBias = F)
     res_selec0$method = "MR-APSS(Cor.SelectionBia=F)"
     
     # saving resuts
     write.table(res, "NC_MRAPSS.MRres", quote=F, col.names = F, append = T,row.names = F)
     
     write.table(res_Omega0, "NC_MRAPSS.MRres", quote=F, col.names = F, append = T,row.names = F)
     
     write.table(res_CI, "NC_MRAPSS.MRres", quote=F, col.names = F, append = T,row.names = F)
     
     write.table(res_selec0, "NC_MRAPSS.MRres", quote=F, col.names = F, append = T,row.names = F)
     
}
  
}

print(proc.time()-start)


apss_res = read.table("NC_MRAPSS.MRres", header = F)
colnames(apss_res) = c("exposure","outcome", "nsnp","pi","nvalid","sigma.sq","tau.sq","b","se","pval",
                        "Method", "cor.Threshold", "IVstrength", "Threshold")
head(apss_res)
write.table(apss_res, "NC_MRAPSS.MRres", quote=F, col.names = T, append = F,row.names = F)

# All the results will be saved to the file "NC_MRAPSS.MRres"
start = proc.time()
ts1 = c("AD", "ASD", "Daytime_Sleepiness", "Height_UKB",  "Intelligence", "RA",      
        "T2D", "Alcohol", "BMI", "Depression", "IBD", "MDD", "SCZ", "Angina", 
        "CAD", "HBP", "Income", "NEB", "Smoking", "Urate", "Anorexia", 
        "CD", "Height_GIANT", "Insomnia", "Neuroticism", "SWB")

ts2 = c("Hair_Light_Brown", "Hair_Dark_Brown",  "Hair_Black",    "Hair_Blonde",  "Tanning")

IV.Threshold = c(5e-06,5e-07, 5e-08)

for( exposure in ts1){
  
  for( outcome in ts2){
    
    cat("Pair: ", exposure,"~", outcome,"\n")
    # read in GWAS summary data for IVs
    clumped = try(read.table(paste0("./MRdat/", exposure,"~",outcome), header = T))
    
    # read in background parameters Omega and C
    C = try(as.matrix(read.table(paste0("./pairs_bg_paras/",exposure, "~", outcome,"_C"),
                                 header = F)))
    
    C = matrix(C[nrow(C), ], 2,2)
    
    Omega = try(as.matrix(read.table(paste0("./pairs_bg_paras/", exposure, "~", outcome,"_Omega"),
                                     header = F)))
      
    Omega = matrix(Omega[nrow(Omega), ], 2, 2)
    
    if(inherits(clumped , 'try-error')) clumped = NULL
    if(inherits(C, 'try-error')) next
    if(inherits(Omega, 'try-error')) next
    if(nrow(clumped) < 4 ) next
    
    # The p-value threshold for selection bias correction
    Threshold =  ifelse(IV.Threshold==5e-05, unique(clumped$Threshold), IV.Threshold)
    
    # MR-APSS
    cat("Run MR-APSS ... \n")
    res = run_APSS_func(clumped = clumped,
                        exposure = exposure,
                        outcome = outcome,
                        C = C,
                        Omega=Omega,
                        IV.Threshold = IV.Threshold,
                        Threshold = Threshold,
                        Cor.SelectionBias = T)
    
     
     # MR-APSS not accounting for selection bias
     cat("Run MR-APSS not accounting for selection bias ... \n")
     res_selec0 = run_APSS_func(clumped = clumped,
                            exposure = exposure,
                            outcome = outcome,
                            C = C,
                            Omega = Omega,
                            IV.Threshold = IV.Threshold,
                            Threshold = Threshold,
                            Cor.SelectionBias = F)
     res_selec0$method = "MR-APSS(Cor.SelectionBia=F)"
     
#      # saving resuts
     write.table(res, "NC_MRAPSS.MRres", quote=F, col.names = F, append = T,row.names = F)
     
     write.table(res_selec0, "NC_MRAPSS.MRres", quote=F, col.names = F, append = T,row.names = F)
     
}
  
}

print(proc.time()-start)
