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

# All the results will be saved to the file "Traits_MRAPSS.MRres"
start = proc.time()

ts1 = c("AD", "ASD", "Daytime_Sleepiness", "Height_UKB",  "Intelligence", "RA",      
        "T2D", "Alcohol", "BMI", "Depression", "IBD", "MDD", "SCZ", "Angina", 
        "CAD", "HBP", "Income", "NEB", "Smoking", "Urate", "Anorexia", 
        "CD", "Height_GIANT", "Insomnia", "Neuroticism", "SWB")

ts2 = ts1

IV.Threshold = 5e-05

for( exposure in ts1){
  
  for( outcome in ts2){
    
    if(exposure==outcome) next
      
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
    
     
     # saving resuts
     write.table(res, "Traits_MRAPSS.MRres", quote=F, col.names = F, append = T,row.names = F)
     
     
}
  
}

print(proc.time()-start)

apss_res = read.table("Traits_MRAPSS.MRres", header = F)
colnames(apss_res) = c("exposure","outcome", "nsnp","pi","nvalid","sigma.sq","tau.sq","b","se","pval",
                        "Method", "cor.Threshold", "IVstrength", "Threshold")
write.table(apss_res, "Traits_MRAPSS.MRres", quote=F, col.names = T, append = F,row.names = F)


