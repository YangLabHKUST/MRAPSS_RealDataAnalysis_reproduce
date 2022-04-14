# TwosampleMR
# install.packages("remotes")
# remotes::install_github("MRCIEU/TwoSampleMR")

# dIVW
# install.packages("devtools")
# devtools::install_github("tye27/mr.divw")

# MRcML
# devtools::install_github("xue-hr/MRcML")

# MRMix
# devtools::install_github("gqi/MRMix")

# 26 exposures
ts1 = c("AD", "ASD", "Daytime_Sleepiness", "Height_UKB",  "Intelligence", "RA",      
        "T2D", "Alcohol", "BMI", "Depression", "IBD", "MDD", "SCZ", "Angina", 
        "CAD", "HBP", "Income", "NEB", "Smoking", "Urate", "Anorexia", 
        "CD", "Height_GIANT", "Insomnia", "Neuroticism", "SWB")

# five negative control outcomes
ts2 = c("Hair_Light_Brown", "Hair_Dark_Brown",  "Hair_Black",    "Hair_Blonde",  "Tanning")

Threshold = 5e-08 # note that the default IV threshold for the eight methods is 5e-08


start0 = proc.time()

for( exposure in ts1 ){
  
  for( outcome in ts2 ){
    
    # reading in data
    cat("*********************************\n")
    cat("Pair: ", exposure,"~", outcome," ")
    dat= try(read.table(paste0("./MRdat/", exposure,"~", outcome), header = T))
    if(inherits(dat , 'try-error')) next
  
    indx = which(dat$pval.exp <= Threshold)
    nsnp = length(indx)
    if(nsnp < 4) next
    
    # IVW
    res.IVW <- try(TwoSampleMR::mr_ivw_fe(dat$b.exp[indx], dat$b.out[indx], dat$se.exp[indx], dat$se.out[indx]))
    cat("IVW: beta.hat=", res.IVW$b, " se=", res.IVW$se, "pval=", res.IVW$pval, "\n")
    
    if(inherits(res.IVW,"try-error")) next
    
    write.table(data.frame(exposure = exposure, outcome = outcome, method =  "IVW", Threshold = Threshold,
                           nsnp = nsnp, beta =  res.IVW$b, se = res.IVW$se, pval = res.IVW$pval),
                "NC_8methods.MRres", quote=F, col.names = F, append = T,row.names = F) 
    
  }
}
print(proc.time()-start0)

start0 = proc.time()

for( exposure in ts1 ){
  
  for( outcome in ts2 ){
    
   # reading in data
    cat("*********************************\n")
    cat("Pair: ", exposure,"~", outcome," ")
    dat= try(read.table(paste0("./MRdat/", exposure,"~", outcome), header = T))
    if(inherits(dat , 'try-error')) next
  
    indx = which(dat$pval.exp <= Threshold)
    nsnp = length(indx)
    if(nsnp < 4) next
    
    ## dIVW
    lambda.opt =  mr.divw::mr.eo(0, dat$b.exp[indx], dat$b.out[indx],
                                  dat$se.exp[indx], dat$se.out[indx],
                                  pval.selection=dat$pval.exp[indx])$lambda.opt
    res.dIVW = mr.divw::mr.divw(dat$b.exp[indx], dat$b.out[indx],
                                dat$se.exp[indx], dat$se.out[indx],
                                pval.selection=dat$pval.exp[indx], lambda=lambda.opt)
    cat("dIVW: beta.hat=", res.dIVW$beta.hat, " se=", res.dIVW$beta.se, "\n")
      
    
    if(inherits(res.dIVW,"try-error")) next
    
    write.table(data.frame(exposure = exposure, outcome = outcome, method =  "dIVW", Threshold = Threshold,
                           nsnp = nsnp, beta = res.dIVW$beta.hat, se = res.dIVW$beta.se, 
                           pval = pchisq(res.dIVW$beta.hat^2/res.dIVW$beta.se^2, 1, lower.tail = F)),
                "NC_8methods.MRres", quote=F, col.names = F, append = T,row.names = F) 
    
  }
}
print(proc.time()-start0)

start0 = proc.time()

for( exposure in ts1 ){
  
  for( outcome in ts2 ){
    
   # reading in data
    cat("*********************************\n")
    cat("Pair: ", exposure,"~", outcome," ")
    dat= try(read.table(paste0("./MRdat/", exposure,"~", outcome), header = T))
    if(inherits(dat , 'try-error')) next
  
    indx = which(dat$pval.exp <= Threshold)
    nsnp = length(indx)
    if(nsnp < 4) next
      
    ##  RAPS
    raps.df = data.frame(beta.exposure = dat$b.exp[indx], beta.outcome = dat$b.out[indx],
                         se.exposure = dat$se.exp[indx], se.outcome = dat$se.out[indx])
    res.raps <- try(suppressWarnings(mr.raps::mr.raps(raps.df, diagnostics=F)))
    cat("RAPS: beta.hat=", res.raps$beta.hat, " se=", res.raps$beta.se, "\n")
    
    if(inherits(res.raps,"try-error")) next
    
    write.table(data.frame(exposure = exposure, outcome = outcome, method =  "RAPS", Threshold = Threshold,
                           nsnp = nsnp, beta =  res.raps$beta.hat, se = res.raps$beta.se,
                           pval = pchisq(res.raps$beta.hat^2/res.raps$beta.se^2, 1, lower.tail = F)),
                "NC_8methods.MRres", quote=F, col.names = F, append = T,row.names = F) 
    

    
  }
}
print(proc.time()-start0)

start0 = proc.time()

for( exposure in ts1 ){
  
  for( outcome in ts2 ){
    
    # reading in data
    cat("*********************************\n")
    cat("Pair: ", exposure,"~", outcome," ")
    dat= try(read.table(paste0("./MRdat/", exposure,"~", outcome), header = T))
    if(inherits(dat , 'try-error')) next
  
    indx = which(dat$pval.exp <= Threshold)
    nsnp = length(indx)
    if(nsnp < 4) next
    
    ## Egger
    res.egger <- try(TwoSampleMR::mr_egger_regression(dat$b.exp[indx], dat$b.out[indx],
                                                      dat$se.exp[indx], dat$se.out[indx]))
    cat("Egger: beta.hat=", res.egger$b, " se=", res.egger$se, "pval=", res.egger$pval, "\n")
    
    
    if(inherits(res.egger,"error")) next
    
    write.table(data.frame(exposure = exposure, outcome = outcome, method =  "Egger", Threshold = Threshold,
                           nsnp = nsnp, beta =  res.egger$b, se = res.egger$se, pval = res.egger$pval),
                "NC_8methods.MRres", quote=F, col.names = F, append = T,row.names = F) 
    
  }
}
print(proc.time()-start0)

# we used the old version of MRMix(2020) in our paper.
# The results for some trait pairs may be slightly different between the old version and the current verison.
start0 = proc.time()

for( exposure in ts1 ){
  
  for( outcome in ts2 ){
    
    # reading in data
    cat("*********************************\n")
    cat("Pair: ", exposure,"~", outcome, " ")
    dat= try(read.table(paste0("./MRdat/", exposure,"~", outcome), header = T))
    if(inherits(dat , 'try-error')) next
  
    indx = which(dat$pval.exp <= Threshold)
    nsnp = length(indx)
    if(nsnp < 4) next
    
    ## MRmix
    res.MRMix = try(MRMix::MRMix(dat[indx,]$b.exp, dat[indx,]$b.out, dat[indx,]$se.exp, dat[indx,]$se.out))
    cat("MRMix: beta.hat=", res.MRMix$theta, " se=", res.MRMix$SE_theta, "pval=", res.MRMix$pvalue_theta, "\n")
    
    if(inherits(res.MRMix,"try-error")) next
    
    write.table(data.frame(exposure = exposure, outcome = outcome, method = "MRMix", Threshold = Threshold,
                           nsnp = nsnp, beta = res.MRMix$theta, se = res.MRMix$SE_theta, pval = res.MRMix$pvalue_theta),
                "NC_8methods.MRres", quote=F, col.names = F, append = T,row.names = F) 
    
  }
}
print(proc.time()-start0)

# Weighted-mode uses bootstrap to estimate the standard error. 
# The standard errors could be slightly different since we have not set the random seed.
start0 = proc.time()

for( exposure in ts1 ){
  
  for( outcome in ts2 ){
    
    # reading in data
    cat("*********************************\n")
    cat("Pair: ", exposure,"~", outcome, " ")
    dat= try(read.table(paste0("./MRdat/", exposure,"~", outcome), header = T))
    if(inherits(dat , 'try-error')) next
  
    indx = which(dat$pval.exp <= Threshold)
    nsnp = length(indx)
    if(nsnp < 4) next
    
    ## MRmix
    res.median <- try(TwoSampleMR::mr_weighted_median(dat$b.exp[indx], dat$b.out[indx], dat$se.exp[indx], dat$se.out[indx]))
    cat("Weighted-median: beta.hat=", res.median$b, " se=", res.median$se, "pval=", res.median$pval, "\n")
    
    if(inherits(res.median,"try-error")) next
    
    write.table(data.frame(exposure = exposure, outcome = outcome, method =  "Weighted-median", Threshold = Threshold,
                           nsnp = nsnp, beta = res.median$b, se = res.median$se, pval = res.median$pval),
                "NC_8methods.MRres", quote=F, col.names = F, append = T,row.names = F) 
    
  }
}
print(proc.time()-start0)

# Weighted-mode uses bootstrap to estimate the standard error. 
# The standard errors could be slightly different since we have not set the random seed.
start0 = proc.time()

for( exposure in ts1 ){
  
  for( outcome in ts2 ){
    
    # reading in data
    cat("*********************************\n")
    cat("Pair: ", exposure,"~", outcome, " ")
    dat= try(read.table(paste0("./MRdat/", exposure,"~", outcome), header = T))
    if(inherits(dat , 'try-error')) next
  
    indx = which(dat$pval.exp <= Threshold)
    nsnp = length(indx)
    if(nsnp < 4) next

    res.mode <- try(TwoSampleMR::mr_weighted_mode(dat$b.exp[indx], dat$b.out[indx], dat$se.exp[indx], dat$se.out[indx]))
    cat("Weighted-mode: beta.hat=", res.mode$b, " se=", res.mode$se, "pval=", res.mode$pval, "\n")
    
    if(inherits(res.mode,"try-error")) next
    
    write.table(data.frame(exposure = exposure, outcome = outcome, method =  "Weighted-mode", Threshold = Threshold,
                           nsnp = nsnp, beta = res.mode$b, se = res.mode$se, pval = res.mode$pval),
                "NC_8methods.MRres", quote=F, col.names = F, append = T,row.names = F) 

    
  }
}
print(proc.time()-start0)


start0 = proc.time()

for( exposure in ts1 ){
  
  for( outcome in ts2 ){
    
   # reading in data
    cat("*********************************\n")
    cat("Pair: ", exposure,"~", outcome, " ")
    dat= try(read.table(paste0("./MRdat/", exposure,"~", outcome), header = T))
    if(inherits(dat , 'try-error')) next
  
    indx = which(dat$pval.exp <= Threshold)
    nsnp = length(indx)
    if(nsnp < 4) next
    
    ## 3. cML-MA
    final_dat = dat[indx,]
    n= min(median(1/final_dat$se.exp^2), median(1/final_dat$se.out^2))
    cML_result = try(suppressWarnings(MRcML::mr_cML_DP(final_dat$b.exp,
                                                       final_dat$b.out,
                                                       final_dat$se.exp,
                                                       final_dat$se.out,
                                                       n = n,
                                                       random_start = 10,
                                                       random_start_pert = 10,
                                                       random_seed = 1,
                                                       num_pert = 200)))
    cat("cML-MA: beta.hat=", cML_result$MA_BIC_DP_theta, " se=", cML_result$MA_BIC_DP_se, "pval=", cML_result$MA_BIC_DP_p, "\n")
    
    
    write.table(data.frame(exposure = exposure, outcome = outcome, method =  "cML-MA", Threshold = Threshold,
                           nsnp = nsnp, beta =  cML_result$MA_BIC_DP_theta, se = cML_result$MA_BIC_DP_se, pval = cML_result$MA_BIC_DP_p),
                "NC_8methods.MRres", quote=F, col.names = F, append = T,row.names = F) 
    
    if(inherits(cML_result,"try-error")) next
    
  }
}
print(proc.time()-start0)

MRmethods_res = read.table("NC_8methods.MRres", header = F)
colnames(MRmethods_res) = c("exposure","outcome","Method", "Threshold", "nsnp", "beta.hat","se", "pval")
write.table(MRmethods_res, "NC_8methods.MRres", quote=F, col.names = T, append = F,row.names = F)

