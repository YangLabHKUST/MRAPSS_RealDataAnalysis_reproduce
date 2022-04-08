library(MRAPSS)
library(readr)

ts1 = c("AD", "ASD", "Daytime_Sleepiness", "Height_UKB",  "Intelligence", "RA",      
        "T2D", "Alcohol", "BMI", "Depression", "IBD", "MDD", "SCZ", "Angina", 
        "CAD", "HBP", "Income", "NEB", "Smoking", "Urate", "Anorexia", 
        "CD", "Height_GIANT", "Insomnia", "Neuroticism", "SWB")

ts2 = c("Hair_Light_Brown", "Hair_Dark_Brown",  "Hair_Black",  "Hair_Blonde",  "Tanning")


for( exposure in ts1){
  for( outcome in ts2){
    
    
    # Start the clock!
    start = proc.time()
    
    # read in formatted GWAS data
    trait1.dir = paste0("./GWAS_26and5_formatted/", exposure)
    trait2.dir = paste0("./GWAS_26and5_formatted/", outcome)
    
    
    cat(exposure,"~",outcome,"\n")
    
    trait1 = readr::read_delim(trait1.dir, delim = "\t", escape_double = FALSE, trim_ws = TRUE, progress = F)
    
    trait2 = readr::read_delim(trait2.dir, delim = "\t", escape_double = FALSE, trim_ws = TRUE, progress = F)
    
    # estimate parameters for the background model
    paras = est_paras(dat1 = trait1,
                      dat2 = trait2,
                      trait1.name = exposure,
                      trait2.name = outcome,
                      h2.fix.intercept = F,
                      LDSC = T,
                      ldscore.dir = "eur_w_ld_chr")
    
    write.table(data.frame(exposure,
                           outcome,
                           nrow(paras$dat),
                           rg = paras$ldsc_res$rg,
                           rg.se = paras$ldsc_res$rg.se,
                           paras$ldsc_res$I[1,1],
                           paras$ldsc_res$I[2,2],
                           rho = paras$ldsc_res$I[1,2],
                           paras$ldsc_res$I.se[1,1],
                           paras$ldsc_res$I.se[2,2],
                           rho = paras$ldsc_res$I.se[1,2],
                           paras$ldsc_res$cov[1,1],  paras$ldsc_res$cov[1,2],  paras$ldsc_res$cov[2,2],
                           paras$ldsc_res$cov.se[1,1],  paras$ldsc_res$cov.se[1,2],  paras$ldsc_res$cov.se[2,2]),
                file = "pairs_ldsc_res", col.names = F, row.names = F, quote=F,append = T)
    
    if(inherits(paras, 'try-error')) next
    
    write.table(matrix(as.vector(paras$Omega), nrow=1), paste0("./pairs_bg_paras/", exposure, "~", outcome, "_Omega"),
                row.names = F, col.names = F, quote = F)
      
    write.table(matrix(as.vector(paras$C), nrow=1), paste0("./pairs_bg_paras/", exposure, "~", outcome, "_C"),
                row.names = F, col.names = F, quote = F)
    
    # clumping      
    cat("Begin clumping ...\n ")
    
    clumped =  clump(paras$dat,
                     IV.Threshold = 5e-05,
                     SNP_col = "SNP",
                     pval_col = "pval.exp",
                     clump_kb = 1000,
                     clump_r2 = 0.001,
                     bfile = "all_1000G_EUR_Phase3",
                     plink_bin = "plink")
    
    # For convenience, users can perform LD clumping through API and do not need to provide a reference panel, i.e.,
    # clumped =  clump(paras$dat,
    #                  IV.Threshold = 5e-05,
    #                  SNP_col = "SNP",
    #                  pval_col = "pval.exp",
    #                  clump_kb = 1000,
    #                  clump_r2 = 0.001)
    

    if(inherits(clumped , 'try-error')) next
    
    
    write.table(clumped, file = paste0("./MRdat/", exposure,"~",outcome), col.names=T, row.names=F, quote=F)
    
    # Stop the clock
    proc.time() - start
    
  }
  
}

