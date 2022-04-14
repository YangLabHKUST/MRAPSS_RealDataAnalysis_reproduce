# devtools::install_version("mixsqp", version = "0.1-97", repos = "http://cran.us.r-project.org")
# devtools::install_version("ashr", version = "2.2-32", repos = "http://cran.us.r-project.org")
# devtools::install_github("jean997/cause@v1.0.0")
# for the latest version, use devtools::install_github("jean997/cause@v1.2.0")

library(readr)
library(cause)
ts1 = c("AD", "ASD", "Daytime_Sleepiness", "Height_UKB",  "Intelligence", "RA",      
        "T2D", "Alcohol", "BMI", "Depression", "IBD", "MDD", "SCZ", "Angina", 
        "CAD", "HBP", "Income", "NEB", "Smoking", "Urate", "Anorexia", 
        "CD", "Height_GIANT", "Insomnia", "Neuroticism", "SWB")

ts2 = c("Hair_Light_Brown", "Hair_Dark_Brown",  "Hair_Black",    "Hair_Blonde",  "Tanning")

Threshold=1e-03

for( exposure in ts1 ){
  
  for( outcome in ts2 ){
    
    start = proc.time()
    
    # read GWAS summary statistics
    cat(exposure,"~",outcome,"\n")
    
    X1 = readr::read_delim(paste0("./GWAS_26and5_formatted/", exposure), delim="\t", escape_double = FALSE, trim_ws = TRUE, progress = F)
    
    X2 = readr::read_delim(paste0("./GWAS_26and5_formatted/", outcome), delim="\t", escape_double = FALSE, trim_ws = TRUE, progress = F)
    
    X1$b = X1$Z/sqrt(X1$N)
    X2$b = X2$Z/sqrt(X2$N)
    X1$se = 1/sqrt(X1$N)
    X2$se = 1/sqrt(X2$N)
    
    X <- try(gwas_merge(X1, X2, 
                        snp_name_cols = c("SNP", "SNP"),
                        beta_hat_cols = c("b", "b"),
                        se_cols = c("se", "se"),
                        A1_cols = c("A1", "A1"),
                        A2_cols = c("A2", "A2")))
    
    if(inherits(X , 'try-error')) next
    
    d0 = X1[, c("SNP", "P")]
    colnames(d0) = c("snp", "pval.exp")
    X0 = merge(X, d0, by="snp")
    
    # clump
    clumped = MRAPSS::clump(X0,
                            IV.Threshold = 1e-03,
                            SNP_col = "snp",
                            pval_col = "pval.exp",
                            clump_kb = 1000,
                            clump_r2 = 0.1,
                            bfile = "/import/home/share/xhu/database/1KG/all_1000G_EUR_Phase3",
                            plink_bin = "/import/home/maxhu/plink/plink")
    
    varlist <- with(X, sample(snp, size=min(nrow(X), 1000000), replace=FALSE))
    params <- try(est_cause_params(X, varlist))
    if(inherits(params , 'try-error')) next
    
    # cause
    if(!is.null(clumped)){
      
      top_ldl_pruned_vars =intersect(as.character(X$snp), as.character(subset(clumped, pval.exp <= Threshold)$snp))
      
      cause_res <- try(cause(X=X, variants = top_ldl_pruned_vars , param_ests = params, force=TRUE))
      
      if(inherits( cause_res , 'try-error')) next
      
      res_elpd <- data.frame(exposure,
                             outcome,
                             Threshold,
                             length(top_ldl_pruned_vars),
                             cause_res$elpd)
      
      res.cause.est = summary(cause_res, ci_size=0.95)
      
      res = data.frame(exposure, outcome,
                       Threshold,length(top_ldl_pruned_vars),
                       matrix(c(res.cause.est$quants[[2]][,1],
                                res.cause.est$quants[[2]][,2],
                                res.cause.est$quants[[2]][,3]), nrow=1))
      
      write.table(res, file="NC_CAUSE_est", append=T,
                  col.names = F, row.names = F,
                  quote = F)
      
      write.table(res_elpd, file="NC_CAUSE_elpd", append=T,
                  col.names = F, row.names = F,
                  quote = F)
      
      rm(top_ldl_pruned_vars)
      rm(res)
      rm(res_elpd)
      rm(res.cause.est)
      rm(cause_res)
      print(proc.time()-start)
      
    }
    
  }
  
}

cause_elpd = unique(read.table("NC_CAUSE_elpd", header = F))
cause_est = unique(read.table("NC_CAUSE_est", header = F))
colnames(cause_elpd) = c("exposure","outcome","Threshold","nsnp","model1","model2","delta_elpd", "se_delta_elpd", "Z")
colnames(cause_est) = c("exposure","outcome","Threshold","nsnp", "b","b_l","b_u","eta","eta_l","eta_u","q","q_l","q_u")
cause_elpd = unique(subset(cause_elpd, model1=="sharing"&model2=="causal"))
cause_elpd$pval = pnorm(cause_elpd$Z)
cause_est = unique(cause_est[, c("exposure","outcome","Threshold","nsnp", "b","b_l","b_u")])
cause_est$se = (cause_est$b_u - cause_est$b_l)/2/1.96
cause_res = unique(merge(unique(cause_elpd[, c("exposure","outcome","pval")]),
                         cause_est[, c("exposure","outcome","Threshold","nsnp", "beta.hat","se")],
                         by=c("exposure","outcome")))
cause_res$Method = "CAUSE"
cause_res = cause_res[, c("exposure","outcome","Method", "Threshold", "nsnp", "beta.hat","se", "pval")]
write.table(cause_res, file="NC_CAUSE.MRres", append=F, col.names = T, row.names = F, quote = F)


