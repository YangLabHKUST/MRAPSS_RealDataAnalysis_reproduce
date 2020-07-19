# CAUSE
# website
library(readr)
library(dplyr)
library(cause)
refPanel.dir = "/home/xianghonghu/1000G_EUR_Phase3_plink/all_1000G_EUR_Phase3"
setwd("/home/share/xhhu/MR-APSS/reproduce")

pairs = read.table("/home/share/xhhu/MR-APSS/reproduce/tested_pairs",header = T)
colnames(pairs) = c("exposure", "outcome")
Threshold = 1e-03

for(i in  a:(a+49)){

  cat('Batch', a , "~", a+49, "pair", i, "\n")

  if(i > nrow(pairs)) break

  # read GWAS summary statistics
  trait1.name = as.character(pairs[i,"exposure"])
  trait2.name = as.character(pairs[i, "outcome"])
  trait1.dir = paste0("/home/share/xhhu/MR-APSS/reproduce/GWAS20_Formmated/", trait1.name)
  trait2.dir = paste0("/home/share/xhhu/MR-APSS/reproduce/GWAS20_Formmated/", trait2.name)

  cat(trait1.name,"~",trait2.name,"\n")

  X1 = readr::read_delim(trait1.dir,delim="\t",
                         escape_double = FALSE,
                         trim_ws = TRUE,
                         progress = F)


  X2 = readr::read_delim(trait2.dir, "\t",
                         escape_double = FALSE,
                         trim_ws = TRUE,
                         progress = F)

  X1$b = X1$Z/sqrt(X1$N)
  X2$b = X2$Z/sqrt(X2$N)
  X1$se = 1/sqrt(X1$N)
  X2$se = 1/sqrt(X2$N)

  X <- try(gwas_merge(X1, X2, snp_name_cols = c("SNP", "SNP"),
                      beta_hat_cols = c("b", "b"),
                      se_cols = c("se", "se"),
                      A1_cols = c("A1", "A1"),
                      A2_cols = c("A2", "A2")))

  if(inherits(X , 'try-error')) next

  d0 = X1[, c("SNP", "P")]
  colnames(d0) = c("snp", "pval.exp")

  X0 = merge(X, d0, by="snp")

  clumped_3 = MRAPSS::clump(subset(X0, pval.exp < 1e-03),
                            SNP_col = "snp",
                            pval_col = "pval.exp",
                            clump_kb = 1000,
                            clump_r2 = 0.1,
                            bfile = "/home/share/1000G_Phase3_all/all_1000G_EUR_Phase3",
                            plink_bin = "plink")


  if(inherits(clumped_3, 'try-error')) {
    clumped_3 = NULL
  }else{
    write.table(clumped_3, file=paste0("./CAUSE_MRdata/",trait1.name,"~",trait2.name, "_cause.MRdat"), col.names = T, row.names = F, quote = F)
  }


  varlist <- with(X, sample(snp, size=min(nrow(X), 1000000), replace=FALSE))
  params <- try(est_cause_params(X, varlist))
  save(params, file=paste0("./CAUSE_MRdata/", trait1.name,"~",trait2.name, ".RData"))

  if(inherits(params , 'try-error')) next

  if(!is.null(clumped_3)){

    if(Threshold == 1e-03) clumped = clumped_3

    top_ldl_pruned_vars =intersect(as.character(X$snp), as.character(subset(clumped, pval.exp < Threshold)$snp))

    cause_res <- try(cause(X=X, variants = top_ldl_pruned_vars , param_ests = params))

    if(inherits( cause_res , 'try-error')) next

    res_elpd <- data.frame(trait1.name,
                           trait2.name,
                           Threshold,
                           length(top_ldl_pruned_vars),
                           cause_res$elpd)

    res.cause.est = summary(cause_res, ci_size=0.95)

    res = data.frame(trait1.name, trait2.name,
                     Threshold,length(top_ldl_pruned_vars),
                     matrix(c(res.cause.est$quants[[2]][,1],
                              res.cause.est$quants[[2]][,2],
                              res.cause.est$quants[[2]][,3]), nrow=1))

    write.table(res, file="./res/CAUSE_est.MRres", append=T,
                col.names = F, row.names = F,
                quote = F)

    write.table(res_elpd, file="./res/CAUSE_elpd.MRres", append=T,
                col.names = F, row.names = F,
                quote = F)

    rm(top_ldl_pruned_vars)
    rm(res)
    rm(res_elpd)
    rm(res.cause.est)
    rm(cause_res)

  }

}


