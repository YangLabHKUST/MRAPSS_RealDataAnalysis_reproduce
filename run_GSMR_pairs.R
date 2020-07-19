run_gsmr <- function(gsmr_data, exposure, outcome, Threshold, ldrho){

  bzx = gsmr_data$b.exp     # effects of the instruments on risk factor
  bzx_se = gsmr_data$se.exp             # standard errors of bzx
  bzx_pval = gsmr_data$pval.exp
  bzy = gsmr_data$b.out     # SNP effects on the disease
  bzy_se = gsmr_data$se.out             # standard errors of bzy
  bzy_pval = gsmr_data$pval.out                    # p-values for bzy

  n_ref = 489                        # Sample size of the reference sample
  gwas_thresh = Threshold*1.000001              # GWAS threshold to select SNPs as the instruments for the GSMR analysis
  single_snp_heidi_thresh = 0.01     # p-value threshold for single-SNP-based HEIDI-outlier analysis
  multi_snp_heidi_thresh = 0.01      # p-value threshold for multi-SNP-based HEIDI-outlier analysis
  nsnps_thresh = 0               # the minimum number of instruments required for the GSMR analysis
  heidi_outlier_flag = T             # flag for HEIDI-outlier analysis
  ld_r2_thresh = 1               # LD r2 threshold to remove SNPs in high LD
  ld_fdr_thresh = 0.05              # FDR threshold to remove the chance correlations between the SNP instruments
  gsmr2_beta = 0                     # 0 - the original HEIDI-outlier method; 1 - the new HEIDI-outlier method that is currently under development
  gsmr_results = gsmr::gsmr(bzx, bzx_se, bzx_pval, bzy, bzy_se, bzy_pval, ldrho, snp_coeff_id, n_ref, heidi_outlier_flag, gwas_thresh, single_snp_heidi_thresh, multi_snp_heidi_thresh, nsnps_thresh, ld_r2_thresh, ld_fdr_thresh, gsmr2_beta)    # GSMR analysis

  MRres = setNames(data.frame(matrix(nrow=1,ncol=8)),
                   c("exposure", "outcome","Threshold",
                     "method", "nsnp", "b", "se","pval"))

  MRres$outcome = outcome
  MRres$exposure = exposure
  MRres$Threshold = Threshold
  MRres$method = "GSMR"
  MRres$nsnp = length(gsmr_results$used_index)
  MRres$b <- c(gsmr_results$bxy)
  MRres$se <- c(gsmr_results$bxy_se)
  MRres$pval <- c(gsmr_results$bxy_pval)

  return(MRres)

}

# run gsmr

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

  clumped = try(read.table(paste0("./MRdata/", trait1.name,"~",trait2.name), header = T))
  if( inherits(clumped, 'try-error'))  next
  snp_coeff_id = try(scan(paste0("./gsmr_data/", trait1.name,"~",trait2.name, ".xmat.gz"), what="", nlines=1))
  snp_coeff = try(read.table(paste0("./gsmr_data/", trait1.name,"~",trait2.name, ".xmat.gz"), header=F, skip=2))
  if( inherits(snp_coeff_id, 'try-error'))  next
  if( inherits(snp_coeff, 'try-error'))  next

  snp_id = Reduce(intersect, list(clumped$SNP, snp_coeff_id))
  clumped = clumped[match(snp_id, clumped$SNP),]
  snp_order = match(snp_id, snp_coeff_id)
  snp_coeff_id = snp_coeff_id[snp_order]
  snp_coeff = snp_coeff[, snp_order]

  if(nrow(clumped)<2) next
  # Calculate the LD correlation matrix
  ldrho = cor(snp_coeff)
  #ldrho = diag(length(snp_coeff))

  # Check the size of the correlation matrix and double-check if the order of the SNPs in the LD correlation matrix is consistent with that in the GWAS summary data
  colnames(ldrho) = rownames(ldrho) = snp_coeff_id

  res = NULL

  for(Threshold in c(5e-05, 5e-06, 5e-07, 5e-08)){
    MRres = try(run_gsmr(clumped, exposure = trait1.name, outcome = trait2.name, Threshold, ldrho))
    if( inherits(MRres, 'try-error'))  next
    res = rbind(MRres, res)
  }

  write.table(res, "./res/GSMR_Rmodi.MRres", row.names = F, col.names = F, quote = F, append = T)
}





