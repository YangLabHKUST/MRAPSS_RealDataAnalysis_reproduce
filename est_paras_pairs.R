library(readr)
library(dplyr)
library(MRAPSS)


refPanel.dir = "/home/xianghonghu/1000G_EUR_Phase3_plink/all_1000G_EUR_Phase3"
ref = data.table::fread("/home/xianghonghu/1000G_EUR_Phase3_plink/all_1000G_EUR_Phase3.bim", header = F)
colnames(ref) = c("CHR","SNP", "V3","BP", "A1", "A2")


setwd("/home/share/xhhu/MR-APSS/reproduce")
pairs = read.table("/home/share/xhhu/MR-APSS/reproduce/tested_pairs",header = T)
colnames(pairs) = c("exposure", "outcome")
loc="/home/share/xhhu/MR-APSS/reproduce"
source("/home/share/xhhu/MR-APSS/Rcode/harmonise.R")
for(i in 1:nrow(pairs)){

   cat('Batch', 1 , "~", nrow(pairs), "pair", i, "\n")

  if(i > nrow(pairs)) break

  time=proc.time()

  trait1.name = as.character(pairs[i,"exposure"])
  trait2.name = as.character(pairs[i, "outcome"])

  trait1.dir = paste0("/home/share/xhhu/MR-APSS/reproduce/GWAS20_Formmated/", trait1.name)
  trait2.dir = paste0("/home/share/xhhu/MR-APSS/reproduce/GWAS20_Formmated/", trait2.name)


  cat(trait1.name,"~",trait2.name,"\n")
#
  trait1 <- readr::read_delim(trait1.dir,delim="\t",
                              escape_double = FALSE,
                              trim_ws = TRUE,
                              progress = F)

  trait2 = readr::read_delim(trait2.dir, "\t",
                             escape_double = FALSE,
                             trim_ws = TRUE,
                             progress = F)


# We have already formatted the datasets
  # trait1 = format_data(trait1,
  #                      snp_col="SNP",
  #                      b_col="b",
  #                      se_col="se",
  #                      freq_col="freq",
  #                      A1_col="A1",
  #                      A2_col="A2",
  #                      p_col="p",
  #                      n_col="n",
  #                      min_freq=0.05)

  # trait2 = format_data(trait2,
  #                   snp_col="SNP",
  #                   b_col="b",
  #                   se_col="se",
  #                   freq_col="freq",
  #                   A1_col="A1",
  #                   A2_col="A2",
  #                   p_col="p",
  #                   n_col="n",
  #                   min_freq=0.05)

  paras = est_paras(dat1 = trait1,
                    dat2 = trait2,
                    trait1.name = trait1.name,
                    trait2.name = trait2.name,
                    ldscore.dir = "~/eur_w_ld_chr")

  write.table(data.frame(trait1.name,
                         trait2.name,
                         nrow(paras$dat),
                         rg = paras$ldsc_res$rg,
                         rg.se = paras$ldsc_res$rg.se,
                         rho = paras$ldsc_res$I[1,2],
                         paras$ldsc_res$cov[1,1],  paras$ldsc_res$cov[1,2],  paras$ldsc_res$cov[2,2],
                         paras$ldsc_res$cov.se[1,1],  paras$ldsc_res$cov.se[1,2],  paras$ldsc_res$cov.se[2,2]),
              file = "ldsc.res", col.names = F, row.names = F, quote=F,append = T)

  if(inherits(paras, 'try-error')) next

  write.table(matrix(as.vector(paras$Omega),nrow=1), paste0("./Omega_est/", trait1.name, "~", trait2.name, "_Omega"), row.names = F,col.names = F, quote = F)
  write.table(matrix(as.vector(paras$Sigma_err), nrow=1), paste0("./Omega_est/", trait1.name, "~", trait2.name, "_Sigma_LD"), row.names = F,col.names = F, quote = F)

  cat("Begin clumping ...\n ")

  clumped =  clump(subset(paras$dat, pval.exp < 5e-04), SNP_col = "SNP",
                   pval_col = "pval.exp",
                   clump_kb = 1000,
                   clump_r2 = 0.001,
                   bfile = "/home/share/1000G_Phase3_all/all_1000G_EUR_Phase3",
                   plink_bin = "plink")

  if(inherits(clumped , 'try-error')) next

  # allele direction of GSMR and reference panle allele direction
  dat_ref = subset(ref, SNP %in% clumped$SNP)

  clumped = clumped[match(dat_ref$SNP, clumped$SNP),]

  clumped = harmonise_ref(clumped, dat_ref)

  write.table(clumped, file = paste0(loc, "/MRdata/", trait1.name,"~",trait2.name), col.names=T, row.names=F, quote=F)

  write.table(clumped[,c(1,2)], paste0(trait1.name,"~",trait2.name, "_snps.allele"), col.names=F, row.names=F, quote=F)

  # Extract the genotype data from a GWAS dataset using GCTA
  system(paste0("/home/xianghonghu/GCTA/gcta64 --bfile ", refPanel.dir, " --extract ", trait1.name,"~",trait2.name, "_snps.allele ", "--update-ref-allele ", trait1.name,"~",trait2.name, "_snps.allele ", "--recode --out ./gsmr_data/", trait1.name,"~",trait2.name))

  rm(clumped)

}




