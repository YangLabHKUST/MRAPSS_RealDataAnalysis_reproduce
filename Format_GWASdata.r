library(MRAPSS)
library(readr)

## 1. AD
AD_raw <- readr::read_delim("./GWAS_26and5_raw/AD", "\t", escape_double = FALSE, trim_ws = TRUE, progress = F)
AD = format_data(AD_raw,
                 snp_col = "MarkerName",
                 b_col = "Beta",
                 se_col = "SE",
                 A1_col = "Effect_allele",
                 A2_col = "Non_Effect_allele",
                 p_col = "Pvalue",
                 n = 54162)

write.table(AD, "./GWAS_26and5_formatted/AD", sep="\t", quote = F, row.names = F, col.names = T)

rm(list=c("AD","AD_raw"))

## 2. Alcohol
Alcohol_raw <- readr::read_delim("./GWAS_26and5_raw/Alcohol", "\t", escape_double = FALSE, trim_ws = TRUE, progress = F)

Alcohol = format_data(Alcohol_raw,
                      snp_col = "MarkerName",
                      b_col = "Beta",
                      se_col = "SE",
                      freq_col = "EAF_A1",
                      A1_col = "A1",
                      A2_col = "A2",
                      p_col = "Pval",
                      n = 414343)

write.table(Alcohol, "./GWAS_26and5_formatted/Alcohol", sep="\t", quote = F, row.names = F, col.names = T)

rm(list=c("Alcohol_raw","Alcohol"))

## 3. Angina
Angina_raw <- readr::read_delim("./GWAS_26and5_raw/Angina", "\t", escape_double = FALSE, trim_ws = TRUE, progress = F)

Angina = format_data(Angina_raw,
                     snp_col = "SNPID_UKB",
                     or_col = "OR",
                     se_col = "SE",
                     freq_col = "MAF_UKB",
                     A1_col = "A1",
                     A2_col = "A2",
                     p_col = "P",
                     n_col = "NMISS",
                     info_col = "INFO_UKB")

write.table(Angina, "./GWAS_26and5_formatted/Angina", sep="\t", quote = F, row.names = F, col.names = T)

rm(list=c("Angina_raw","Angina"))

## 4. Anorexia
Anorexia_raw <- readr::read_delim("./GWAS_26and5_raw/Anorexia", "\t", escape_double = FALSE, trim_ws = TRUE, progress = F)

Anorexia = format_data(Anorexia_raw,
                       snp_col = "ID",
                       b_col = "BETA",
                       se_col = "SE",
                       A1_col = "ALT",
                       A2_col = "REF",
                       p_col = "PVAL",
                       ncase_col  = "NCAS",
                       ncontrol_col = "NCON",
                       info_col = "IMPINFO")

write.table(Anorexia, "./GWAS_26and5_formatted/Anorexia", sep = "\t", quote = F, row.names = F, col.names = T)

rm(list=c("Anorexia_raw","Anorexia"))


## 5. ASD
ASD_raw <- readr::read_delim("./GWAS_26and5_raw/ASD", "\t", escape_double = FALSE, trim_ws = TRUE, progress = F)

ASD = format_data(ASD_raw,
                  snp_col = "snp",
                  b_col = "b",
                  se_col = "StdErr",
                  freq_col = "freq_A1",
                  A1_col = "A1",
                  A2_col = "A2",
                  p_col = "p",
                  n_col = "N")

write.table(ASD, "./GWAS_26and5_formatted/ASD", sep = "\t", quote = F, row.names = F, col.names = T)

rm(list=c("ASD_raw","ASD"))

## 6. BMI
BMI_raw <- readr::read_delim("./GWAS_26and5_raw/BMI", "\t", escape_double = FALSE, trim_ws = TRUE, progress = F)

BMI = format_data(BMI_raw,
                  snp_col = "SNPID_UKB",
                  b_col = "BETA",
                  se_col = "SE",
                  freq_col = "MAF_UKB",
                  A1_col = "A1",
                  A2_col = "A2",
                  p_col = "P",
                  n_col = "NMISS",
                  info_col = "INFO_UKB")

write.table(BMI, "./GWAS_26and5_formatted/BMI", sep = "\t", quote = F, row.names = F, col.names = T)

rm(list=c("BMI_raw","BMI"))

## 7. CAD
CAD_raw <- readr::read_delim("./GWAS_26and5_raw/CAD", "\t", escape_double = FALSE, trim_ws = TRUE, progress = F)

CAD = format_data(CAD_raw,
                  snp_col = "markername",
                  b_col = "beta",
                  se_col = "se",
                  freq_col = "effect_allele_freq",
                  A1_col = "effect_allele",
                  A2_col = "noneffect_allele",
                  p_col = "p_dgc",
                  n = 184305)

write.table(CAD, "./GWAS_26and5_formatted/CAD", sep = "\t", quote = F, row.names = F, col.names = T)

rm(list=c("CAD_raw","CAD"))

## 8. CD (Crohn Disease)
# CD_rsid: CD annotated with rs-number and freq

CD_raw <- readr::read_delim("./GWAS_26and5_raw/CD_rsid", "\t", escape_double = FALSE, trim_ws = TRUE, progress = F)

CD = format_data(CD_raw,
                 snp_col = "rsid",
                 b_col = "Effect",
                 se_col = "StdErr",
                 freq_col = "freq",
                 A1_col = "Allele1",
                 A2_col = "Allele2",
                 p_col = "P.value",
                 n = 40266)

write.table(CD, "./GWAS_26and5_formatted/CD", sep = "\t", quote = F, row.names = F, col.names = T)

rm(list=c("CD_raw","CD"))

## 9. DaytimeSleepiness
dat_raw <- readr::read_delim("./GWAS_26and5_raw/DaytimeSleepiness", "\t", escape_double = FALSE, trim_ws = TRUE, progress = F)

dat = format_data(dat_raw,
                  snp_col = "SNP",
                  b_col = "BETA",
                  se_col = "SE",
                  freq_col = "A1FREQ",
                  A1_col = "ALLELE1",
                  A2_col = "ALLELE0",
                  info_col = "INFO",
                  p_col = "P",
                  n = 452071)

write.table(dat, "./GWAS_26and5_formatted/Daytime_Sleepiness",  sep="\t", quote = F, row.names = F, col.names = T)

rm(list=c("dat_raw","dat"))

## 10. Depression
dat_raw <- readr::read_delim("./GWAS_26and5_raw/Depression", "\t", escape_double = FALSE, trim_ws = TRUE, progress = F)

dat = MRAPSS::format_data(dat_raw,
                          snp_col = "RSID",
                          z_col = "Z",
                          freq_col = "MAF_UKB",
                          A1_col = "A1",
                          A2_col = "A2",
                          p_col = "P",
                          n_col = "N",
                          info_col = "INFO_UKB")

write.table(dat, "./GWAS_26and5_formatted/Depression", sep = "\t", quote = F, row.names = F, col.names = T)

rm(list=c("dat_raw","dat"))

## 11-14. Hair_Black, Hair_Blonde, Hair_Dark_Brown, Hair_Light_Brown
for (trait in c("Hair_Black", "Hair_Blonde", "Hair_Dark_Brown", "Hair_Light_Brown")){
  
  dat_raw <- readr::read_delim(paste0("./GWAS_26and5_raw/", trait), delim = "\t", escape_double = FALSE, trim_ws = TRUE, progress = F)
  
  dat = format_data(dat_raw,
                    snp_col = "SNPID_UKB",
                    or_col = "OR",
                    se_col = "SE",
                    freq_col = "MAF_UKB",
                    A1_col = "A1",
                    A2_col = "A2",
                    p_col = "P",
                    n_col = "NMISS",
                    info_col = "INFO_UKB")
  
  write.table(dat, paste0("./GWAS_26and5_formatted/", trait), sep="\t", quote = F, row.names = F, col.names = T)
  
  rm(list=c("dat_raw", "dat"))
  
}


## 15. HBP (High Blood Pressure)

HBP_raw <- readr::read_delim("./GWAS_26and5_raw/HBP", "\t", escape_double = FALSE, trim_ws = TRUE, progress = F)

HBP = MRAPSS::format_data(HBP_raw,
                          snp_col = "SNPID_UKB",
                          or_col = "OR",
                          se_col = "SE",
                          freq_col = "MAF_UKB",
                          A1_col = "A1",
                          A2_col = "A2",
                          p_col = "P",
                          n_col = "NMISS",
                          info_col = "INFO_UKB")

write.table(HBP, "./GWAS_26and5_formatted/HBP", sep="\t", quote = F, row.names = F, col.names = T)

rm(list=c("HBP_raw", "HBP"))

## 16. Height (GIANT)

Height_raw <- readr::read_delim("./GWAS_26and5_raw/Height_GIANT", "\t", escape_double = FALSE, trim_ws = TRUE, progress = F)
Height_raw = Height_raw[!is.na(Height_raw$Freq.Allele1.HapMapCEU), ]
Height = format_data(Height_raw,
                     snp_col = "MarkerName",
                     b_col = "b",
                     se_col = "SE",
                     freq_col = "Freq.Allele1.HapMapCEU",
                     A1_col = "Allele1",
                     A2_col = "Allele2",
                     p_col = "p",
                     n_col = "N")

write.table(Height, "./GWAS_26and5_formatted/Height_GIANT", sep="\t", quote = F, row.names = F, col.names = T)

rm(list=c("Height_raw","Height"))

## 17. Height (UKBB)

Height_UKB_raw <- readr::read_delim("./GWAS_26and5_raw/Height_UKB", "\t", escape_double = FALSE,
                                    trim_ws = TRUE, progress = F)

Height_UKB = format_data(Height_UKB_raw,
                         snp_col = "SNPID_UKB",
                         b_col = "BETA",
                         se_col = "SE",
                         freq_col = "MAF_UKB",
                         A1_col = "A1",
                         A2_col = "A2",
                         p_col = "P",
                         n_col = "NMISS",
                         info_col = "INFO_UKB")

write.table(Height_UKB, "./GWAS_26and5_formatted/Height_UKB", sep = "\t", quote = F, row.names = F, col.names = T)

rm(list=c("Height_UKB_raw","Height_UKB"))

# 18. IBD (Inflammatory Bowel Disease)

IBD_raw <- readr::read_delim("./GWAS_26and5_raw/IBD_rsid", "\t", escape_double = FALSE, trim_ws = TRUE, progress = F)

IBD = format_data(IBD_raw,
                  snp_col = "rsid",
                  b_col = "Effect",
                  se_col = "StdErr",
                  freq_col = "freq",
                  A1_col = "Allele1",
                  A2_col = "Allele2",
                  p_col = "P.value",
                  n = 59957)

write.table(IBD, "./GWAS_26and5_formatted/IBD", sep="\t", quote = F, row.names = F, col.names = T)

rm(list=c("IBD_raw","IBD"))

## 19. Income

Income_raw <- readr::read_delim("./GWAS_26and5_raw/Income", " ", escape_double = FALSE, trim_ws = TRUE, progress = F)

Income = format_data(Income_raw,
                     snp_col = "SNP",
                     b_col = "Beta",
                     se_col = "Standard_Error_of_Beta",
                     A1_col = "Effect_Allele",
                     A2_col = "Non_effect_Allele",
                     p_col = "P",
                     n=286301)

write.table(Income, "./GWAS_26and5_formatted/Income", sep = "\t", quote = F, row.names = F, col.names = T)

rm(list=c("Income_raw","Income"))

## 20. Insomnia

Insomnia_raw <- readr::read_delim("./GWAS_26and5_raw/Insomnia", " ", escape_double = FALSE, trim_ws = TRUE, progress = F)

Insomnia = format_data(Insomnia_raw,
                       snp_col = "SNP",
                       b_col = "BETA_INSOMNIA",
                       se_col = "SE_INSOMNIA",
                       freq_col = "A1FREQ",
                       A1_col = "ALLELE1",
                       A2_col = "ALLELE0",
                       info_col = "INFO",
                       p_col = "P_INSOMNIA",
                       n = 453379)

write.table(Insomnia, "./GWAS_26and5_formatted/Insomnia", sep = "\t", quote = F, row.names = F, col.names = T)

rm(list=c("Insomnia_raw","Insomnia"))

## 21. Intelligence

Intelligence_raw <- readr::read_delim("./GWAS_26and5_raw/Intelligence", "\t", escape_double = FALSE, trim_ws = TRUE, progress = F)

Intelligence = format_data(Intelligence_raw,
                           snp_col = "SNPID_UKB",
                           b_col = "BETA",
                           se_col = "SE",
                           freq_col = "MAF_UKB",
                           A1_col = "A1",
                           A2_col = "A2",
                           p_col = "P",
                           n_col = "NMISS",
                           info_col = "INFO_UKB")

write.table(Intelligence, "./GWAS_26and5_formatted/Intelligence", sep = "\t", quote = F, row.names = F, col.names = T)

rm(list=c("Intelligence_raw","Intelligence"))

## 22. MDD
MDD_raw <- readr::read_delim("./GWAS_26and5_raw/MDD", "\t", escape_double = FALSE, trim_ws = TRUE, progress = F)

MDD = MRAPSS::format_data(MDD_raw,
                          snp_col = "SNPID_UKB",
                          or_col = "OR",
                          se_col = "SE",
                          freq_col = "MAF_UKB",
                          A1_col = "A1",
                          A2_col = "A2",
                          p_col = "P",
                          n_col = "NMISS",
                          info_col = "INFO_UKB")

write.table(MDD, "./GWAS_26and5_formatted/MDD", sep = "\t", quote = F, row.names = F, col.names = T)

rm(list=c("MDD_raw","MDD"))

## 23. NEB
NEB_raw <- readr::read_delim("./GWAS_26and5_raw/NEB", "\t", escape_double = FALSE, trim_ws = TRUE, progress = F)

NEB = MRAPSS::format_data(NEB_raw,
                          snp_col = "SNPID",
                          freq_col = "Freq_HapMap",
                          A1_col = "A1",
                          A2_col = "A2",
                          z_col = "Zscore",
                          p_col = "Pvalue",
                          n = 343072)

write.table(NEB, file = "./GWAS_26and5_formatted/NEB", sep = "\t", quote = F, row.names = F, col.names = T)

rm(list=c("NEB_raw","NEB"))

## 24. Neuroticism
dat_raw <- readr::read_delim("./GWAS_26and5_raw/Neuroticism", "\t", escape_double = FALSE, trim_ws = TRUE, progress = F)

dat = MRAPSS::format_data(dat_raw,
                          snp_col = "SNPID_UKB",
                          b_col = "BETA",
                          se_col = "SE",
                          freq_col = "MAF_UKB",
                          A1_col = "A1",
                          A2_col = "A2",
                          p_col = "P",
                          n_col = "NMISS",
                          info_col = "INFO_UKB")

write.table(dat, "./GWAS_26and5_formatted/Neuroticism", sep = "\t", quote = F, row.names = F, col.names = T)

rm(list=c("dat_raw","dat"))

## 25. RA
# RA_freq: RA annotated with a column of freq
# Note that the freq_col is not a required column. If freq_col is not available, it will skip the step of QC with freq.

RA_raw <- readr::read_delim("./GWAS_26and5_raw/RA_freq", "\t", escape_double = FALSE, trim_ws = TRUE, progress = F)

RA = format_data(RA_raw,
                  snp_col = "SNPID",
                  or_col = "OR(A1)",
                  A1_col = "A1",
                  A2_col = "A2",
                  freq_col = "freq",
                  p_col = "P-val",
                  n = 58284)

write.table(RA, "./GWAS_26and5_formatted/RA", sep = "\t", quote = F, row.names = F, col.names = T)

rm(list=c("RA_raw","RA"))


## 26. SCZ

SCZ_raw <- readr::read_delim("./GWAS_26and5_raw/SCZ_rsid", "\t", escape_double = FALSE, trim_ws = TRUE, progress = F)

SCZ = format_data(SCZ_raw,
                  snp_col = "rsid",
                  or_col = "OR",
                  se_col = "SE",
                  freq_col = "Freq.A1",
                  A1_col = "A1",
                  A2_col = "A2",
                  p_col = "P",
                  n=105318)

write.table(SCZ, "./GWAS_26and5_formatted/SCZ", sep = "\t", quote = F, row.names = F, col.names = T)

rm(list=c("SCZ_raw","SCZ"))

## 27. Smoking
Smoking_raw <- readr::read_delim("./GWAS_26and5_raw/Smoking", "\t", escape_double = FALSE, trim_ws = TRUE, progress = F)
Smoking_raw = Smoking_raw[!is.na(Smoking_raw$AF), ]
Smoking = format_data(Smoking_raw,
                      snp_col = "RSID",
                      b_col = "BETA",
                      se_col = "SE",
                      freq_col = "AF",
                      A1_col = "ALT",
                      A2_col = "REF",
                      p_col = "PVALUE",
                      n_col = "N")

write.table(Smoking , "./GWAS_26and5_formatted/Smoking", sep = "\t", quote = F, row.names = F, col.names = T)

rm(list=c("Smoking_raw","Smoking"))


## 28. SWB

SWB_raw <- readr::read_delim("./GWAS_26and5_raw/SWB", "\t", escape_double = FALSE, trim_ws = TRUE, progress = F)

SWB = MRAPSS::format_data(SWB_raw,
                          snp_col = "MarkerName",
                          freq_col = "EAF",
                          A1_col = "A1",
                          A2_col = "A2",
                          b_col = "Beta",
                          se_col = "SE",
                          p_col = "Pval",
                          n = 298420)

write.table(SWB, file = "./GWAS_26and5_formatted/SWB", sep = "\t", quote = F, row.names = F, col.names = T)

rm(list=c("SWB_raw","SWB"))

## 29. T2D

T2D_raw <- readr::read_delim("./GWAS_26and5_raw/T2D_rsid", "\t", escape_double = FALSE, trim_ws = TRUE, progress = F)

T2D = MRAPSS::format_data(T2D_raw,
                          snp_col = "rsid",
                          b_col = "Beta",
                          se_col = "SE",
                          freq_col = "EAF",
                          A1_col = "EA",
                          A2_col = "NEA",
                          p_col = "Pvalue",
                          n_col = "Neff")

write.table(T2D, "./GWAS_26and5_formatted/T2D", sep = "\t", quote = F, row.names = F, col.names = T)

rm(list=c("T2D_raw","T2D"))

## 30. Tanning

dat_raw <- readr::read_delim("./GWAS_26and5_raw/Tanning", delim = "\t", escape_double = FALSE, trim_ws = TRUE, progress = F)

dat = MRAPSS::format_data(dat_raw,
                          snp_col = "SNPID_UKB",
                          b_col = "BETA",
                          se_col = "SE",
                          freq_col = "MAF_UKB",
                          A1_col = "A1",
                          A2_col = "A2",
                          p_col = "P",
                          n_col = "NMISS",
                          info_col = "INFO_UKB")

write.table(dat, "./GWAS_26and5_formatted/Tanning", sep = "\t", quote = F, row.names = F, col.names = T)

rm(list=c("dat_raw", "dat"))

## 31. Urate

Urate_raw <- readr::read_delim("./GWAS_26and5_raw/Urate", ",", escape_double = FALSE, trim_ws = TRUE, progress = F)

Urate = format_data(Urate_raw,
                    snp_col = "MarkerName",
                    b_col = "beta",
                    se_col = "se",
                    A1_col = "A1",
                    A2_col = "A2",
                    p_col = "p_gc",
                    n_col = "n_total")

write.table(Urate, "./GWAS_26and5_formatted/Urate", sep = "\t", quote = F, row.names = F, col.names = T)

rm(list=c("Urate_raw","Urate"))
