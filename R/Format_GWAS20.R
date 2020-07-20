setwd("/home/share/xhhu/MR-APSS/reproduce/GWAS20")
library(MRAPSS)
#1 "AD.txt"
AD_raw <- readr::read_delim("AD.txt", "\t", escape_double = FALSE, trim_ws = TRUE, progress = F)
AD = format_data(AD_raw,
                 snp_col = "MarkerName",
                 b_col = "Beta",
                 se_col = "SE",
                 A1_col = "Effect_allele",
                 A2_col = "Non_Effect_allele",
                 p_col = "Pvalue",
                 n = 54162)

write.table(AD, "/home/share/xhhu/MR-APSS/reproduce/GWAS20_Formmated/AD",sep="\t", quote = F, row.names = F, col.names = T)

rm(list=c("AD","AD_raw"))

#2 "Alcohol.txt"
Alcohol_raw <- readr::read_delim("Alcohol.txt", "\t",
                                 escape_double = FALSE,
                                 trim_ws = TRUE, progress = F)

Alcohol = format_data(Alcohol_raw,
                      snp_col = "RSID",
                      b_col = "BETA",
                      se_col = "SE",
                      freq_col = "AF",
                      A1_col = "ALT",
                      A2_col = "REF",
                      p_col = "PVALUE",
                      n_col = "N")

write.table(Alcohol, "/home/share/xhhu/MR-APSS/reproduce/GWAS20_Formmated/Alcohol",sep="\t", quote = F, row.names = F, col.names = T)

rm(list=c("Alcohol_raw","Alcohol"))


##3 "Alcohol_ukb.txt"
Alcohol_ukb_raw <- readr::read_delim("Alcohol_ukb.txt", "\t",
                                     escape_double = FALSE,
                                     trim_ws = TRUE, progress = F)

Alcohol_ukb = format_data(Alcohol_ukb_raw,
                          snp_col = "MarkerName",
                          b_col = "Beta",
                          se_col = "SE",
                          freq_col = "EAF_A1",
                          A1_col = "A1",
                          A2_col = "A2",
                          p_col = "Pval",
                          n = 414343)

write.table(Alcohol_ukb, "/home/share/xhhu/MR-APSS/reproduce/GWAS20_Formmated/Alcohol_ukb",sep="\t", quote = F, row.names = F, col.names = T)

rm(list=c("Alcohol_ukb_raw","Alcohol_ukb"))

# 4 "AMD.txt"
AMD_raw <- readr::read_delim("AMD.txt", "\t",
                             escape_double = FALSE,
                             trim_ws = TRUE, progress = F)

AMD_raw$sign = ifelse(AMD_raw$Overall=="+", 1, -1)

AMD_raw$Z = AMD_raw$sign * abs(qnorm(AMD_raw$GC.Pvalue/2, 0, 1))

AMD = format_data(AMD_raw,
                  snp_col = "Marker",
                  A1_col = "Allele1",
                  A2_col = "Allele2",
                  z_col = "Z",
                  p_col = "GC.Pvalue",
                  ncase_col = "Ncases",
                  ncontrol_col = "Ncontrols")

write.table(AMD, "/home/share/xhhu/MR-APSS/reproduce/GWAS20_Formmated/AMD",sep="\t", quote = F, row.names = F, col.names = T)

rm(list=c("AMD_raw","AMD"))

# 5 "Angina.txt"
Angina_raw <- readr::read_delim("Angina.txt", "\t",
                                escape_double = FALSE,
                                trim_ws = TRUE, progress = F)

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

write.table(Angina, "/home/share/xhhu/MR-APSS/reproduce/GWAS20_Formmated/Angina",sep="\t", quote = F, row.names = F, col.names = T)
rm(list=c("Angina_raw","Angina"))


# 6 "Asthma.txt"
Asthma_raw <- readr::read_delim("Asthma.txt", "\t",
                                escape_double = FALSE,
                                trim_ws = TRUE, progress = F)

Asthma = format_data(Asthma_raw,
                     snp_col = "SNPID_UKB",
                     or_col = "OR",
                     se_col = "SE",
                     freq_col = "MAF_UKB",
                     A1_col = "A1",
                     A2_col = "A2",
                     p_col = "P",
                     n_col = "NMISS",
                     info_col = "INFO_UKB")

write.table(Asthma, "/home/share/xhhu/MR-APSS/reproduce/GWAS20_Formmated/Asthma",sep="\t", quote = F, row.names = F, col.names = T)
rm(list=c("Asthma_raw","Asthma"))


# 7. "BMI.txt"
BMI_raw <- readr::read_delim("BMI.txt", "\t", escape_double = FALSE,
                             trim_ws = TRUE, progress = F)
BMI = format_data(BMI_raw,
                  snp_col = "SNP",
                  b_col = "b",
                  se_col = "se",
                  freq_col = "Freq1.Hapmap",
                  A1_col = "A1",
                  A2_col = "A2",
                  p_col = "p",
                  n_col = "N")

write.table(BMI, "/home/share/xhhu/MR-APSS/reproduce/GWAS20_Formmated/BMI",sep="\t", quote = F, row.names = F, col.names = T)

rm(list=c("BMI_raw","BMI"))

# 8. "BMI_ukb.txt"
BMI_ukb_raw <- readr::read_delim("BMI_ukb.txt", "\t", escape_double = FALSE,
                                 trim_ws = TRUE, progress = F)
BMI_ukb = format_data(BMI_ukb_raw,
                      snp_col = "SNPID_UKB",
                      b_col = "BETA",
                      se_col = "SE",
                      freq_col = "MAF_UKB",
                      A1_col = "A1",
                      A2_col = "A2",
                      p_col = "P",
                      n_col = "NMISS",
                      info_col = "INFO_UKB")
write.table(BMI_ukb, "/home/share/xhhu/MR-APSS/reproduce/GWAS20_Formmated/BMI_ukb",sep="\t", quote = F, row.names = F, col.names = T)

rm(list=c("BMI_ukb_raw","BMI_ukb"))

# 9. "CAD.txt"
CAD_raw <- readr::read_delim("CAD.txt", "\t", escape_double = FALSE, trim_ws = TRUE, progress = F)
CAD = format_data(CAD_raw,
                  snp_col = "markername",
                  b_col = "beta",
                  se_col = "se",
                  freq_col = "effect_allele_freq",
                  A1_col = "effect_allele",
                  A2_col = "noneffect_allele",
                  p_col = "p_dgc",
                  n = 184305)
write.table(CAD, "/home/share/xhhu/MR-APSS/reproduce/GWAS20_Formmated/CAD",sep="\t", quote = F, row.names = F, col.names = T)
rm(list=c("CAD_raw","CAD"))

# 10. "Height.txt"
Height_raw <- readr::read_delim("Height.txt", "\t", escape_double = FALSE, trim_ws = TRUE, progress = F)
Height = format_data(Height_raw,
                     snp_col = "MarkerName",
                     b_col = "b",
                     se_col = "SE",
                     freq_col = "Freq.Allele1.HapMapCEU",
                     A1_col = "Allele1",
                     A2_col = "Allele2",
                     p_col = "p",
                     n_col = "N")

write.table(Height, "/home/share/xhhu/MR-APSS/reproduce/GWAS20_Formmated/Height",sep="\t", quote = F, row.names = F, col.names = T)
rm(list=c("Height_raw","Height"))

# 11.  "Height_ukb.txt"
Height_ukb_raw <- readr::read_delim("Height_ukb.txt", "\t", escape_double = FALSE,
                                    trim_ws = TRUE, progress = F)
Height_ukb = format_data(Height_ukb_raw,
                         snp_col = "SNPID_UKB",
                         b_col = "BETA",
                         se_col = "SE",
                         freq_col = "MAF_UKB",
                         A1_col = "A1",
                         A2_col = "A2",
                         p_col = "P",
                         n_col = "NMISS",
                         info_col = "INFO_UKB")
write.table(Height_ukb, "/home/share/xhhu/MR-APSS/reproduce/GWAS20_Formmated/Height_ukb",sep="\t", quote = F, row.names = F, col.names = T)
rm(list=c("Height_ukb_raw","Height_ukb"))


# 12. "Income.txt"
Income_raw <- readr::read_delim("Income.txt", " ", escape_double = FALSE, trim_ws = TRUE, progress = F)
Income = format_data(Income_raw,
                     snp_col = "SNP",
                     b_col = "Beta",
                     se_col = "Standard_Error_of_Beta",
                     A1_col = "Effect_Allele",
                     A2_col = "Non_effect_Allele",
                     p_col = "P",
                     n=286301)
write.table(Income, "/home/share/xhhu/MR-APSS/reproduce/GWAS20_Formmated/Income",sep="\t", quote = F, row.names = F, col.names = T)

rm(list=c("Income_raw","Income"))

# 13.
Intelligence_raw <- readr::read_delim("Intelligence.txt", "\t", escape_double = FALSE, trim_ws = TRUE, progress = F)
Intelligence = format_data(Intelligence_raw,
                           snp_col = "rsid",
                           b_col = "Beta",
                           se_col = "SE",
                           A1_col = "ref",
                           A2_col = "alt",
                           p_col = "p_value",
                           n=78308)
write.table(Intelligence, "/home/share/xhhu/MR-APSS/reproduce/GWAS20_Formmated/Intelligence",sep="\t", quote = F, row.names = F, col.names = T)
rm(list=c("Intelligence_raw","Intelligence"))

# 14.  "Intelligence_ukb.txt"
Intelligence_ukb_raw <- readr::read_delim("Intelligence_ukb.txt", "\t", escape_double = FALSE, trim_ws = TRUE, progress = F)
Intelligence_ukb = format_data(Intelligence_ukb_raw,
                               snp_col = "SNPID_UKB",
                               b_col = "BETA",
                               se_col = "SE",
                               freq_col = "MAF_UKB",
                               A1_col = "A1",
                               A2_col = "A2",
                               p_col = "P",
                               n_col = "NMISS",
                               info_col = "INFO_UKB")

write.table(Intelligence_ukb, "/home/share/xhhu/MR-APSS/reproduce/GWAS20_Formmated/Intelligence_ukb",sep="\t", quote = F, row.names = F, col.names = T)

rm(list=c("Intelligence_ukb_raw","Intelligence_ukb"))

# 15. "MouthUlcers.txt"
MouthUlcers_raw <- readr::read_delim("MouthUlcers.txt", "\t", escape_double = FALSE, trim_ws = TRUE, progress = F)

MouthUlcers = format_data(MouthUlcers_raw,
                          snp_col = "SNP",
                          b_col = "LOG_OR",
                          se_col = "LOG_SE",
                          freq_col = "A1FREQ",
                          A1_col = "ALLELE1",
                          A2_col = "ALLELE0",
                          p_col = "P_BOLT_LMM_INF",
                          ncase_col = "N_CASES",
                          ncontrol_col = "N_CONTROLS",
                          info_col = "INFO")

write.table(MouthUlcers, "/home/share/xhhu/MR-APSS/reproduce/GWAS20_Formmated/MouthUlcers",sep="\t", quote = F, row.names = F, col.names = T)
rm(list=c("MouthUlcers_raw","MouthUlcers"))

# 16. "Never_Smoking.txt"
Never_Smoking_ukb_raw <- readr::read_delim("Never_Smoking_ukb.txt", "\t", escape_double = FALSE, trim_ws = TRUE, progress = F)

Never_Smoking_ukb = format_data(Never_Smoking_ukb_raw,
                            snp_col = "SNPID_UKB",
                            or_col = "OR",
                            se_col = "SE",
                            freq_col = "MAF_UKB",
                            A1_col = "A1",
                            A2_col = "A2",
                            p_col = "P",
                            n_col = "NMISS",
                            info_col = "INFO_UKB")

write.table(Never_Smoking_ukb, "/home/share/xhhu/MR-APSS/reproduce/GWAS20_Formmated/Never_Smoking_ukb",sep="\t", quote = F, row.names = F, col.names = T)
rm(list=c("Never_Smoking_ukb_raw","Never_Smoking_ukb"))

# 17. "SCZ.txt"
#awk -F '\t' 'NR == FNR {f2[$1,$4] = $2; next}($3,$4) in f2 {print f2[$3,$4]"\t"$0}' /home/share/1000G_Phase3_all/all_1000G_EUR_Phase3.bim SCZ.txt >  SCZ1.txt
#sed  -i '1i\rsid\tSNP\tFreq.A1\tCHR\tBP\tA1\tA2\tOR\tSE\tP' SCZ1.txt

SCZ_raw <- readr::read_delim("SCZ1.txt", "\t", escape_double = FALSE, trim_ws = TRUE, progress = F)

SCZ = format_data(SCZ_raw,
                  snp_col = "rsid",
                  or_col = "OR",
                  se_col = "SE",
                  freq_col = "Freq.A1",
                  A1_col = "A1",
                  A2_col = "A2",
                  p_col = "P",
                  n=105318)

write.table(SCZ, "/home/share/xhhu/MR-APSS/reproduce/GWAS20_Formmated/SCZ",sep="\t", quote = F, row.names = F, col.names = T)

rm(list=c("SCZ_raw","SCZ"))

#18. "Smoking.txt"
Smoking_raw <- readr::read_delim("Smoking.txt", "\t", escape_double = FALSE, trim_ws = TRUE, progress = F)

Smoking = format_data(Smoking_raw,
                  snp_col = "RSID",
                  b_col = "BETA",
                  se_col = "SE",
                  freq_col = "AF",
                  A1_col = "ALT",
                  A2_col = "REF",
                  p_col = "PVALUE",
                  n_col = "N")

write.table(Smoking , "/home/share/xhhu/MR-APSS/reproduce/GWAS20_Formmated/Smoking",sep="\t", quote = F, row.names = F, col.names = T)
rm(list=c("Smoking_raw","Smoking"))

#19. "T2D.txt"
T2D_raw <- readr::read_delim("T2D.txt", " ", escape_double = FALSE, trim_ws = TRUE, progress = F)
T2D = format_data(T2D_raw,
                  snp_col = "SNP",
                  b_col = "b",
                  se_col = "se",
                  freq_col = "frq_A1",
                  A1_col = "A1",
                  A2_col = "A2",
                  p_col = "P",
                  n_col = "N")
write.table(T2D, "/home/share/xhhu/MR-APSS/reproduce/GWAS20_Formmated/T2D",sep="\t", quote = F, row.names = F, col.names = T)
rm(list=c("T2D_raw","T2D"))

#20 "Urate.csv"
Urate_raw <- readr::read_delim("Urate.csv", ",", escape_double = FALSE, trim_ws = TRUE, progress = F)
Urate = format_data(Urate_raw,
                    snp_col = "MarkerName",
                    b_col = "beta",
                    se_col = "se",
                    A1_col = "A1",
                    A2_col = "A2",
                    p_col = "p_gc",
                    n_col = "n_total")
write.table(Urate, "/home/share/xhhu/MR-APSS/reproduce/GWAS20_Formmated/Urate",sep="\t", quote = F, row.names = F, col.names = T)
rm(list=c("Urate_raw","Urate"))

exp = c("Alcohol","Alcohol_ukb","BMI","BMI_ukb","Height","Height_ukb","Intelligence","Intelligence_ukb","Smoking","Never_Smoking_ukb")
out = c("AMD","Asthma","Angina","T2D","Urate","SCZ","MouthUlcers","CAD","AD","Income")
pairs = merge(exp,out)
colnames(pairs) = c("exposure","outcome")
write.table(pairs, "/home/share/xhhu/MR-APSS/reproduce/tested_pairs", row.names = F, col.names = T, quote = F)

