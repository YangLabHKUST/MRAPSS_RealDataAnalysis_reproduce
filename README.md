# MRAPSS RealData analysis

The downlowad links for 20 GWAS summay data used in MR-APSS paper can be found at the "Table-GWAS_sources.xlsx".

The file "tested_pairs" are the trait pairs involved in analysis.

The R code for processing datasets and implementing other compared MR methods are in "R" folder.

Once you have downloaded the GWAS summay data, you can used the code "Format_GWAS20.R" to format the GWAS datasets, 

then you can use "est_paras_pairs.R" to estimate the nuisance parameters for all the trait pairs and do IV selection and LD clumping. 

The R code file begin with "run_" calls functions from R packages of MR methods to do MR tests.

The results are given in "Results" folder.
