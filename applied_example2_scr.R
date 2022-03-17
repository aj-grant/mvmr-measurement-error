################################################################################
#Script used for the second applied example
################################################################################

library(MendelianRandomization)
library(TwoSampleMR)
library(tidyverse)

set.seed(20200127)
#Load datasets
ST8 = read.csv('ST8.csv', header = T)
ST4 = read.csv('ST4.csv', header = T)
ST16 = read.csv('ST16.csv', header = T)
ST13 = read.csv('ST13.csv', header = T)

#BMI
bmi_all = read.table('BMI Cluster/Datasets downloaded/bmi_all_snp.txt', header = T,
                     colClasses = c("character", "numeric", "numeric", "character", "character", "numeric", "numeric", "numeric"))

#Smoking
smo_all = read.table('BMI Cluster/Datasets downloaded/smo_edu_out.txt', header = T,
                     colClasses = c("character", "numeric", "numeric", "character", "numeric", "numeric", "numeric"))

#SBP
sbp_all = read.table('BMI Cluster/Datasets downloaded/sbp_edu_out.txt', header = T,
                     colClasses = c("character", "numeric", "numeric", "numeric"))
var_rsid_map = read.table(file = 'BMI Cluster/Datasets downloaded/var_rsid_map2.txt',
                          colClasses = c("character", "character", "character"),  header = TRUE, row.names = NULL)
names(var_rsid_map) = c("variant", "SNP", "alt")
sbp_all = inner_join(sbp_all, var_rsid_map, by = "variant")

ST8$p = 2 * (1 - pnorm(abs(ST8$gx/ST8$gx_se)))
ST8_clump = clump_data(format_data(ST8, snp_col = "SNP", beta_col = "gx", se_col = "gx_se", effect_allele_col = "ea",
                                   pval_col = "p"), clump_r2 = 0.001)
ST8_r001 = data.frame("SNP" = ST8_clump$SNP)
ST8_r001 = inner_join(ST8_r001, ST8, by = "SNP")
edu_chd_est_clump = MendelianRandomization::mr_ivw(mr_input(bx = ST8_r001$gx, bxse = ST8_r001$gx_se, by = ST8_r001$gy, byse = ST8_r001$gy_se))
edu_chd_est_me_clump = mrest_me(mr_input(bx = ST8_r001$gx, bxse = ST8_r001$gx_se, by = ST8_r001$gy, byse = ST8_r001$gy_se))
c(edu_chd_est_clump$Estimate, edu_chd_est_clump$CILower, edu_chd_est_clump$CIUpper)
c(edu_chd_est_me_clump$thest, edu_chd_est_me_clump$thest - qnorm(0.975) *sqrt(edu_chd_est_me_clump$Var[1, 1]),
  edu_chd_est_me_clump$thest + qnorm(0.975) *sqrt(edu_chd_est_me_clump$Var[1, 1]))

#MV Analysis
edu_snps = data.frame("SNP" = ST8$SNP, "EA" = ST8$ea, "edu_beta" = ST8$gx, "edu_se" = ST8$gx_se)
edu_snps$p = 2 * (1 - pnorm(abs(edu_snps$edu_beta / edu_snps$edu_se)))
bmi_snps = data.frame("SNP" = ST4$SNP, "EA" = ST4$ea, "edu_beta" = ST4$gz, "edu_se" = ST4$gz_se)
bmi_snps$p = 2 * (1 - pnorm(abs(bmi_snps$edu_beta / bmi_snps$edu_se)))
smo_snps = data.frame("SNP" = ST16$SNP, "EA" = ST16$ea, "edu_beta" = ST16$gz, "edu_se" = ST16$gz_se)
smo_snps$p = 2 * (1 - pnorm(abs(smo_snps$edu_beta / smo_snps$edu_se)))
sbp_snps = data.frame("SNP" = ST13$SNP, "EA" = ST13$ea, "edu_beta" = ST13$gz, "edu_se" = ST13$gz_se)
sbp_snps$p = 2 * (1 - pnorm(abs(sbp_snps$edu_beta / sbp_snps$edu_se)))
all_snps = rbind(edu_snps, bmi_snps, smo_snps, sbp_snps)
all_fm = format_data(all_snps, snp_col = "SNP", beta_col = "edu_beta", se_col = "edu_se", effect_allele_col = "EA",
                     pval_col = "p")
all_clump = clump_data(all_fm, clump_r2 = 0.001)
all_dat = all_clump[, c("SNP", "effect_allele.exposure", "beta.exposure", "se.exposure")]
names(all_dat) = c("SNP", "EA", "Beta_edu", "SE_edu")

all_dat = inner_join(all_dat, bmi_all[, c("SNP", "Tested_Allele", "BETA", "SE")], by = "SNP")
names(all_dat)[5:7] = c("EA_bmi", "Beta_bmi", "SE_bmi")
all_dat = inner_join(all_dat, smo_all[, c("SNP", "EFFECT_ALLELE", "BETA", "SE")], by = "SNP")
names(all_dat)[8:10] = c("EA_smo", "Beta_smo", "SE_smo")
all_dat = inner_join(all_dat, sbp_all[, c("SNP", "alt", "beta", "se")], by = "SNP")
names(all_dat)[11:13] = c("EA_sbp", "Beta_sbp", "SE_sbp")

all_dat = all_dat[!duplicated(all_dat$SNP), ]

chd_all = rbind(data.frame("SNP" = ST8$SNP, "EA_chd" = ST8$ea, "Beta_chd" = ST8$gy, "SE_chd" = ST8$gy_se),
                data.frame("SNP" = ST4$SNP, "EA_chd" = ST4$ea, "Beta_chd" = ST4$gy, "SE_chd" = ST4$gy_se),
                data.frame("SNP" = ST16$SNP, "EA_chd" = ST16$ea, "Beta_chd" = ST16$gy, "SE_chd" = ST16$gy_se),
                data.frame("SNP" = ST13$SNP, "EA_chd" = ST13$ea, "Beta_chd" = ST13$gy, "SE_chd" = ST13$gy_se))

all_dat = inner_join(all_dat, chd_all[!duplicated(chd_all$SNP), ], by = "SNP")

all_dat$Beta_bmi = all_dat$Beta_bmi * (2*as.numeric(all_dat$EA==all_dat$EA_bmi)-1)
all_dat$Beta_smo = all_dat$Beta_smo * (2*as.numeric(all_dat$EA==all_dat$EA_smo)-1)
all_dat$Beta_sbp = all_dat$Beta_sbp * (2*as.numeric(all_dat$EA==all_dat$EA_sbp)-1)
all_dat$Beta_chd = all_dat$Beta_chd * (2*as.numeric(all_dat$EA==all_dat$EA_chd)-1)

all_dat = all_dat[, c("SNP", "EA", "Beta_edu", "SE_edu", "Beta_bmi", "SE_bmi", "Beta_smo", "SE_smo", "Beta_sbp",
                      "SE_sbp", "Beta_chd", "SE_chd")]

edu_chd_est_mv = mr_mvivw(mr_mvinput(bx = cbind(all_dat$Beta_edu, all_dat$Beta_bmi, all_dat$Beta_smo, all_dat$Beta_sbp),
                                     bxse = cbind(all_dat$SE_edu, all_dat$SE_bmi, all_dat$SE_smo, all_dat$SE_sbp),
                                     by = all_dat$Beta_chd, byse = all_dat$SE_chd))
edu_chd_est_mv_me = mrest_me(mr_mvinput(bx = cbind(all_dat$Beta_edu, all_dat$Beta_bmi, all_dat$Beta_smo, all_dat$Beta_sbp),
                                        bxse = cbind(all_dat$SE_edu, all_dat$SE_bmi, all_dat$SE_smo, all_dat$SE_sbp),
                                        by = all_dat$Beta_chd, byse = all_dat$SE_chd))
c(edu_chd_est_mv$Estimate[1], edu_chd_est_mv$CILower[1], edu_chd_est_mv$CIUpper[1])
c(edu_chd_est_mv_me$thest[1], edu_chd_est_mv_me$thest[1] - qnorm(0.975) * sqrt(edu_chd_est_mv_me$Var[1, 1]),
  edu_chd_est_mv_me$thest[1] + qnorm(0.975) * sqrt(edu_chd_est_mv_me$Var[1, 1]))

1 - edu_chd_est_mv$Estimate[1] / edu_chd_est_clump$Estimate
1 - edu_chd_est_mv_me$thest[1] / edu_chd_est_me_clump$thest

#CIs with delta method
Est_ivw = 1 - edu_chd_est_mv$Estimate[1] / edu_chd_est_clump$Estimate
rho1 = 0
SE_ivw = sqrt(edu_chd_est_mv$StdError[1]^2/edu_chd_est_clump$Estimate^2 +
                edu_chd_est_mv$Estimate[1]^2 * edu_chd_est_clump$StdError^2 / edu_chd_est_clump$Estimate^4 -
                2 * edu_chd_est_mv$Estimate[1] * edu_chd_est_mv$StdError[1] * edu_chd_est_clump$StdError * rho1 / edu_chd_est_clump$Estimate^3)
c(Est_ivw - qnorm(0.975) * SE_ivw, Est_ivw + qnorm(0.975) * SE_ivw)

Est_mle = 1 - edu_chd_est_mv_me$thest[1] / edu_chd_est_me_clump$thest
rho2 = 0
SE_mle = sqrt(edu_chd_est_mv_me$Var[1, 1]/edu_chd_est_me_clump$thest^2 +
                edu_chd_est_mv_me$thest[1]^2 * edu_chd_est_me_clump$Var[1, 1] / edu_chd_est_me_clump$thest^4 -
                2 * edu_chd_est_mv_me$thest[1] * sqrt(edu_chd_est_mv_me$Var[1, 1]) * sqrt(edu_chd_est_me_clump$Var[1, 1]) * rho2 / edu_chd_est_me_clump$thest^3)
c(Est_mle - qnorm(0.975) * SE_mle, Est_mle + qnorm(0.975) * SE_mle)
