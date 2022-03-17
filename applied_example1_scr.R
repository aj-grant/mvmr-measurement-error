################################################################################
#Script used for the first applied example
################################################################################

library(TwoSampleMR)
library(tidyverse)
library(MendelianRandomization)
library(RColorBrewer)
library(ggstance)

bmi_all = read.table('bmi.giant-ukbb.meta-analysis.combined.23May2018.txt', header = TRUE)
bmi_all$SNP = sapply(1:length(bmi_all$SNP), function(j){substr(bmi_all$SNP[j], 1, str_locate(bmi_all$SNP[j], ":")-1)})
bmi_sig = bmi_all[bmi_all$P < 5e-8, ]

whr_all = read.table('whr.giant-ukbb.meta-analysis.combined.23May2018.txt', header = TRUE)
whr_all$SNP = sapply(1:length(whr_all$SNP), function(j){substr(whr_all$SNP[j], 1, str_locate(whr_all$SNP[j], ":")-1)})
whr_sig = whr_all[whr_all$P < 5e-8, ]

comb_snps = data.frame("SNP" = unique(c(bmi_sig$SNP, whr_sig$SNP)))

bmi_comb = inner_join(comb_snps, bmi_all, by = "SNP")
bmi_comb_fm = format_data(bmi_comb, snp_col = "SNP", beta_col = "BETA", se_col = "SE", effect_allele_col = "Tested_Allele",
                          other_allele_col = "Other_Allele", pval_col = "P", "chr_col" = "CHR", pos_col = "POS")
comb_clump = clump_data(bmi_comb_fm, clump_r2 = 0.001)

bmi_snps_fm = format_data(bmi_sig, snp_col = "SNP", beta_col = "BETA", se_col = "SE", effect_allele_col = "Tested_Allele",
                          other_allele_col = "Other_Allele", pval_col = "P", "chr_col" = "CHR", pos_col = "POS")
bmi_clump = clump_data(bmi_snps_fm, clump_r2 = 0.001)

whr_snps_fm = format_data(whr_sig, snp_col = "SNP", beta_col = "BETA", se_col = "SE", effect_allele_col = "Tested_Allele",
                          other_allele_col = "Other_Allele", pval_col = "P", "chr_col" = "CHR", pos_col = "POS")
whr_clump = clump_data(whr_snps_fm, clump_r2 = 0.001)

bmi_dat = bmi_clump[, c("SNP", "effect_allele.exposure", "beta.exposure", "se.exposure")]
names(bmi_dat) = c("SNP", "EA", "beta_bmi", "se_bmi")

whr_dat = whr_clump[, c("SNP", "effect_allele.exposure", "beta.exposure", "se.exposure")]
names(whr_dat) = c("SNP", "EA", "beta_whr", "se_whr")

all_dat = comb_clump[, c("SNP", "effect_allele.exposure", "beta.exposure", "se.exposure")]
names(all_dat) = c("SNP", "EA", "beta_bmi", "se_bmi")
all_dat = inner_join(all_dat, whr_all[, c("SNP", "BETA", "SE")], by = "SNP")
names(all_dat) = c("SNP", "EA", "beta_bmi", "se_bmi", "beta_whr", "se_whr")

################################################################################
#Combined SNP list
snps_all = full_join(data.frame("SNP" = all_dat$SNP), data.frame("SNP" = all_dat_m$SNP), by = "SNP")
snps_all = full_join(snps_all, data.frame("SNP" = all_dat_f$SNP), by = "SNP")

pheno_reps = round(length(snps_all$SNP) / 100)
A = lapply(1:pheno_reps, function(j){
  ph = phenoscanner(snpquery = snps_all[(100*(j-1)+1):(100*j), ], pvalue = 1)
  ph$results
})
A[[(pheno_reps+1)]] = phenoscanner(snpquery = snps_all[(pheno_reps * 100 + 1):(pheno_reps * 100 + dim(snps_all)[1] %% 100), ],
                                   pvalue = 1)$results
pheno_all = rbind(A[[1]], A[[2]])
if (pheno_reps > 3){
  for (j in 3:(pheno_reps+1)){
    pheno_all = rbind(pheno_all, A[[j]])
  }
} else {
  pheno_all = rbind(bmi_pheno_all, A[[3]])
  pheno_all = rbind(bmi_pheno_all, A[[4]])
}
rm(A)

pheno_reps = (length(bmi_dat$SNP) - length(bmi_dat$SNP) %% 100)/100
A = lapply(1:pheno_reps, function(j){
  ph = phenoscanner(snpquery = bmi_dat$SNP[(100*(j-1)+1):(100*j)], pvalue = 1)
  ph$results
})
A[[(pheno_reps+1)]] = phenoscanner(snpquery = bmi_dat$SNP[(pheno_reps * 100 + 1):(pheno_reps * 100 + dim(bmi_dat)[1] %% 100)],
                                   pvalue = 1)$results
pheno_bmi = rbind(A[[1]], A[[2]])
if (pheno_reps > 3){
  for (j in 3:(pheno_reps+1)){
    pheno_bmi = rbind(pheno_bmi, A[[j]])
  }
} else {
  pheno_bmi= rbind(pheno_bmi, A[[3]])
  pheno_bmi = rbind(pheno_bmi, A[[4]])
}
rm(A)

pheno_reps = (length(whr_dat$SNP) - length(whr_dat$SNP) %% 100)/100
A = lapply(1:pheno_reps, function(j){
  ph = phenoscanner(snpquery = whr_dat$SNP[(100*(j-1)+1):(100*j)], pvalue = 1)
  ph$results
})
A[[(pheno_reps+1)]] = phenoscanner(snpquery = whr_dat$SNP[(pheno_reps * 100 + 1):(pheno_reps * 100 + dim(whr_dat)[1] %% 100)],
                                   pvalue = 1)$results
pheno_whr = rbind(A[[1]], A[[2]])
if (pheno_reps > 3){
  for (j in 3:(pheno_reps+1)){
    pheno_whr = rbind(pheno_whr, A[[j]])
  }
} else {
  pheno_whr= rbind(pheno_whr, A[[3]])
  pheno_whr = rbind(pheno_whr, A[[4]])
}

chd_dat = pheno_all %>% filter(dataset == "CARDIoGRAMplusC4D_CHD_Mixed_2015") %>% select(c("rsid", "a1", "beta", "se"))
chd_dat$beta = as.numeric(chd_dat$beta)
chd_dat$se = as.numeric(chd_dat$se)
names(chd_dat) = c("SNP", "EA_chd", "beta_chd", "se_chd")

chd_bmi_dat = pheno_bmi %>% filter(dataset == "CARDIoGRAMplusC4D_CHD_Mixed_2015") %>% select(c("rsid", "a1", "beta", "se"))
chd_bmi_dat$beta = as.numeric(chd_bmi_dat$beta)
chd_bmi_dat$se = as.numeric(chd_bmi_dat$se)
names(chd_bmi_dat) = c("SNP", "EA_chd", "beta_chd", "se_chd")

chd_whr_dat = pheno_whr %>% filter(dataset == "CARDIoGRAMplusC4D_CHD_Mixed_2015") %>% select(c("rsid", "a1", "beta", "se"))
chd_whr_dat$beta = as.numeric(chd_whr_dat$beta)
chd_whr_dat$se = as.numeric(chd_whr_dat$se)
names(chd_whr_dat) = c("SNP", "EA_chd", "beta_chd", "se_chd")

################################################################################
#Harmonise
bmi_dat = inner_join(bmi_dat, chd_bmi_dat, by = "SNP")
bmi_dat$beta_chd = bmi_dat$beta_chd * (as.numeric(as.character(bmi_dat$EA) == as.character(bmi_dat$EA_chd))*2-1)
bmi_dat = bmi_dat[, -which(names(bmi_dat)=="EA_chd")]

whr_dat = inner_join(whr_dat, chd_whr_dat, by = "SNP")
whr_dat$beta_chd = whr_dat$beta_chd * (as.numeric(as.character(whr_dat$EA) == as.character(whr_dat$EA_chd))*2-1)
whr_dat = whr_dat[, -which(names(whr_dat)=="EA_chd")]

all_dat = inner_join(all_dat, chd_dat, by = "SNP")
all_dat$beta_chd = all_dat$beta_chd * (as.numeric(as.character(all_dat$EA) == as.character(all_dat$EA_chd))*2-1)
all_dat = all_dat[, -7]

################################################################################
#MR Analyses
mrob_bmi = mr_input(snps = bmi_dat$SNP, bx = bmi_dat$beta_bmi, bxse = bmi_dat$se_bmi,
                    by = bmi_dat$beta_chd, byse = bmi_dat$se_chd)
est_bmi_ivw = mr_ivw(mrob_bmi)

mrob_whr = mr_input(snps = whr_dat$SNP, bx = whr_dat$beta_whr, bxse = whr_dat$se_whr,
                    by = whr_dat$beta_chd, byse = whr_dat$se_chd)
est_whr_ivw = mr_ivw(mrob_whr)

set.seed(20201201)
mrob = mr_mvinput(snp = all_dat$SNP, bx = cbind(all_dat$beta_bmi, all_dat$beta_whr),
                  bxse = cbind(all_dat$se_bmi, all_dat$se_whr), by = all_dat$beta_chd, byse = all_dat$se_chd)
est_ivw = mr_mvivw(mrob)
est_me = mrest_me(mrob)
est_me_cor = mrest_me_cor(mrob, corX = matrix(c(1, 0.433, 0.433, 1), nrow = 2))
round(cbind(est_me$thest, sqrt(c(est_me$Var[1, 1], est_me$Var[2, 2])),
            est_me$thest - qnorm(0.975, 0, 1) * sqrt(c(est_me$Var[1, 1], est_me$Var[2, 2])),
            est_me$thest + qnorm(0.975, 0, 1) * sqrt(c(est_me$Var[1, 1], est_me$Var[2, 2]))), 3)
round(cbind(est_me_cor$thest, sqrt(c(est_me_cor$Var[1, 1], est_me_cor$Var[2, 2])),
            est_me_cor$thest - qnorm(0.975, 0, 1) * sqrt(c(est_me_cor$Var[1, 1], est_me_cor$Var[2, 2])),
            est_me_cor$thest + qnorm(0.975, 0, 1) * sqrt(c(est_me_cor$Var[1, 1], est_me_cor$Var[2, 2]))), 3)

plot_dat = data.frame(
  "Estimate" = c(est_bmi_ivw$Estimate, est_whr_ivw$Estimate, est_ivw$Estimate, est_me$thest, est_me_cor$thest),
  "CILower" = c(est_bmi_ivw$CILower, est_whr_ivw$CILower,
                est_ivw$CILower, est_me$thest - qnorm(0.975, 0, 1) * sqrt(diag(est_me$Var)),
                est_me_cor$thest - qnorm(0.975, 0, 1) * sqrt(diag(est_me_cor$Var))),
  "CIUpper" = c(est_bmi_ivw$CIUpper, est_whr_ivw$CIUpper, est_ivw$CIUpper,
                est_me$thest + qnorm(0.975, 0, 1) * sqrt(diag(est_me$Var)),
                est_me_cor$thest + qnorm(0.975, 0, 1) * sqrt(diag(est_me_cor$Var))),
  "Exposure" = rep(c("BMI", "WHR"), 4),
  "Method" = c(rep("Univ.", 2), rep("IVW", 2), rep("MLE", 2), rep("MLE (cor)", 2)))

plot_dat$Method = factor(plot_dat$Method, levels = c("Univ.", "IVW", "MLE", "MLE (cor)"))

ggplot(plot_dat2, aes(y = Exposure, x = Estimate)) +
  geom_pointrangeh(aes(xmin = CILower, xmax = CIUpper, color = Method, shape = Method),
                  position = position_dodge2(width = 0.5, reverse = T)) +
  geom_vline(xintercept = 0, lty = 2) + theme_classic() +
  scale_color_manual(values = brewer.pal(4, "Dark2")[c(4, 1, 2, 3)]) +
  scale_shape_manual(values = c(18, 16, 17, 15)) +
  theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.text = element_text(size = 8),
        legend.title = element_text(size = 8)) + xlab("Estimate (95% CI)")
