library(TwoSampleMR)
## main method ##
#################
# exposure
exp_path <- 'fetalIV2_clumped.tsv'
BW_exp <- read_exposure_data(filename = exp_path, 
                             clump = FALSE,
                             sep= "\t",
                             phenotype_col = "birth weight",
                             snp_col = "SNP",
                             beta_col = "beta.exposure",
                             se_col = "se.exposure",
                             effect_allele_col ="effect_allele.exposure",
                             other_allele_col = "other_allele.exposure",
                             eaf_col = "eaf.exposure",
                             pval_col = "pval.exposure"
)

# import outcome GWAS
# PT GWAS, Donwloaded from XXX
# outcome
menarche_out <-read_outcome_data(filename = "fetalIV1-out.txt", 
                                 snps = BW_exp$SNP,
                                 sep = "\t",
                                 phenotype_col = "puberty timing",     
                                 snp_col = "SNP",
                                 beta_col = "beta",
                                 se_col = "se",
                                 effect_allele_col = "A1",
                                 other_allele_col = "A2",
                                 eaf_col = "eaf",
                                 pval_col = "pval"
)
set.seed(2022052002)

# harmonize
mydata2 <- harmonise_data(exposure_dat=BW_exp,
                          outcome_dat=menarche_out,
                          action= 2
)

# mr 
mr_results2 <- mr(mydata2, 
                  parameters =default_parameters(),
                  method_list=c('mr_ivw_mre',
                                'mr_egger_regression',
                                'mr_weighted_median'))
mr_results2

## sensitivity method ##
#################
# MR-PRESSO
library(MRPRESSO)
mr_pre2 <- mr_presso(BetaOutcome ="beta.outcome", 
                     BetaExposure = "beta.exposure", 
                     SdOutcome ="se.outcome", 
                     SdExposure = "se.exposure", 
                     OUTLIERtest = TRUE,
                     DISTORTIONtest = TRUE, 
                     data = mydata2, 
                     NbDistribution = 1000, 
                     SignifThreshold = 0.1,
                     seed=202205202)
mr_pre2

## Sensitivity test ##
#################
# horizontal pleiotropy
het2 <- mr_heterogeneity(mydata2) 
het2
pleio2 <- mr_pleiotropy_test(mydata2)
pleio2

# change label
mydata2$exposure[mydata2$exposure == "exposure"] <- "birth weight"
mydata2$outcome[mydata2$outcome == "outcome"] <- "puberty timing" 
mr_results2$exposure[mr_results2$exposure == "exposure"] <- "birth weight"
mr_results2$outcome[mr_results2$outcome == "outcome"] <- "puberty timing" 

# lOO plots
res_loo <- mr_leaveoneout(mydata2) 
LOO2 <- mr_leaveoneout_plot(res_loo)
LOO2[[1]]

# Forest plots
res_single <- mr_singlesnp(mydata2,all_method = 
                             c("mr_ivw", 
                               "mr_egger_regression",
                               "mr_weighted_median"))
FOR2 <- mr_forest_plot(res_single)
FOR2[[1]]

#scatter plots
SCA2 <- mr_scatter_plot(mr_results2,mydata2)
SCA2[[1]]

#Funnel plots
res_single <- mr_singlesnp(mydata2,all_method = 
                             c("mr_ivw", 
                               "mr_egger_regression",
                               "mr_weighted_median"))
FUN2 <- mr_funnel_plot(res_single)
FUN2[[1]]
