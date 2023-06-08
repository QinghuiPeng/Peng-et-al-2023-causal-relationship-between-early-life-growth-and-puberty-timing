library(TwoSampleMR)
## main method ##
#################
# exposure
exp_path <- 'CBMIIV2_clumped.tsv'
CBMI_exp <- read_exposure_data(filename = exp_path, 
                             clump = FALSE,
                             sep= "\t",
                             phenotype_col = "Childhood BMI",
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
menarche_out <-read_outcome_data(filename = "CBMIIV1-out.txt", 
                                 snps = CBMI_exp$SNP,
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
set.seed(2022052004)

# harmonize
mydata4 <- harmonise_data(exposure_dat=CBMI_exp,
                          outcome_dat=menarche_out,
                          action= 1
)

# mr 
mr_results4 <- mr(mydata4, 
                  parameters =default_parameters(),
                  method_list=c('mr_ivw_mre',
                                'mr_egger_regression',
                                'mr_weighted_median'))
mr_results4

## sensitivity method ##
#################
# MR-PRESSO
library(MRPRESSO)
mr_pre4 <- mr_presso(BetaOutcome ="beta.outcome", 
                     BetaExposure = "beta.exposure", 
                     SdOutcome ="se.outcome", 
                     SdExposure = "se.exposure", 
                     OUTLIERtest = TRUE,
                     DISTORTIONtest = TRUE, 
                     data = mydata4, 
                     NbDistribution = 1000, 
                     SignifThreshold = 0.1,
                     seed=202205204)
mr_pre4

## Sensitivity test ##
#################
# horizontal pleiotropy
het4 <- mr_heterogeneity(mydata4) 
het4
pleio4 <- mr_pleiotropy_test(mydata4)
pleio4

# change label
mydata4$exposure[mydata4$exposure == "exposure"] <- "childhood BMI"
mydata4$outcome[mydata4$outcome == "outcome"] <- "puberty timing" 
mr_results4$exposure[mr_results4$exposure == "exposure"] <- "childhood BMI"
mr_results4$outcome[mr_results4$outcome == "outcome"] <- "puberty timing" 

# lOO plots
res_loo <- mr_leaveoneout(mydata4) 
LOO4 <- mr_leaveoneout_plot(res_loo)
LOO4[[1]]

# Forest plots
res_single <- mr_singlesnp(mydata4,all_method = 
                             c("mr_ivw", 
                               "mr_egger_regression",
                               "mr_weighted_median"))
FOR4 <- mr_forest_plot(res_single)
FOR4[[1]]

#scatter plots
SCA4 <- mr_scatter_plot(mr_results4,mydata4)
SCA4[[1]]

#Funnel plots
res_single <- mr_singlesnp(mydata4,all_method = 
                             c("mr_ivw", 
                               "mr_egger_regression",
                               "mr_weighted_median"))
FUN4 <- mr_funnel_plot(res_single)
FUN4[[1]]
