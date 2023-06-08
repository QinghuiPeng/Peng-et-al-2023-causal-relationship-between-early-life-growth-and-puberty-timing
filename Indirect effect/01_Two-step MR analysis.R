library(reshape)
library(TwoSampleMR)
## Two-step：step one ##
#################
## main method ##
# exposure
exp_path <- 'fetalIV1_clumped.tsv'
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
# outcome
CBMI_out <-read_outcome_data(filename = "CBMI2019.txt",   #outcome data#
                                 snps = BW_exp$SNP,
                                 sep = "\t",
                                 phenotype_col = "Phenotype",     
                                 snp_col = "SNP",
                                 beta_col = "beta",
                                 se_col = "se",
                                 eaf_col = "eaf",
                                 effect_allele_col = "A1",
                                 other_allele_col = "A2",
                                 pval_col = "pval"
)
set.seed(2022052005)
mydata5 <- harmonise_data(exposure_dat=BW_exp,
                         outcome_dat=CBMI_out,
                         action= 1
)
mr_results5<- mr(mydata5, 
                parameters =default_parameters(),
                method_list=c('mr_ivw_mre',
                              'mr_egger_regression',
                              'mr_weighted_median'))
mr_results5

## sensitivity method ##
#################
# MR-PRESSO
library(MRPRESSO)
mr_pre5 <- mr_presso(BetaOutcome ="beta.outcome", 
                     BetaExposure = "beta.exposure", 
                     SdOutcome ="se.outcome", 
                     SdExposure = "se.exposure", 
                     OUTLIERtest = TRUE,
                     DISTORTIONtest = TRUE, 
                     data = mydata5, 
                     NbDistribution = 1000, 
                     SignifThreshold = 0.1)
mr_pre5

## Sensitivity test ##
#################
# horizontal pleiotropy
het5 <- mr_heterogeneity(mydata5) 
het5
pleio5 <- mr_pleiotropy_test(mydata5)
pleio5

# change label
mydata5$exposure[mydata5$exposure == "exposure"] <- "birth weight"
mydata5$outcome[mydata5$outcome == "outcome"] <- "childhood BMI" 
mr_results5$exposure[mr_results5$exposure == "exposure"] <- "birth weight"
mr_results5$outcome[mr_results5$outcome == "outcome"] <- "childhood BMI" 

# lOO plots
res_loo <- mr_leaveoneout(mydata5) 
LOO5 <- mr_leaveoneout_plot(res_loo)
LOO5[[1]]

# Forest plots
res_single <- mr_singlesnp(mydata5,all_method = 
                             c("mr_ivw", 
                               "mr_egger_regression",
                               "mr_weighted_median"))
FOR5 <- mr_forest_plot(res_single)
FOR5[[1]]

#scatter plots
SCA5 <- mr_scatter_plot(mr_results5,mydata5)
SCA5[[1]]

#Funnel plots
res_single <- mr_singlesnp(mydata5,all_method = 
                             c("mr_ivw", 
                               "mr_egger_regression",
                               "mr_weighted_median"))
FUN5 <- mr_funnel_plot(res_single)
FUN5[[1]]

## Two-step：step two ## Same estimate of causal effect as CBMI for PT (liberal analysis)
#################

## indirect effect ##
#################
# Pval（multivariate delta method）#
data1 <- read.csv("indirect.csv", 
                  header=T)
data1$se <- sqrt(data1$b1^2*data1$se2^2+data1$b2^2*data1$se1^2)
data1$beta <- data1$b1*data1$b2
data1$L <- data1$beta-1.96*data1$se
data1$U <- data1$beta+1.96*data1$se
data1$Pval <- 2*pnorm(abs(data1$beta/data1$se),lower.tail = FALSE)
data1