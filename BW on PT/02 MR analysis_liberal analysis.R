library(TwoSampleMR)
## main method ##
#################
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
# PT GWAS, Donwloaded from XXX
# add eaf col
library(reshape)
dat1 <-read.table("PT2019.txt",header=T)  
dat2 <-read.table("fetalIV1-exp.txt",header=T)
dat3 <- merge(dat2,dat1,by="SNP")
dat3$A1.y <- toupper(dat3$A1.y)
dat3$A2.y <- toupper(dat3$A2.y)
dat3$Minor_Allele <- toupper(dat3$Minor_Allele)
dat3$m <- 1-dat3$eaf
dat3$maf <- ifelse(dat3$eaf < 0.5,dat3$eaf,dat3$m)
dat3$n <- 1-dat3$maf
dat3$eaf.y <- ifelse(dat3$A1.y==dat3$Minor_Allele,dat3$maf,dat3$n)
dat3 <- subset(dat3,select=-c(A1.x,A2.x,beta.x,pval.x,se.x,eaf,m)) 
dat4 <- rename(dat3,c(A1.y = "A1",A2.y = "A2",beta.y="beta",se.y="se",pval.y="pval",eaf.y="eaf"))
# save results
write.table(dat4,file="fetalIV1-out.txt",sep = "\t",quote=FALSE,row.names = FALSE)

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
set.seed(2022052001)

# harmonize
mydata1 <- harmonise_data(exposure_dat=BW_exp,
                          outcome_dat=menarche_out,
                          action= 2
)

# mr 
mr_results1 <- mr(mydata1, 
                  parameters =default_parameters(),
                  method_list=c('mr_ivw_mre',
                                'mr_egger_regression',
                                'mr_weighted_median'))
mr_results1

## sensitivity method ##
#################
# MR-PRESSO
library(MRPRESSO)
mr_pre1 <- mr_presso(BetaOutcome ="beta.outcome", 
                     BetaExposure = "beta.exposure", 
                     SdOutcome ="se.outcome", 
                     SdExposure = "se.exposure", 
                     OUTLIERtest = TRUE,
                     DISTORTIONtest = TRUE, 
                     data = mydata1, 
                     NbDistribution = 1000, 
                     SignifThreshold = 0.1,
                     seed=202205201)
mr_pre1

## Sensitivity test ##
#################
# horizontal pleiotropy
het1 <- mr_heterogeneity(mydata1) 
het1
pleio1 <- mr_pleiotropy_test(mydata1)
pleio1

# change label
mydata1$exposure[mydata1$exposure == "exposure"] <- "birth weight"
mydata1$outcome[mydata1$outcome == "outcome"] <- "puberty timing" 
mr_results1$exposure[mr_results1$exposure == "exposure"] <- "birth weight"
mr_results1$outcome[mr_results1$outcome == "outcome"] <- "puberty timing" 

# lOO plots
res_loo <- mr_leaveoneout(mydata1) 
LOO1 <- mr_leaveoneout_plot(res_loo)
LOO1[[1]]

# Forest plots
res_single <- mr_singlesnp(mydata1,all_method = 
                             c("mr_ivw", 
                               "mr_egger_regression",
                               "mr_weighted_median"))
FOR1 <- mr_forest_plot(res_single)
FOR1[[1]]

#scatter plots
SCA1 <- mr_scatter_plot(mr_results1,mydata1)
SCA1[[1]]

#Funnel plots
res_single <- mr_singlesnp(mydata1,all_method = 
                             c("mr_ivw", 
                               "mr_egger_regression",
                               "mr_weighted_median"))
FUN1 <- mr_funnel_plot(res_single)
FUN1[[1]]
