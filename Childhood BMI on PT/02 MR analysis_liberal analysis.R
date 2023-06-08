library(TwoSampleMR)
## main method ##
#################
# exposure
exp_path <- 'CBMIIV1_clumped.tsv'
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
# add eaf col
library(reshape)
dat1 <-read.table("PT2019.txt",header=T)  
dat2 <-read.table("CBMIIV1-exp.txt",header=T)
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
write.table(dat4,file="CBMIIV1-out.txt",sep = "\t",quote=FALSE,row.names = FALSE)

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
set.seed(2022052003)

# harmonize
mydata3 <- harmonise_data(exposure_dat=CBMI_exp,
                          outcome_dat=menarche_out,
                          action= 1
)

# mr 
mr_results3 <- mr(mydata3, 
                  parameters =default_parameters(),
                  method_list=c('mr_ivw_mre',
                                'mr_egger_regression',
                                'mr_weighted_median'))
mr_results3

## sensitivity method ##
#################
# MR-PRESSO
library(MRPRESSO)
mr_pre3 <- mr_presso(BetaOutcome ="beta.outcome", 
                     BetaExposure = "beta.exposure", 
                     SdOutcome ="se.outcome", 
                     SdExposure = "se.exposure", 
                     OUTLIERtest = TRUE,
                     DISTORTIONtest = TRUE, 
                     data = mydata3, 
                     NbDistribution = 1000, 
                     SignifThreshold = 0.1,
                     seed=202205203)
mr_pre3

## Sensitivity test ##
#################
# horizontal pleiotropy
het3 <- mr_heterogeneity(mydata3) 
het3
pleio3 <- mr_pleiotropy_test(mydata3)
pleio3

# change label
mydata3$exposure[mydata3$exposure == "exposure"] <- "childhood BMI"
mydata3$outcome[mydata3$outcome == "outcome"] <- "puberty timing" 
mr_results3$exposure[mr_results3$exposure == "exposure"] <- "childhood BMI"
mr_results3$outcome[mr_results3$outcome == "outcome"] <- "puberty timing" 

# lOO plots
res_loo <- mr_leaveoneout(mydata3) 
LOO3 <- mr_leaveoneout_plot(res_loo)
LOO3[[1]]

# Forest plots
res_single <- mr_singlesnp(mydata3,all_method = 
                             c("mr_ivw", 
                               "mr_egger_regression",
                               "mr_weighted_median"))
FOR3 <- mr_forest_plot(res_single)
FOR3[[1]]

#scatter plots
SCA3 <- mr_scatter_plot(mr_results3,mydata3)
SCA3[[1]]

#Funnel plots
res_single <- mr_singlesnp(mydata3,all_method = 
                             c("mr_ivw", 
                               "mr_egger_regression",
                               "mr_weighted_median"))
FUN3 <- mr_funnel_plot(res_single)
FUN3[[1]]
