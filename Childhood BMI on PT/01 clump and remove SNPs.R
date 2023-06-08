library(devtools)
library(TwoSampleMR)
## Liberal analysis ##
# outcome --> puberty timing, GWAS download from xxx #
# exposure --> birth weight, SNPS were extract  from  XXX #
#################
# format exposure
CBMI_exp1 <- read_exposure_data(filename = 'CBMIIV1-exp.txt',   #exposure data#
                              clump = FALSE,
                              sep= "\t",
                              phenotype_col = "Childhood BMI",
                              snp_col = "SNP",
                              beta_col = "beta",
                              se_col = "se",
                              effect_allele_col ="A1",
                              other_allele_col = "A2",
                              eaf_col = "eaf",
                              pval_col = "pval"
                              )
# clumping
CBMI_exp <- clump_data(CBMI_exp1
                     ,clump_r2=0.001
                     ,clump_kb=10000)
# save results
write.table(CBMI_exp, file='CBMIIV1_clumped.tsv', sep='\t', quote=F, row.names=F)


## conservative analysis ##
# outcome --> puberty timing, GWAS download from xxx #
# exposure --> birth weight, SNPS were extract  from  XXX #
#################
# format exposure
CBMI_exp1 <- read_exposure_data(filename = 'CBMIIV1-exp.txt',   #exposure data#
                              clump = FALSE,
                              sep= "\t",
                              phenotype_col = "Childhood BMI",
                              snp_col = "SNP",
                              beta_col = "beta",
                              se_col = "se",
                              effect_allele_col ="A1",
                              other_allele_col = "A2",
                              eaf_col = "eaf",
                              pval_col = "pval"
)
# else exclude outliers identified by MR-PRESSO (rs10493544,rs539515,rs6567160)
ex1 <- which(CBMI_exp1$SNP=="rs10493544")
a <- CBMI_exp1[-ex1,]
ex2 <- which(a$SNP=="rs539515")
b <- a[-ex2,]
ex3 <- which(b$SNP=="rs6567160")
CBMI_exp2 <- b[-ex3,]
# clumping
CBMI_exp <- clump_data(CBMI_exp2
                     ,clump_r2=0.001
                     ,clump_kb=10000)
# save results
write.table(CBMI_exp, file='CBMIIV2_clumped.tsv', sep='\t', quote=F, row.names=F)