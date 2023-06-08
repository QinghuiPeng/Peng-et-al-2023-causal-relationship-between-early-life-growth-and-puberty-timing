library(devtools)
library(TwoSampleMR)
## Liberal analysis ##
# outcome --> puberty timing, GWAS download from xxx #
# exposure --> birth weight, SNPS were extract  from  XXX #
#################
# format exposure
BW_exp1 <- read_exposure_data(filename = 'fetalIV1-exp.txt',   #exposure data#
                              clump = FALSE,
                              sep= "\t",
                              phenotype_col = "birth weight",
                              snp_col = "SNP",
                              beta_col = "beta",
                              se_col = "se",
                              effect_allele_col ="A1",
                              other_allele_col = "A2",
                              eaf_col = "eaf",
                              pval_col = "pval",
                              samplesize_col = "n",
                              id_col = "birth weight"
                              )
# exclude rare variants (eaf＜0.01 & 1-eaf＞0.99, rs138715366)
# exclude variants located in imprinted genes (rs11042596), ref from https://www.geneimprint.com/site/genes-by-species
# exclude variants had maternal genetic effect (MTA：rs10872678，MNTA：rs560887), ref from https://www.nature.com/articles/s41588-021-00896-x
ex1 <- which(BW_exp1$SNP=="rs138715366")
a <- BW_exp1[-ex1,]
ex2 <- which(a$SNP=="rs11042596")
b <- a[-ex2,]
ex3 <- which(b$SNP=="rs10872678")
c <- b[-ex3,]
ex4 <- which(c$SNP=="rs560887")
BW_exp2 <- c[-ex4,]
# clumping
BW_exp <- clump_data(BW_exp2
                     ,clump_r2=0.001
                     ,clump_kb=10000)
# save results
write.table(BW_exp, file='fetalIV1_clumped.tsv', sep='\t', quote=F, row.names=F)


## conservative analysis ##
# outcome --> puberty timing, GWAS download from xxx #
# exposure --> birth weight, SNPS were extract  from  XXX #
#################
# format exposure
BW_exp1 <- read_exposure_data(filename = 'fetalIV1-exp.txt',   #exposure data#
                              clump = FALSE,
                              sep= "\t",
                              phenotype_col = "birth weight",
                              snp_col = "SNP",
                              beta_col = "beta",
                              se_col = "se",
                              effect_allele_col ="A1",
                              other_allele_col = "A2",
                              eaf_col = "eaf",
                              pval_col = "pval",
                              samplesize_col = "n",
                              id_col = "birth weight"
)
# else exclude variants related with factors lying in potential pleiotropy pathways from birth weight to puberty timing (rs1012167,rs11708067), ref from http://www.phenoscanner.medschl.cam.ac.uk/
# else exclude outliers identified by MR-PRESSO (rs1482852,rs2551347,rs35261542)
ex1 <- which(BW_exp1$SNP=="rs138715366")
a <- BW_exp1[-ex1,]
ex2 <- which(a$SNP=="rs11042596")
b <- a[-ex2,]
ex3 <- which(b$SNP=="rs10872678")
c <- b[-ex3,]
ex4 <- which(c$SNP=="rs560887")
d <- c[-ex4,]
ex5 <- which(d$SNP=="rs1012167")
e <- d[-ex5,]
ex6 <- which(e$SNP=="rs11708067")
f <- e[-ex6,]
ex7 <- which(f$SNP=="rs1482852")
g <- f[-ex7,]
ex8 <- which(g$SNP=="rs2551347")
h <- g[-ex8,]
ex9 <- which(h$SNP=="rs35261542")
BW_exp2 <- h[-ex9,]
# clumping
BW_exp <- clump_data(BW_exp2
                     ,clump_r2=0.001
                     ,clump_kb=10000)
# save results
write.table(BW_exp, file='fetalIV2_clumped.tsv', sep='\t', quote=F, row.names=F)