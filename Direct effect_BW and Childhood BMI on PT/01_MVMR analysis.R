library(reshape)
library(MVMR)
## MVMR analysis ##
# outcome --> puberty timing
# exposure --> birth weight and childhood BMI
#################
BWCBMI <-read.table("BWCBMIPT.txt",header=T)
head(BWCBMI)
set.seed(2022052006)
BWCBMI_MVMR <- format_mvmr(BXGs = BWCBMI[,c(4,6)],
                           BYG = BWCBMI[,2],
                           seBXGs = BWCBMI[,c(5,7)],
                           seBYG = BWCBMI[,3],
                           RSID = BWCBMI[,1])
head(BWCBMI_MVMR)
sres <- strength_mvmr(r_input = BWCBMI_MVMR, gencov = 0)
pres <- pleiotropy_mvmr(r_input =BWCBMI_MVMR, gencov = 0)
res <- ivw_mvmr(r_input = BWCBMI_MVMR)