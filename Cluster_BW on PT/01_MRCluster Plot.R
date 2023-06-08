library(ggplot2)
library(ISLR)
## MR-Cluster analysis ##
#################
library(mrclust)
data1 <-read.table("MRcluster.txt",header=T)
rsid = data1$SNP
bx = data1$bx
bxse = data1$bxse
by = data1$by
byse = data1$byse
ratio_est = by/bx
ratio_est_se = byse/abs(bx)

# Running an MR-Clust analysis 
set.seed(2022052007)
res_em = mr_clust_em(theta = ratio_est, theta_se = ratio_est_se, bx = bx,
                     by = by, bxse = bxse, byse = byse, obs_names = rsid)
names(res_em)

# Allocation probabilities and the best clustering of the data
names(res_em$results)
head(res_em$results$all,n=24)  # allocated to "all" clusters 
head(res_em$results$best,n=24)  # allocated to the "best" cluster
#Identified clusters
relab = c(2,1);
clsts = res_em$plots$two_stage$data$cluster_class;
clsts[clsts=="Junk"]<-"Null"
clst_class = rep("Null", length(clsts));
clst_class_num = rep(3, length(clsts));
for(i in unique(clsts)){
  if(i!="Null"){
    clst_class[clsts==as.character(i)] = as.character(relab[as.numeric(i)]);
    clst_class_num[clsts==as.character(i)]= relab[as.numeric(i)];
  }
}
res_em$results$best$cluster_class = clst_class;
res_em$results$best$cluster = clst_class_num;
res_em$plots$two_stage = two_stage_plot(res_em$results$best, bx, by, bxse, byse, rsid)
library(dplyr)
res_em$plots$two_stage$data = res_em$plots$two_stage$data %>% 
  filter(probability >=0.7) %>%
  # select(-cluster_class) %>%
  # left_join(res_em1_junk$plots$two_stage$data %>% select(observation, cluster_class)) %>% 
  library(tidyverse)
res_em$plots$two_stage$data <- dplyr::mutate(res_em$plots$two_stage$data, text_label =
                  ifelse(observation == "rs35261542","CDKAL1",
                         ifelse(observation == "rs6575803","MIR2392",
                                ifelse(observation == "rs2551347","KLHL29",
                                       ifelse(observation == "rs6930558","NMBR",
                                              ifelse(observation == "rs7819593","ZFPM2",""
                                                     ))))))
head(res_em$plots$two_stage$data$text_label,n=64)

# Running plot
plot_out = res_em$plots$two_stage + theme_bw() +
  #ggplot(aes(color = ~cluster_class,  size = ~probability)) + 
  ggplot2::xlab("Genetic association with birth weight") + 
  ggplot2::ylab("Genetic association with puberty timing") + ggplot2::ggtitle("") + 
  ggplot2::theme(axis.text=element_text(size=14),axis.title=element_text(size=18),
                 legend.text=element_text(size=14), legend.title=element_text(size=18),legend.key.size = unit(2, 'lines'))
plot_out = plot_out + scale_color_manual(values = c('#666666','#3399CC',"#FFCC66")) +
  ggplot2::geom_point(aes(text = res_em$plots$two_stage$data$observation, color = res_em$plots$two_stage$data$cluster_class))+
  #ggplot2::geom_line(aes(color = res_em$plots$two_stage$data$cluster_class, linetype = 2)) +
  ggrepel::geom_text_repel(mapping = aes(label = res_em$plots$two_stage$data$text_label,linetype=res_em$plots$two_stage$data$cluster_class), size = 6,nudge_y = -0.003, nudge_x = -0.003, min.segment.length = 1)

plot_out
