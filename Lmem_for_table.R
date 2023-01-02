### lmem code to analysis labels, which put in one common table 

library(reshape2)
library(data.table)
library(ggplot2)
library(lme4)
library(emmeans)
library(lmerTest)
library(stringi)
library(stringr)
library(dplyr)
library(purrr)
library(tidyverse)
library(scales)
library(optimx)
library(remotes)
library(permutes)
library(buildmer)
library(Matrix)
options(scipen = 999)

data_path<-'/Users/kristina/Documents/stc/df_900_300_mean'


read_plus <- function(flnm) {
  read_csv(flnm) %>% 
    mutate(filename = flnm)
}
######### Looad dataframes #########
tbl_with_sources <-list.files(data_path,pattern = "*.csv", full.names = T) %>% 
  map_df(~fread(.))

setDT(tbl_with_sources)

label<- unique(tbl_with_sources$label)

########## Remove outliers (more than 2 sigma) ########

data_beta_sum = tbl_with_sources %>%
  group_by(subject, label) %>%
  summarise(beta_mean = mean(beta_power),
            beta_sd = sd(beta_power)) %>%
  ungroup() %>%
  mutate(beta_high = beta_mean + (2 * beta_sd)) %>%
  mutate(beta_low = beta_mean - (2 * beta_sd))

data_accuracy_clean = tbl_with_sources %>%
  inner_join(data_beta_sum) %>%
  filter(beta_power < beta_high) %>%
  filter(beta_power > beta_low)


###### create df for stc plots with mean beta power ########
mean_df_step_1<- tbl_with_sources%>% group_by(subject, trial_type,label,round) %>% 
  dplyr::summarise(mean_beta = mean(beta_power))
mean_df_step_2<- mean_df_step_1%>% group_by(subject, trial_type,label) %>% 
  dplyr::summarise(mean_beta = mean(mean_beta))
mean_df_step_3<- mean_df_step_2%>% group_by(trial_type,label) %>% 
  dplyr::summarise(mean_beta = mean(mean_beta))
setwd("/Users/kristina/Documents/stc/lmem_label")
######## create data with differences in beta power choose needed condition #######
hp<- filter(mean_df_step_3, trial_type=="norisk")
lp<- filter(mean_df_step_3, trial_type=="risk")
diff<- as.data.frame(lp$mean_beta-hp$mean_beta)
colnames(diff)<- "mean_beta"
write.csv(diff, "dfff_lp_minus_label_900_300.csv")





############ LMEM with post-hocs ############

p_vals1 <- data.table()
for (l in 1:length(labels)){
  temp<- subset(data_accuracy_clean,label == labels[l])
  
  print(temp)
  m <-  lmer(beta_power ~ trial_type + (1|subject), data = temp)
  summary(m)
  
  #Tuk<-data.table(summary(emmeans(regrid(ref_grid(m, transform = "response")), pairwise ~ trial_type, adjust = 'tukey',lmer.df = "satterthwaite"))$contrasts)
  Tuk <- data.table(summary(emmeans(m, pairwise ~ trial_type, adjust = 'tukey',lmer.df = "satterthwaite",lmerTest.limit=2498162))$contrasts)
  Tuk[,contrast:=gsub(' - ','_',contrast)]
  Tuk[,p.value:=format(p.value, digits = 3)]
  #Tuk[,contrast:=paste0(trial_type,'-',contrast)]
  #Tuk[,contrast:=paste0(trial_type,'-',feedback_cur,'-',contrast)] ####### if you need triple interaction #########
  columns <- c('contrast','p.value')
  Tuk <- Tuk[,..columns]
  Tuk$interval <- "mean_beta"
  Tuk <- dcast(Tuk,formula = interval~contrast,value.var = 'p.value')
  Tuk$label <- unique(temp$label)
  #Tuk$short_filename <- l
  p_vals1 <- rbind(p_vals1,Tuk)
}
############ LMEM with post-hocs trial_type*feedback ############
p_vals <- data.table()
for (l in 1:length(labels)){
  temp<- subset(data_accuracy_clean,label == labels[l])
  
  m <-  lmer(beta_power~ trial_type*feedback_cur+(1|subject), data = temp,REML=FALSE)
  summary(m)
  
  #Tuk<-data.table(summary(emmeans(regrid(ref_grid(m, transform = "response")), pairwise ~ trial_type, adjust = 'tukey',lmer.df = "satterthwaite"))$contrasts)
  Tuk <- data.table(summary(emmeans(m, pairwise ~ feedback_cur|trial_type, adjust = 'tukey',lmer.df = "satterthwaite",lmerTest.limit=2498162))$contrasts)
  Tuk[,contrast:=gsub(' - ','_',contrast)]
  Tuk[,p.value:=format(p.value, digits = 3)]
  Tuk[,contrast:=paste0(trial_type,'-',contrast)]
  #Tuk[,contrast:=paste0(trial_type,'-',feedback_cur,'-',contrast)] ####### if you need triple interaction #########
  columns <- c('contrast','p.value')
  Tuk <- Tuk[,..columns]
  Tuk$interval <- "mean_beta"
  Tuk <- dcast(Tuk,formula = interval~contrast,value.var = 'p.value')
  Tuk$label <- unique(temp$label)
  #Tuk$short_filename <- l
  p_vals <- rbind(p_vals,Tuk)
}

p_vals[, risk-negative_positive_fdr:=p.adjust(`risk-negative_positive`, method = 'fdr')]

############# LMEM main effects ###########

p_vals <- data.table()
for (l in 1:length(labels)){
  temp<- subset(data_accuracy_clean,label == labels[l])
  
  print(temp)
  m <-  lmer(beta_power ~ trial_type + (1|subject), data = temp,)
  summary(m)
  an <- anova(m)
  #print(an)
  an <- data.table(an,keep.rownames = TRUE)
  an_cols <- c('rn','Pr(>F)') 
  an <- an[, ..an_cols]
  an$`Pr(>F)` <- format(an$`Pr(>F)`, digits = 3)
  an$interval <- "beta_power"
  an$interval <- gsub('beta power','',an$interval)
  an <- dcast(an,formula = interval~rn,value.var = 'Pr(>F)')
  an$label <- unique(temp$label)
  p_vals <- rbind(p_vals,an)
}
setwd("/Users/kristina/Documents/stc/lmem_label")
write.csv(p_vals1, "anova_label_900_300.csv")


