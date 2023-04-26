### lmem code to analysis labels, which put in one common table 

library(reshape2)
library(data.table)
library(ggplot2)
library(ggpubr)
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

data_path<-'/Users/kristina/Documents/stc/df_1100_1500_khan'


read_plus <- function(flnm) {
  read_csv(flnm) %>% 
    mutate(filename = flnm)
} #read csv files
######### Load dataframes #########
tbl_with_sources <-list.files(data_path,pattern = "*.csv", full.names = T) %>% 
  map_df(~fread(.))
#tbl_with_sources<- na.omit(tbl_with_sources)
setDT(tbl_with_sources)
labels<- unique(tbl_with_sources$label)

autists<- filter(tbl_with_sources, subject %in% c('P301','P304', 'P307','P312','P313', 'P314',
                                                  'P316','P321','P322','P323','P324',
                                                  'P325','P326','P327','P328','P329','P333',
                                                  'P334','P335', 'P341','P342'))
autists$group<- "autists"
normal<- filter(tbl_with_sources, subject %in% c('P001','P004','P019','P021','P022','P034',"P035",
                                                 'P040','P044','P047','P048','P053','P058','P059','P060','P061',
                                                 'P063','P064','P065','P067'))

normal$group<- "normal"

df<- rbind(autists,normal)

lp_hp<- filter(df, trial_type=="norisk"|trial_type=="risk")

######### if you need, remove outliers #######
data_beta_sum = lp_hp %>%
  group_by(subject, label) %>%
  summarise(beta_mean = mean(beta_power),
            beta_sd = sd(beta_power)) %>%
  ungroup() %>%
  mutate(beta_high = beta_mean + (2 * beta_sd)) %>%
  mutate(beta_low = beta_mean - (2 * beta_sd))

data_accuracy_clean = lp_hp %>%
  inner_join(data_beta_sum) %>%
  filter(beta_power < beta_high) %>%
  filter(beta_power > beta_low)




emm_options(lmerTest.limit = 10000)
emm_options(pbkrtest.limit = 9000)


######## LMEM MAIN EFFECTS #########
p_vals <- data.table()
for (l in 1:length(labels)){
  temp<- subset(lp_hp,label == labels[l])
  
  print(temp)
  m <-  lmer(beta_power ~ trial_type*group*feedback_cur + (1|subject), data = temp)
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

setwd("/Users/kristina/Documents/stc/lmem_label/dfs")
write.csv(p_vals, "lmem_label_1100_1500.csv")



###### create table with marginal means based on LMEM moodel #######
marginal <- data.table()

for (l in 1:length(labels)){
  temp<- subset(lp_hp,label == labels[l])
  m <-  lmer(beta_power~ trial_type*group*feedback_cur+(1|subject), data = temp)
  summary(m)
  marginal_em <- emmeans(m, ~ as.factor(feedback_cur|trial_type|group), level = 0.99,lmer.df = "satterthwaite")
  marginal_em<- as.data.frame(marginal_em)
  print(marginal_em)
  
  marginal_em$label <- unique(temp$label)
  #Tuk$emmeans <- 
  marginal <- rbind(marginal,marginal_em)
}


###### trial_type|group ######
hp_aut<- filter(marginal, trial_type=="norisk" &  group=="autists")
lp_aut<-filter(marginal, trial_type=="risk" & group=="autists")
hp_norm<- filter(marginal, trial_type=="norisk" &  group=="normal")
lp_norm<-filter(marginal, trial_type=="risk" & group=="normal")

#    between group #
lp_norm_aut<- lp_aut$emmean - lp_norm$emmean
hp_norm_aut<- hp_aut$emmean - hp_norm$emmean

#     ingroup #
lp_hp_aut<- lp_aut$emmean-hp_aut$emmean
lp_hp_norm<- lp_norm$emmean-hp_norm$emmean

###### trial_type|group|feedback ######
#hp_aut<- filter(marginal, trial_type=="norisk" &  group=="autists")
lp_aut_pos<-filter(marginal, trial_type=="risk" & group=="autists" & feedback_cur=="positive")
lp_norm_pos<-filter(marginal, trial_type=="risk" & group=="normal"& feedback_cur=="positive")

lp_aut_neg<-filter(marginal, trial_type=="risk" & group=="autists" & feedback_cur=="negative")
lp_norm_neg<-filter(marginal, trial_type=="risk" & group=="normal" & feedback_cur=="negative")

#     ingroup #
lp_norm<- lp_norm_neg$emmean - lp_norm_pos$emmean
lp_aut<- lp_aut_neg$emmean - lp_aut_pos$emmean
#    between group #
lp_neg_aut_norm <- lp_aut_neg$emmean - lp_norm_neg$emmean
lp_pos_aut_norm <- lp_aut_pos$emmean - lp_norm_pos$emmean


####### between group POST-HOC ######
p_vals <- data.table()
for (l in 1:length(labels)){
  temp<- subset(lp_hp,label == labels[l])
  
  m <-  lmer(beta_power~ trial_type*group*feedback_cur+(1|subject), data = temp)
  summary(m)
  #Tuk<-data.table(summary(emmeans(regrid(ref_grid(m, transform = "response")), pairwise ~ trial_type, adjust = 'tukey',lmer.df = "satterthwaite"))$contrasts)
  Tuk <- data.table(summary(emmeans(m, pairwise ~ as.factor(group|trial_type|feedback_cur), adjust = 'tukey',lmer.df = "satterthwaite",lmerTest.limit=2498162))$contrasts)
  #Tuk[,contrast:=gsub(' - ','_',contrast)]
  Tuk[,p.value:=format(p.value, digits = 3)]
  #Tuk[,contrast:=paste0(group,'-',contrast)]
  Tuk[,contrast:=paste0(feedback_cur,'-',trial_type,'-',contrast)] ####### if you need triple interaction #########
  columns <- c('contrast','p.value')
  Tuk <- Tuk[,..columns]
  Tuk$interval <- "mean_beta"
  Tuk <- dcast(Tuk,formula = interval~contrast,value.var = 'p.value')
  Tuk$label <- unique(temp$label)
  p_vals <- rbind(p_vals,Tuk)
}
p_vals<- cbind(p_vals,lp_neg_aut_norm )
p_vals<- cbind(p_vals,lp_pos_aut_norm )
setwd("/Users/kristina/Documents/stc/lmem_label/dfs")
write.csv(p_vals, "tukey_label_1100_1500_betweengroup.csv")



####### between group post-hoc's ######
p_vals <- data.table()
for (l in 1:length(labels)){
  temp<- subset(lp_hp,label == labels[l])
  
  m <-  lmer(beta_power~ trial_type*group*feedback_cur+(1|subject), data = temp)
  summary(m)
  #Tuk<-data.table(summary(emmeans(regrid(ref_grid(m, transform = "response")), pairwise ~ trial_type, adjust = 'tukey',lmer.df = "satterthwaite"))$contrasts)
  Tuk <- data.table(summary(emmeans(m, pairwise ~ as.factor(feedback_cur|trial_type|group), adjust = 'tukey',lmer.df = "satterthwaite",lmerTest.limit=2498162))$contrasts)
  #Tuk[,contrast:=gsub(' - ','_',contrast)]
  Tuk[,p.value:=format(p.value, digits = 3)]
  #Tuk[,contrast:=paste0(group,'-',contrast)]
  Tuk[,contrast:=paste0(trial_type,'-',group,'-',contrast)] ####### if you need triple interaction #########
  columns <- c('contrast','p.value')
  Tuk <- Tuk[,..columns]
  Tuk$interval <- "mean_beta"
  Tuk <- dcast(Tuk,formula = interval~contrast,value.var = 'p.value')
  Tuk$label <- unique(temp$label)
  p_vals <- rbind(p_vals,Tuk)
}

p_vals<- cbind(p_vals,lp_aut)
p_vals<- cbind(p_vals,lp_norm)

setwd("/Users/kristina/Documents/stc/lmem_label/dfs")
write.csv(p_vals, "tukey_label_1100_1500_ingroup.csv")

