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
options(scipen = 999)
options(scipen = 0)

setwd("/Users/kristina/Documents/stc/lmem_label/lmem_sources")
setwd('/Users/kristina/Documents/stc/lmem_label/df_lmem_label')
data_path<- "/Users/kristina/Documents/stc/lmem_label/df_lmem_label"
out_path<- "/Users/kristina/Documents/stc/lmem_label/lmem_sources"
df<- read.csv('/Users/kristina/Documents/stc/lmem_label/tbl_new_int.csv', sep=",")
df[,1]<-NULL

##### Merge all df with label at one df #########
tbl <-list.files(pattern = "*.csv",full.names = T) 
read_plus <- function(flnm) {
  read_csv(flnm) %>% 
    mutate(filename = flnm)
}

tbl_with_sources <-list.files(pattern = "*.csv", full.names = T) %>% 
  map_df(~read_plus(.))
setDT(tbl_with_sources)
tbl_with_sources[,1]<-NULL
tbl_with_sources$trial_number<-NULL

df<- tbl_with_sources[!duplicated(tbl_with_sources),]
temp<-df
files <- list.files(pattern = ".csv") 

setDT(df)



write.csv(df, "tbl_with_sources.csv")
data <- data_frame(filename = tbl) %>% mutate(file_contents = map(tbl, ~ read_csv(file.path(data_path, .))))        
                   



files <- data.table(full_filename=list.files(data_path, pattern = '*.csv', full.names = T))
files$short_filename <- list.files(data_path, pattern = '*.csv', full.names = F)
emm_options(pbkrtest.limit = 8000)
p_vals <- data.table()
emm_options(lmerTest.limit = 2498162)
cols <- colnames(df)[grep('[0-9]',colnames(df))]
label<- unique(df$filename)

for (l in 1:length(label)){
  temp<- subset(df,filename == label[l])
  labell<- unique(temp$filename)
  for (j in cols) {
      m <- lmer(get(j) ~ trial_type * feedback_cur + (1|subject), data = temp)
      #Tuk<-data.table(summary(emmeans(regrid(ref_grid(m, transform = "response")), pairwise ~ trial_type, adjust = 'tukey',lmer.df = "satterthwaite"))$contrasts)
      Tuk <- data.table(summary(emmeans(m, pairwise ~ feedback_cur|trial_type, adjust = 'tukey',lmer.df = "satterthwaite",lmerTest.limit=2498162))$contrasts)
      #print(summary(Tuk))
      #print(Tuk)
    
      Tuk[,contrast:=gsub(' - ','_',contrast)]
      Tuk[,p.value:=format(p.value, digits = 3)]
      Tuk[,contrast:=paste0(trial_type,'-',contrast)]
      columns <- c('contrast','p.value')
      Tuk <- Tuk[,..columns]
      Tuk$interval <- j
      Tuk$interval <- gsub('beta power','',Tuk$interval)
      Tuk <- dcast(Tuk,formula = interval~contrast,value.var = 'p.value')
      Tuk$filename <- labell
      #Tuk$short_filename <- l
      p_vals <- rbind(p_vals,Tuk)
    }
}


write.csv(p_vals  "p_feedback.csv")

for (l in 1:length(label)){
  temp<- subset(df,filename == label[l])
  labell<- unique(temp$filename)
  for (j in cols) {
      m <- lmer(get(j) ~ trial_type * feedback_cur + (1|subject), data = temp)
      #Tuk<-data.table(summary(emmeans(regrid(ref_grid(m, transform = "response")), pairwise ~ trial_type, adjust = 'tukey',lmer.df = "satterthwaite"))$contrasts)
      Tuk <- data.table(summary(emmeans(m, pairwise ~ trial_type, adjust = 'tukey',lmer.df = "satterthwaite",lmerTest.limit=2498162))$contrasts)
      #print(summary(Tuk))
      #print(Tuk)
    
      #Tuk[,contrast:=gsub(' - ','_',contrast)]
      Tuk[,p.value:=format(p.value, digits = 3)]
      Tuk[,contrast:=paste0(trial_type,'-',contrast)]
      columns <- c('contrast','p.value')
      Tuk <- Tuk[,..columns]
      Tuk$interval <- j
      Tuk$interval <- gsub('beta power','',Tuk$interval)
      Tuk <- dcast(Tuk,formula = interval~contrast,value.var = 'p.value')
      Tuk$filename <- labell
      #Tuk$short_filename <- l
      p_vals <- rbind(p_vals,Tuk)
    }
}


write.csv(p_vals  "p_trial_type.csv")





