library(tidyverse)
library(kableExtra)
library(gtsummary)
library(expss)
library(haven)
library(sjlabelled)
library(readxl)
library(gtools)
library(tableone)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(hrbrthemes)
library(xgboost)
library(finalfit)
library(caret)
library(pROC) 
library(groupdata2)
library(SHAPforxgboost)
source("/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/abcd_study_with_kat/code/ABCD_cog_resilience/shap.R")

#-- Load in unimputed data
load("/Users/Kat/Dropbox/HSPH_T32/Potential projects/ABCD melatonin/temp_data/data_all.Rda")


#-- Load in imputed data
load("/Users/Kat/Dropbox/HSPH_T32/Potential projects/ABCD melatonin/temp_data/d_imp_mel_use.Rdata")


#-- Load in data dictionary - to be able to label variables
dict <- read.csv("/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/abcd_study_with_kat/temp_dat/data_dict_short_labels.csv")

setwd("/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/release5/core/")

#-- Create an ID/familyID/eventname scaffold such that each person has an observation for bl, 1, 2, 3-year follow-up events


bl <- read.csv("./abcd-general/abcd_p_demo.csv") %>%
  filter(eventname=="baseline_year_1_arm_1") %>%
  dplyr::select(src_subject_id)

age9 <- bl %>% mutate(age =  9)
age10 <- bl %>% mutate(age =  10)
age11 <- bl %>% mutate(age =  11)
age12 <- bl %>% mutate(age =  12)
age13 <- bl %>% mutate(age =  13)
age14 <- bl %>% mutate(age =  14)
age15 <- bl %>% mutate(age =  15)

scaffold <- rbind(age9,age10,age11,age12,age13,age14,age15)

#-- Demo: age
ages <- read.csv("./abcd-general/abcd_y_lt.csv") %>%
  mutate(
    age = round(interview_age/12),
    #-- Round 16-year-olds to 15 - the bin is too small
    age = ifelse(age==16,15,age)
  ) %>%
  dplyr::select(src_subject_id, eventname, age)

summary(as.factor(ages$age))

#-- Link scaffold to eventnames by age

scaffold_ages <- unique(merge(x=scaffold,y=ages,by=c("src_subject_id", "age"),all.x=T))


#-- Medications across follow-up

#-- search for dietary supplements in medication lists
mel_toMatch <- c("Melatonin","melatonin","Melatonex")
vit_toMatch <- c("Vitamin","vitamin","VITAMIN")

ph_p_meds <- read.csv("./physical-health/ph_p_meds.csv") %>%
  mutate(
    med_use = case_when(
      (eventname != "3_year_follow_up_y_arm_1" & eventname != "4_year_follow_up_y_arm_1" & brought_medications %in% c(NA_real_, 2, 4)) |
        ((eventname == "3_year_follow_up_y_arm_1" | eventname == "4_year_follow_up_y_arm_1") & brought_medications_1yr_p %in% c(NA_real_, 777, 2)) ~ NA_real_,
      (eventname != "3_year_follow_up_y_arm_1" & eventname != "4_year_follow_up_y_arm_1" & brought_medications == 3) |
        ((eventname == "3_year_follow_up_y_arm_1" | eventname == "4_year_follow_up_y_arm_1") & brought_medications_1yr_p == 1) ~ 1,
      T ~ 0),
    mel_use = apply(., 1, function(x)as.integer(any(grep(paste(mel_toMatch,collapse="|"),x)))),
    #Recode as missing if med_use missing
    mel_use = case_when(
      is.na(med_use) ~ NA_real_,
      T ~ mel_use),
    vit_use = apply(., 1, function(x)as.integer(any(grep(paste(vit_toMatch,collapse="|"),x)))),
    #Recode as missing if med_use missing
    vit_use = case_when(
      is.na(med_use) ~ NA_real_,
      T ~ vit_use),
    visit = case_when(eventname == "baseline_year_1_arm_1" ~ 0,
                      eventname == "1_year_follow_up_y_arm_1" ~ 1,
                      eventname == "2_year_follow_up_y_arm_1" ~ 2,
                      eventname == "3_year_follow_up_y_arm_1" ~ 3,
                      eventname == "4_year_follow_up_y_arm_1" ~ 4)
  ) %>%
  group_by(src_subject_id) %>%
  mutate(
    mel_use_years = cumsum(mel_use)
  ) %>%
  dplyr::select(src_subject_id,eventname,visit,med_use,mel_use,vit_use)

summary(as.factor(ph_p_meds$visit))

ph_p_meds_scaff <- merge(x=scaffold_ages, y=ph_p_meds, by=c("src_subject_id","eventname"),all.x=T) %>%
  dplyr::select(-c("eventname")) %>%
  arrange(src_subject_id,age,desc(mel_use))

ph_p_meds_scaff <- unique(ph_p_meds_scaff)

#-- Look at proportions of users by age
prop.table(table(ph_p_meds_scaff$age,ph_p_meds_scaff$mel_use), margin=1)

#-- Run a generalized linear mixed effect model to assess the association with age
#-- Square the z-stat to get a Chi-sq test (1df, bc age continuous)
library(lme4)
m <- glmer(mel_use ~ age + (1 | src_subject_id), data = ph_p_meds_scaff, 
           family = binomial, control = glmerControl(optimizer = "bobyqa"),
           nAGQ = 10)
summary(m)


#-- Transpose med use data - multiple obs per age occur - fix
ph_p_meds_wide <- ph_p_meds_scaff %>%
  pivot_wider(names_from = c(age,visit),
              values_from = c(med_use, vit_use, mel_use)) %>%
  mutate(mel_use_9 = pmax(mel_use_9_0,na.rm=T),
         mel_use_10 = pmax(mel_use_10_0,mel_use_10_1,na.rm=T),
         mel_use_11 = pmax(mel_use_11_0,mel_use_11_1,mel_use_11_2,mel_use_11_3,na.rm=T),
         mel_use_12 = pmax(mel_use_12_1,mel_use_12_2,mel_use_12_3,mel_use_12_4,na.rm=T),
         mel_use_13 = pmax(mel_use_13_2,mel_use_13_3,mel_use_13_4,na.rm=T),
         mel_use_14 = pmax(mel_use_14_2,mel_use_14_3,mel_use_14_4,na.rm=T),
         mel_use_15 = pmax(mel_use_15_3,mel_use_15_4,na.rm=T),
         mel_use = case_when(mel_use_9==1 | mel_use_10==1 | mel_use_11==1 |
                             mel_use_12==1 | mel_use_13 ==1 | 
                             mel_use_14==1 | mel_use_15==1 ~ 1,
                             is.na(mel_use_9) & is.na(mel_use_10) & is.na(mel_use_11) & is.na(mel_use_12) & is.na(mel_use_13) & is.na(mel_use_14) & is.na(mel_use_15) ~ NA_real_,
                             T ~ 0)) %>%
  dplyr::select(src_subject_id,mel_use_9,mel_use_10,
                mel_use_11,mel_use_12,mel_use_13,mel_use_14,mel_use_15,mel_use)

table(ph_p_meds_wide$mel_use)

#-- Table - comparison of covariates by any mel use, any vit use

names(data_all)

predictors <- data_all %>% 
  dplyr::select(names(data_all)[c(1,2,23:150)])

d <- merge(x=predictors,y=ph_p_meds_wide,by="src_subject_id",all=T) 
names(d)

summary(as.factor(d$mel_use))


#-- Figure - prevalence of mel use, vit use, by covariates, over time

calc_prop_est_ci <- function(v){
  est=mean(v,na.rm=T)
  se=sqrt(est*(1-est)/length(v))
  #Express in percent
  ci_l=100*(est-1.96*se)
  ci_h=100*(est+1.96*se)
  est=100*est
  return(cbind(est,se,ci_l,ci_h))
} 

#-- Overall prevalence across age

prevs <- data.frame(t(apply(d[,131:138],MARGIN=2,FUN=calc_prop_est_ci))) %>%
  tibble::rownames_to_column(.,var="age") %>%
  mutate(age = substring(age,9),
         age = case_when(
           age=="" ~ "Any 9-16",
           age=="15" ~ "15-16",
           T ~ age),
         group="All")

names(prevs) <- c("age","est","se","ci_l","ci_h","group")
prevs

#--- Plot prevalence proportions
forplot<-prevs
forplot$visit <- factor(prevs$age, levels=c("9","10","11","12","13","14","15-16",
                                          "Any 9-16"))

prev_plot <- ggplot(prevs, aes(x=factor(age, levels=c("9","10","11","12","13","14","15-16",
                                                          "Any 9-16")),y=est,group=group)) + 
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin=ci_l, ymax=ci_h), width=.1, 
                linetype="solid",position = position_dodge(0.8)) + 
  labs(x="Age", y="Prevalence (%)") +
  ylim(0,20) +
  scale_fill_ipsum() +
  theme(text=element_text(size=11),
        axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
        panel.background = element_rect(fill = 'white'),
        panel.grid = element_line(color="gray", size=0.1)) 
prev_plot



#-- Stratify by race/ethnicity

table(d$race_ethnicity)

#-- Run an Chi-sq test comparing use by sex
chisq.test(d$race_ethnicity,d$mel_use)


d_imp_w <- d %>% filter(re_nhw==1)
d_imp_b <- d %>% filter(re_nhb==1)
d_imp_h <- d %>% filter(re_h==1)
d_imp_a <- d %>% filter(re_nha==1)
d_imp_o <- d %>% filter(re_o==1)


prevs_w <- data.frame(t(apply(d_imp_w[,c(131:138)],MARGIN=2,FUN=calc_prop_est_ci))) %>%
  mutate(race_eth = "White")%>%
  tibble::rownames_to_column(.,var="age") %>%
  mutate(age = substring(age,9),
         age = case_when(
           age=="" ~ "Any 9-16",
           age=="15" ~ "15-16",
           T ~ age)
  )
prevs_b <- data.frame(t(apply(d_imp_b[,c(131:138)],MARGIN=2,FUN=calc_prop_est_ci))) %>%
  mutate(race_eth = "Black")%>%
  tibble::rownames_to_column(.,var="age") %>%
  mutate(age = substring(age,9),
         age = case_when(
           age=="" ~ "Any 9-16",
           age=="15" ~ "15-16",
           T ~ age)
  )
prevs_h <- data.frame(t(apply(d_imp_h[,c(131:138)],MARGIN=2,FUN=calc_prop_est_ci))) %>%
  mutate(race_eth = "Hispanic")%>%
  tibble::rownames_to_column(.,var="age") %>%
  mutate(age = substring(age,9),
         age = case_when(
           age=="" ~ "Any 9-16",
           age=="15" ~ "15-16",
           T ~ age)
  )
prevs_a <- data.frame(t(apply(d_imp_a[,c(131:138)],MARGIN=2,FUN=calc_prop_est_ci))) %>%
  mutate(race_eth = "Asian")%>%
  tibble::rownames_to_column(.,var="age") %>%
  mutate(age = substring(age,9),
         age = case_when(
           age=="" ~ "Any 9-16",
           age=="15" ~ "15-16",
           T ~ age)
  )
prevs_o <- data.frame(t(apply(d_imp_o[,c(131:138)],MARGIN=2,FUN=calc_prop_est_ci))) %>%
  mutate(race_eth = "Other")%>%
  tibble::rownames_to_column(.,var="age") %>%
  mutate(age = substring(age,9),
         age = case_when(
           age=="" ~ "Any 9-16",
           age=="15" ~ "15-16",
           T ~ age)
  )

prevs_re <- rbind(prevs_w,prevs_b,prevs_h,prevs_a,prevs_o)

prevs_filler <- prevs_re %>% filter(age=="9") %>%
  mutate(
    age="",
    X1=0,X2=0,X3=0,X4=0
  )

prevs_re <- rbind(prevs_re,prevs_filler)


names(prevs_re) <- c("age","est","se","ci_l","ci_h","race_eth")

prevs_re

#--- Plot prevalences of med, vit, melatonin use
prev_plot_re <- ggplot(prevs_re, aes(x=factor(age, 
                                     levels=c("9","10","11","12","13","14","15-16","","Any 9-16")), y=est, 
                                     fill=as.factor(race_eth))) + 
  geom_bar(data=prevs_re, stat="identity", position=position_dodge(), colour="gray") +
  geom_errorbar(data=prevs_re, aes(ymin=ci_l, ymax=ci_h), width=.1, 
                linetype="solid",position = position_dodge(0.8)) + 
  #geom_text(aes(label=paste0(round(est,1),"%")), position=position_dodge(width=0.9), vjust=-3) +
  labs(x="Age", y="Prevalence (%)", fill="Race/ethnicity") +
  theme(text=element_text(size=11),
        axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
        panel.background = element_rect(fill = 'white'),
        panel.grid = element_line(color="gray", size=0.1))
prev_plot_re + scale_fill_brewer(palette = "Pastel2")


#-- Figure without "other" race
prevs_re <- prevs_re %>% filter(race_eth != "Other")

#--- Plot prevalences of med, vit, melatonin use
prev_plot_re <- ggplot(prevs_re, aes(x=factor(age, 
                                              levels=c("9","10","11","12","13","14","15-16","","Any 9-16")), y=est, 
                                     fill=as.factor(race_eth))) + 
  geom_bar(data=prevs_re, stat="identity", position=position_dodge(), colour="gray") +
  geom_errorbar(data=prevs_re, aes(ymin=ci_l, ymax=ci_h), width=.1, 
                linetype="solid",position = position_dodge(0.8)) + 
  #geom_text(aes(label=paste0(round(est,1),"%")), position=position_dodge(width=0.9), vjust=-3) +
  labs(x="Age", y="Prevalence (%)", fill="Race/ethnicity") +
  theme(text=element_text(size=11),
        axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
        panel.background = element_rect(fill = 'white'),
        panel.grid = element_line(color="gray", size=0.1))
prev_plot_re + scale_fill_brewer(palette = "Pastel2")




#------------------------------------------------------------------------------#
#-- Add labels and measurement details to the predictor data - for eMethods
pred_names <- data.frame(include)
names(pred_names) <- "Variable"
predictors_labeled <- merge(x=pred_names, y=dict, by="Variable", all.x=T)
write.csv(predictors_labeled, file="/Users/Kat/Dropbox/HSPH_T32/Potential projects/ABCD melatonin/results/predictors_labeled.csv")


#------------------------------------------------------------------------------#
#-- Predict mel_use --#

names(d)
#----------------------------------------#
#-- Create folds with kids from the same family kept together
d$rel_family_id <- as.factor(d$rel_family_id)

set.seed(123)
inc <- fold(d,id_col = "rel_family_id",k = 10)
summary(inc$.folds)
foldid = as.numeric(inc$.folds)
summary(foldid)

#---------------------------------------#
#-- Variable Importance - SHAP values --#

#-- Cross-validated (10 folds)
set.seed(123)
names(inc)

#-- drop ICE(income only) - collinear with ICE(race+income)
include <- names(inc[,c(3:23,25,27:40,46:130)])
include

d_flpov <- inc %>%
  dplyr::select(all_of(include), mel_use, .folds) %>%
  filter(!is.na(mel_use)) # 11864 out of 11868 observations

summary(d_flpov$mel_use)


kfolds=10
cv_err = rep(0, kfolds)
auc = rep(0,kfolds)
shap_values = list(shap_score = d_flpov[0,c(1:length(include))], mean_shap_score=NA_real_)
for(i in 1:kfolds){
  in_data <- filter(d_flpov, .folds!=i) 
  labels = in_data$mel_use
  in_data <- as.matrix(in_data[,c(1:length(include))])
  out_data <- filter(d_flpov, .folds==i)
  out_labels = out_data$mel_use
  out_data <- as.matrix(out_data[,c(1:length(include))])
  
  dtrain <- xgb.DMatrix(data = in_data,label = labels) 
  
  params <- list(booster = "gbtree", objective = "binary:logistic", 
                 eta=0.3, gamma=0.1, max_depth=2, 
                 min_child_weight=1, 
                 subsample=1, colsample_bytree=1)
  
  fit <- xgboost(params = params, 
                 data = in_data, 
                 label=labels, 
                 nrounds = 50, 
                 maximize = F, 
                 eval="error")
  preds <- predict(fit, newdata=out_data)
  roc_test <- roc(out_labels, preds, algorithm = 2) 
  auc[i] <- roc_test$auc
  shap_result = shap.score.rank(xgb_model = fit, 
                                X_train =out_data,
                                shap_approx = F)
  
  shap_values$shap_score<-rbind(shap_values$shap_score,shap_result$shap_score)
  shap_values$mean_shap_score <- cbind(shap_values$mean_shap_score,shap_result$mean_shap_score)
  
  err <- out_labels - preds
  mse <- mean(err^2)
  # Record the RMSE
  cv_err[i] <- sqrt(mse)
}
shap_values$mean_shap_score <- rowMeans(cbind(shap_values$mean_shap_score),na.rm=T)

quantile(shap_values$mean_shap_score,c(0.7,0.8,0.9,0.95,1))
mean(auc)

## Prepare data for top N variables
shap_long = shap.prep(shap = shap_values,
                      X_train = d_flpov[,c(1:length(include))], 
                      top_n = length(include))
imp_vars_fluid <- data.frame(shap_values$mean_shap_score) %>%
  tibble::rownames_to_column(., "Variable")
imp_vars_fluid <- merge(x=imp_vars_fluid, y=dict, by="Variable",all.x=T) %>%
  dplyr::select(Short_Label,Variable,shap_values.mean_shap_score) %>%
  arrange(desc(shap_values.mean_shap_score))


names(imp_vars_fluid) <- c("Short_Label","Variable","inc_mel_SHAP")
save(imp_vars_fluid, file="/Users/Kat/Dropbox/HSPH_T32/Potential projects/ABCD melatonin/results/mel_use_SHAP_nonmissing_UPDATE_jan16.Rda")


## Plot var importance based on SHAP
var_importance(shap_values, top_n=length(include))

## Plot shap overall metrics
# plot.shap.summary(data_long = shap_long)


source("/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/abcd_study_with_kat/code/ABCD_cog_resilience/xgb_plot_shap_1.R")

xgb.plot.shap.1(data = as.matrix(d_flpov[,c(1:length(include))]), # input data
                features = c(imp_vars_fluid$Variable[c(1:10)]), 
                model=fit,
                xlab = c(imp_vars_fluid$Short_Label[c(1:10)]),
                n_col = 2, # layout option
                plot_loess = T # add red line to plot
)

imp_vars_fluid$Variable[c(1:10)]

#-- Standardize continuous variables 
inc <- inc %>%
  mutate_at(c("sds_p_ss_dims","reshist_addr1_p1vlnt","cbcl_scr_syn_thought_r",
              "asr_scr_somatic_t","crpbi_y_ss_caregiver","reshist_addr1_coi_r_he_nat"), ~(scale(.) %>% as.vector)) 

glm_cvxgboost = glm(data=inc, mel_use ~ sds_p_ss_dims+famhx_ss_moth_prob_prf_p+
                      reshist_addr1_p1vlnt+as.factor(race_ethnicity)+
                      reshist_addr1_seg_ice_inc_bw+cbcl_scr_syn_thought_r+
                      asr_scr_somatic_t+ reshist_addr1_coi_r_he_nat+crpbi_y_ss_caregiver+inr, family="binomial")

obj<-summary(glm_cvxgboost)$coefficients
obj


mel_or <- data.frame(obj)
mel_or <- tibble::rownames_to_column(mel_or, "Variable")
names(mel_or) <- c("Variable","Estimate","Std..Error","z.value","p")

mel_or <- mel_or %>%
  mutate(
    OR=round(exp(Estimate),2),
    ci_l = round(exp(Estimate - 1.96*Std..Error),2),
    ci_u = round(exp(Estimate + 1.96*Std..Error),2),
    sig = case_when(
      p < 0.0001 ~ "***",
      p < 0.01 ~ "**",
      p < 0.05 ~ "*",
      T~""
    ),
    mel_or_ci = paste0(OR," (",ci_l, ",", ci_u, ")",sig)
  ) %>%
  dplyr::select(Variable,mel_or_ci)

mel_or_label <- merge(x=mel_or, y=dict, by="Variable", all.x=T) %>%
  dplyr::select(Short_Label, Variable, mel_or_ci)

mel_or_label_shap <- merge(x=mel_or_label,y=imp_vars_fluid,by=c("Variable","Short_Label"),all.x=T)

pred <- predict(glm_cvxgboost, newdata=inc, type="response")
roc_test <- roc(inc$mel_use, pred, algorithm = 2) 
roc_test$auc

save(mel_or_label_shap,file="/Users/Kat/Dropbox/HSPH_T32/Potential projects/ABCD melatonin/results/mel_use_coef_nonmissing_update_Jan16.Rda")
write.csv(mel_or_label_shap, file="/Users/Kat/Dropbox/HSPH_T32/Potential projects/ABCD melatonin/results/mel_use_coef_nonmissing_update_Jan16.csv")


#--------------------------------------------------------------------------------#
#-- Predict initiation of use post baseline 

#-- Load in baseline mel use variable (mel_use_0)
bl_use <- data_all %>% dplyr::select(src_subject_id,mel_use_0)


#-- Merge with study data set
inc <- merge(x=inc, y=bl_use, by="src_subject_id", all.x=T)
#-- Remove those who used mel at baseline
inc <- inc %>% filter(mel_use_0 != 1)
summary(inc$mel_use)

names(inc)

#-- drop ICE(income only) - collinear with ICE(race+income)
include <- names(inc[,c(3:23,25,27:40,46:130)])
include

d_flpov <- inc %>%
  dplyr::select(all_of(include), mel_use, .folds) %>%
  filter(!is.na(mel_use)) # 11374 out of 11375 observations

summary(d_flpov$mel_use)


kfolds=10
cv_err = rep(0, kfolds)
auc = rep(0,kfolds)
shap_values = list(shap_score = d_flpov[0,c(1:length(include))], mean_shap_score=NA_real_)
for(i in 1:kfolds){
  in_data <- filter(d_flpov, .folds!=i) 
  labels = in_data$mel_use
  in_data <- as.matrix(in_data[,c(1:length(include))])
  out_data <- filter(d_flpov, .folds==i)
  out_labels = out_data$mel_use
  out_data <- as.matrix(out_data[,c(1:length(include))])
  
  dtrain <- xgb.DMatrix(data = in_data,label = labels) 
  
  params <- list(booster = "gbtree", objective = "binary:logistic", 
                 eta=0.3, gamma=0.1, max_depth=2, 
                 min_child_weight=1, 
                 subsample=1, colsample_bytree=1)
  
  fit <- xgboost(params = params, 
                 data = in_data, 
                 label=labels, 
                 nrounds = 50, 
                 maximize = F, 
                 eval="error")
  preds <- predict(fit, newdata=out_data)
  roc_test <- roc(out_labels, preds, algorithm = 2) 
  auc[i] <- roc_test$auc
  shap_result = shap.score.rank(xgb_model = fit, 
                                X_train =out_data,
                                shap_approx = F)
  
  shap_values$shap_score<-rbind(shap_values$shap_score,shap_result$shap_score)
  shap_values$mean_shap_score <- cbind(shap_values$mean_shap_score,shap_result$mean_shap_score)
  
  err <- out_labels - preds
  mse <- mean(err^2)
  # Record the RMSE
  cv_err[i] <- sqrt(mse)
}
shap_values$mean_shap_score <- rowMeans(cbind(shap_values$mean_shap_score),na.rm=T)

quantile(shap_values$mean_shap_score,c(0.7,0.8,0.9,0.95,1))
mean(auc)

## Prepare data for top N variables
shap_long = shap.prep(shap = shap_values,
                      X_train = d_flpov[,c(1:length(include))], 
                      top_n = length(include))
imp_vars_fluid <- data.frame(shap_values$mean_shap_score) %>%
  tibble::rownames_to_column(., "Variable")
imp_vars_fluid <- merge(x=imp_vars_fluid, y=dict, by="Variable",all.x=T) %>%
  dplyr::select(Short_Label,Variable,shap_values.mean_shap_score) %>%
  arrange(desc(shap_values.mean_shap_score))


names(imp_vars_fluid) <- c("Short_Label","Variable","inc_mel_SHAP")
save(imp_vars_fluid, file="/Users/Kat/Dropbox/HSPH_T32/Potential projects/ABCD melatonin/results/mel_incident_use_SHAP_nonmissing_update_Jan16.Rda")


## Plot var importance based on SHAP
var_importance(shap_values, top_n=length(include))

## Plot shap overall metrics
# plot.shap.summary(data_long = shap_long)


source("/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/abcd_study_with_kat/code/ABCD_cog_resilience/xgb_plot_shap_1.R")

xgb.plot.shap.1(data = as.matrix(d_flpov[,c(1:length(include))]), # input data
                features = c(imp_vars_fluid$Variable[c(1:10)]), 
                model=fit,
                xlab = c(imp_vars_fluid$Short_Label[c(1:10)]),
                n_col = 2, # layout option
                plot_loess = T # add red line to plot
)

imp_vars_fluid$Variable[c(1:10)]

#-- Standardize adult crime
inc <- inc %>%
  mutate_at(c("sds_p_ss_dims","reshist_addr1_p1vlnt","asr_scr_somatic_t",
              "reshist_addr1_coi_r_he_nat","parent_subst_rules"), ~(scale(.) %>% as.vector)) 

glm_cvxgboost = glm(data=inc, mel_use ~ sds_p_ss_dims+famhx_ss_moth_prob_prf_p+
                      reshist_addr1_p1vlnt+as.factor(race_ethnicity)+
                      reshist_addr1_seg_ice_inc_bw+asr_scr_somatic_t+
                      pubertal_sex_p+reshist_addr1_coi_r_he_nat+inr+parent_subst_rules, family="binomial")

obj<-summary(glm_cvxgboost)$coefficients
obj


mel_or <- data.frame(obj)
mel_or <- tibble::rownames_to_column(mel_or, "Variable")
names(mel_or) <- c("Variable","Estimate","Std..Error","z.value","p")

mel_or <- mel_or %>%
  mutate(
    OR=round(exp(Estimate),2),
    ci_l = round(exp(Estimate - 1.96*Std..Error),2),
    ci_u = round(exp(Estimate + 1.96*Std..Error),2),
    sig = case_when(
      p < 0.0001 ~ "***",
      p < 0.01 ~ "**",
      p < 0.05 ~ "*",
      T~""
    ),
    mel_or_ci = paste0(OR," (",ci_l, ",", ci_u, ")",sig)
  ) %>%
  dplyr::select(Variable,mel_or_ci)

mel_or_label <- merge(x=mel_or, y=dict, by="Variable", all.x=T) %>%
  dplyr::select(Short_Label, Variable, mel_or_ci)

pred <- predict(glm_cvxgboost, newdata=inc, type="response")
roc_test <- roc(inc$mel_use, pred, algorithm = 2) 
roc_test$auc

save(mel_or_label,file="/Users/Kat/Dropbox/HSPH_T32/Potential projects/ABCD melatonin/results/incident_mel_use_coef_nonmissing_update_Jan16.Rda")
write.csv(mel_or_label, file="/Users/Kat/Dropbox/HSPH_T32/Potential projects/ABCD melatonin/results/incident_mel_use_coef_nonmissing_update_Jan16.csv")


