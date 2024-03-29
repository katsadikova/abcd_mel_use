#----------------------------------------------------------------------------#
#--- SHAP_any_mel_imp_mult.R
#--- Date: 1/17/2024, update 3/8/2024
#----------------------------------------------------------------------------#

library(tidyverse)
library(kableExtra)
library(gtsummary)
library(expss)
library(haven)
library(sjlabelled)
library(readxl)
library(gtools)
library(tableone)
library(corrplot)
library(reshape2)
library(mlogit)
library(SuperLearner)
library(xgboost)
library(finalfit)
library(pROC) 
library(groupdata2)
library(SHAPforxgboost)
library(mice)
source("/Users/Kat/Dropbox/HSPH_T32/Potential projects/abcd_study_with_kat/code/ABCD_cog_resilience/shap.R")

#-- Load in imputed data
load("/Users/Kat/Dropbox/HSPH_T32/Potential projects/ABCD melatonin/temp_data/d_imp_mel_use_UPDATE_Mar5.Rdata")

setwd("/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/ABCD/release5/core")

#-- Load in data dictionary - to be able to label variables
dict <- read.csv("/Users/Kat/Dropbox/HSPH_T32/Potential projects/abcd_study_with_kat/temp_dat/Predictors_to_merge.csv")

names(mice::complete(Mice))


#-- Create an ID/familyID/eventname scaffold such that each person has an observation for each age
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

#-- Run an Chi-sq test comparing use by age
chisq.test(ph_p_meds_scaff$age,ph_p_meds_scaff$mel_use)


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

summary(ph_p_meds_wide$mel_use)
names(ph_p_meds_wide)

#-- Table - comparison of covariates by any mel use, any vit use

names(complete(Mice))

#-----------------------------------------#
#-- Run analyses within each imputation --#
#-----------------------------------------#

#-- Create a list for the output
shap_means <- list()

for (j in c(1:30)){
  
  #-- Pull out relevant predictor variables
  predictors <- complete(Mice,j) %>% 
    dplyr::select(names(complete(Mice))[c(1,2,23:149,152)])
  
  #s <- data.frame(summary(predictors))
  #write.csv(s, file="/Users/Kat/Dropbox/HSPH_T32/Potential projects/ABCD melatonin/results/predictor_summaries_imp_updated_Jan17.csv")
  
  #-- Merge predictor data with melatonin use data by subject id
  d <- merge(x=predictors,y=ph_p_meds_wide,by="src_subject_id",all=T) %>% mutate(
    #-- Windsorize breastfeeding (99th percentile = 36 mo)
    devhx_18_p = ifelse(devhx_18_p>36,36,devhx_18_p),
  )
  names(d)
  
  #-- Merge in baseline melatonin use & define post-baseline mel use!
  bl_use <- complete(Mice,j) %>% dplyr::select(src_subject_id,mel_use_0,mel_use_1,mel_use_2,mel_use_3,mel_use_4)
  d <- merge(x=d, y=bl_use, by="src_subject_id", all.x=T) %>%
    mutate(
      mel_use = case_when(
        mel_use_1 == 1 | mel_use_2 == 1 | mel_use_3 == 1 | mel_use_4 == 1 ~ 1,
        is.na(mel_use_1) & is.na(mel_use_2) & is.na(mel_use_3) & is.na(mel_use_4) ~ NA_real_,
        T ~ 0
      )
    )
  
  #-- Among those who did not use melatonin at baseline, who started using melatonin in follow-up?
  #-- 11044 out of 11868 obsevations
  inc <- d %>% filter(mel_use_0==0 & !is.na(mel_use)) #-- if don't have baseline information, proportion with missing mel_use increases
  summary(as.factor(inc$mel_use))
  
  
  #---------------------------------------------------------------#
  #-- Create folds with kids from the same family kept together --#
  inc$rel_family_id <- as.factor(inc$rel_family_id)
  
  set.seed(123)
  inc <- fold(inc,id_col = "rel_family_id",k = 10)
  summary(inc$.folds)
  foldid = as.numeric(inc$.folds)
  summary(foldid)
  
  #---------------------------------------#
  #-- Variable Importance - SHAP values --#
  
  #-- Cross-validated (10 folds)
  set.seed(408123)
  names(inc)
  
  #-- drop ICE(income only) - collinear with ICE(race+income)
  include <- names(inc[,c(3:23,25,27:40,46:130)])
  include
  
  #-- Add labels and measurement details to the predictor data - for eMethods
  pred_names <- data.frame(include)
  names(pred_names) <- "Variable"
  predictors_labeled <- merge(x=pred_names, y=dict, by="Variable", all.x=T)
  #write.csv(predictors_labeled, file="/Users/Kat/Dropbox/HSPH_T32/Potential projects/ABCD melatonin/results/predictors_labeled_imp.csv")
  
  d_flpov <- inc %>%
    dplyr::select(all_of(include), mel_use, .folds) %>%
    filter(!is.na(mel_use)) # 11864 out of 11868 observations
  
  summary(as.factor(d_flpov$mel_use))
  
  
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
    mutate(imp=j) %>%
    dplyr::select(Short_Label,Variable,shap_values.mean_shap_score,imp) %>%
    arrange(desc(shap_values.mean_shap_score))
  names(imp_vars_fluid) <- c("Short_Label","Variable","inc_mel_SHAP","imp")
  
  shap_means[[j]] <- imp_vars_fluid
} 


#-----------------------------------------------------------------------------------#
#-- Bind the imputations together, take average of SHAP values across imputations --#
#-----------------------------------------------------------------------------------#

shap_means_all<-data.frame(bind_rows(shap_means))
shap_means_all <- shap_means_all %>% 
  group_by(Short_Label,Variable) %>%
  summarise(shap_mean = mean(inc_mel_SHAP))  %>%
  arrange(desc(shap_mean))
save(shap_means_all, file="/Users/Kat/Dropbox/HSPH_T32/Potential projects/ABCD melatonin/results/shap_mean_incident_mel_use_30imp.Rda")
load("/Users/Kat/Dropbox/HSPH_T32/Potential projects/ABCD melatonin/results/shap_mean_incident_mel_use_30imp.Rda")
shap_means_all$Variable[c(1:10)]

#------------------------------------------------------------------------------------#
#-- Run logistic regression with top 10 predictors within each imputation, combine --#
#------------------------------------------------------------------------------------#

#--- convert to long-from data set to add the "max_problems" variable
long <- complete(Mice, action="long", include=TRUE)
long <- long %>%
  #-- Standardize continuous predictors
  dplyr::mutate(across(.cols = c(sds_p_ss_dims,reshist_addr1_p1vlnt,
                                 asr_scr_somatic_t,asr_scr_aggressive_t,
                                 screentime2_p_hours,parent_su_comm_risk,
                                 cbcl_scr_syn_anxdep_r),
                       .fns = ~scale(.x, center = T, scale = T)[, 1],
                       .names = '{.col}')) 

names(long)

#--- convert back to mids object
imp <- as.mids(long)

#-- Run logistic regression with the top 10 selected predictors within each of the 30 imputations

coefs <- list()
aucs <- list()

for (j in c(1:30)){
  
  #-- Pull out relevant predictor variables
  predictors <- complete(imp,j)
  
  
  #-- Merge predictor data with melatonin use data by subject id
  d <- data.frame(merge(x=predictors,y=ph_p_meds_wide,by="src_subject_id",all=T))
  print(names(d))
  
  glm_cvxgboost <- glm(mel_use ~ sds_p_ss_dims+as.factor(race_ethnicity)+famhx_ss_moth_prob_prf_p+
                         reshist_addr1_p1vlnt+reshist_addr1_seg_ice_inc_bw+asr_scr_somatic_t+
                         asr_scr_aggressive_t+screentime2_p_hours+parent_su_comm_risk+
                         cbcl_scr_syn_anxdep_r, family="binomial",data=d)
  
  coef <- data.frame(summary(glm_cvxgboost)$coefficients) %>% mutate(imp=j)
  coef <- tibble::rownames_to_column(coef, "Variable")
  coef
  coefs[[j]] <- coef
  
  
  pred <- predict(glm_cvxgboost, newdata=d, type="response")
  roc_test <- roc(d$mel_use, pred, algorithm = 2) 
  roc_test$auc
  aucs[[j]] <- roc_test$auc 
}
coefs_all<-data.frame(bind_rows(coefs))
names(coefs_all) <- c("Variable","Estimate","SE","z","p_value","imp")

ests_mean <-coefs_all %>%
  group_by(Variable) %>%
  summarise(n_imp=n(),
            Est_mean=mean(Estimate),
            Vw = mean(SE*SE),
            Vb = sum((Estimate - Est_mean)^2) / 20,
            SE = sqrt(Vw+Vb+Vb/n_imp),
            p = 2*pnorm(abs(Est_mean / SE), lower.tail=F)) %>%
  mutate(
    OR=round(exp(Est_mean),2),
    ci_l = round(exp(Est_mean - 1.96*SE),2),
    ci_u = round(exp(Est_mean + 1.96*SE),2),
    ci = paste0("(",ci_l, ",", ci_u, ")"),
    p=round(p,4)
  ) %>%
  dplyr::select(Variable,OR,ci,p)

ests_mean


mel_or_label <- merge(x=ests_mean, y=dict, by="Variable", all.x=T) %>%
  dplyr::select(Short_Label, Variable, OR,ci,p)



#-------------------------------------------------------------#
#-- Put Table 1a together (with SHAP)
tab1b <- merge(x=shap_means_all,y=mel_or_label,by=c("Short_Label","Variable"),all.y=T)
save(tab1b, file="/Users/Kat/Dropbox/HSPH_T32/Potential projects/ABCD melatonin/results/mel_use_tab1b_imp_update_Mar9.Rda")
write.csv(tab1b, file="/Users/Kat/Dropbox/HSPH_T32/Potential projects/ABCD melatonin/results/mel_use_tab1b_imp_update_Mar9.csv")

aucs_mean<-mean(unlist(aucs))
aucs_mean
