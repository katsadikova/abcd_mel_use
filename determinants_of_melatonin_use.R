#-- Determinants of melatonin use

library(tidyverse)
library(kableExtra)
library(gtsummary)
library(expss)
library(haven)
library(sjlabelled)
library(readxl)
library(gtools)
library(tableone)
library(mice)
library(HIMA)
library(corrplot)
library(reshape2)
library(mclust)
library(mlogit)
library(zoo)
## Parallel backend for foreach (also loads foreach and parallel; includes doMC)
library(doParallel)
## Reproducible parallelization
library(doRNG)


setwd("/Users/Kat/Library/CloudStorage/OneDrive-HarvardUniversity/VDI/ABCD/release5/core/")

#-- Create an ID/familyID/eventname scaffold such that each person has an observation for bl, 1, 2, 3-year follow-up events

bl <- read.csv("./abcd-general/abcd_y_lt.csv") %>%
  filter(eventname=="baseline_year_1_arm_1") %>%
  dplyr::select(src_subject_id, rel_family_id, eventname)

yr1 <- bl %>% mutate(eventname = "1_year_follow_up_y_arm_1")
yr2 <- bl %>% mutate(eventname = "2_year_follow_up_y_arm_1")
yr3 <- bl %>% mutate(eventname = "3_year_follow_up_y_arm_1")
yr4 <- bl %>% mutate(eventname = "4_year_follow_up_y_arm_1")

scaffold <- rbind(bl,yr1,yr2,yr3,yr4)

summary(as.factor(scaffold$eventname))

#-- Medications across follow-up

#-- search for dietary supplements in medication lists
toMatch <- c("Vitamin","vitamin","VITAMIN")

ph_p_meds <- read.csv("./physical-health/ph_p_meds.csv") %>%
  mutate(
    med_use = case_when(
      (eventname != "3_year_follow_up_y_arm_1" & eventname != "4_year_follow_up_y_arm_1" & brought_medications %in% c(NA_real_, 2, 4)) |
      ((eventname == "3_year_follow_up_y_arm_1" | eventname == "4_year_follow_up_y_arm_1") & brought_medications_1yr_p %in% c(NA_real_, 777, 2)) ~ NA_real_,
      (eventname != "3_year_follow_up_y_arm_1" & eventname != "4_year_follow_up_y_arm_1" & brought_medications == 3) |
      ((eventname == "3_year_follow_up_y_arm_1" | eventname == "4_year_follow_up_y_arm_1") & brought_medications_1yr_p == 1) ~ 1,
      T ~ 0),
    mel_use = apply(., 1, function(x)as.integer(any(grep("Melatonin",x)))),
    #Recode as missing if med_use missing
    mel_use = case_when(
      is.na(med_use) ~ NA_real_,
      T ~ mel_use),
    vit_use = apply(., 1, function(x)as.integer(any(grep(paste(toMatch,collapse="|"),x)))),
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
  dplyr::select(src_subject_id,eventname,visit,med_use,mel_use,vit_use,mel_use_years)

prop.table(table(ph_p_meds$visit,ph_p_meds$med_use,useNA="ifany"),margin=1)
prop.table(table(ph_p_meds$visit,ph_p_meds$vit_use,useNA="ifany"),margin=1)
prop.table(table(ph_p_meds$visit,ph_p_meds$mel_use,useNA="ifany"),margin=1)


#-- Merge long-form meds data onto scaffold such that we have a better idea about who used/didn't use/is missing med use data
ph_p_meds_scaf <- merge(scaffold,ph_p_meds,by=c("src_subject_id","eventname"),all.x=T) %>%
  mutate(
    visit = case_when(eventname == "baseline_year_1_arm_1" ~ 0,
                      eventname == "1_year_follow_up_y_arm_1" ~ 1,
                      eventname == "2_year_follow_up_y_arm_1" ~ 2,
                      eventname == "3_year_follow_up_y_arm_1" ~ 3,
                      eventname == "4_year_follow_up_y_arm_1" ~ 4)
  ) %>%
  arrange(src_subject_id,visit) %>%
  dplyr::select(-c(eventname)) 

#-- Transpose med use data
ph_p_meds_wide <- ph_p_meds_scaf %>%
  pivot_wider(names_from = visit,
              values_from = c(med_use, vit_use, mel_use, mel_use_years)) 

#----------------------------------------------------------------------------------------
#-- Pull together predictor data set

#-- Anthropometrics - BMI
ph_y_anthro <- read.csv("./physical-health/ph_y_anthro.csv") %>%
  mutate(
    #-- Correct 2 height measurements - in feet rather than inches
    anthroheightcalc = ifelse(anthroheightcalc<5,anthroheightcalc*12,anthroheightcalc),
    #- BMI calculated from pounds and inches: weight (lb) / [height (in)]2 x 703
    bmi = anthroweightcalc / (anthroheightcalc^2) * 703
  ) %>%
  #-- Select baseline visit
  filter(eventname=="baseline_year_1_arm_1") %>%
  dplyr::select(src_subject_id,anthroweightcalc,anthroheightcalc, bmi) # a handful of kiddos with very low weight

summary(ph_y_anthro$bmi)

#-- School functioning
ce_y_srpf <- read.csv("./culture-environment/ce_y_srpf.csv") %>%
  #-- Select baseline visit
  filter(eventname=="baseline_year_1_arm_1") %>%
  dplyr::select(src_subject_id, srpf_y_ss_ses)

#-- Adversity - Neighborhood safety and crime
ce_p_nsc <- read.csv("./culture-environment/ce_p_nsc.csv") %>%
  #-- Select baseline visit
  filter(eventname=="baseline_year_1_arm_1") %>%
  dplyr::select(src_subject_id, nsc_p_ss_mean_3_items)

#-- Adversity - Family environment - conflict scale, parent reported
ce_p_fes <- read.csv("./culture-environment/ce_p_fes.csv") %>%
  #-- Select baseline visit
  filter(eventname=="baseline_year_1_arm_1") %>%
  dplyr::select(src_subject_id, fes_p_ss_fc)

#-- Adversity - Family environment - conflict scale, child reported
ce_y_fes <- read.csv("./culture-environment/ce_y_fes.csv") %>%
  #-- Select baseline visit
  filter(eventname=="baseline_year_1_arm_1") %>%
  dplyr::select(src_subject_id, fes_y_ss_fc)

#-- Adversity - Emotional Neglect - Parent & secondary caregiver behavior
ce_y_crpbi <- read.csv("./culture-environment/ce_y_crpbi.csv") %>%
  #-- Select baseline visit
  filter(eventname=="baseline_year_1_arm_1") %>%
  dplyr::select(src_subject_id, crpbi_y_ss_parent,crpbi_y_ss_caregiver)

#-- Adversity - Physical Neglect -Parental monitoring
ce_y_pm <- read.csv("./culture-environment/ce_y_pm.csv") %>%
  #-- Select baseline visit
  filter(eventname=="baseline_year_1_arm_1") %>%
  dplyr::select(src_subject_id, pmq_y_ss_mean)


#-- COI
led_l_coi <- read.csv("./linked-external-data/led_l_coi.csv") %>%
  #-- Select baseline visit
  filter(eventname=="baseline_year_1_arm_1") %>%
  dplyr::select(src_subject_id, reshist_addr1_coi_r_ed_nat, 
                reshist_addr1_coi_r_he_nat, reshist_addr1_coi_r_se_nat)


#-- Urbanization
led_l_denspop <- read.csv("./linked-external-data/led_l_denspop.csv") %>%
  #-- Select baseline visit
  filter(eventname=="baseline_year_1_arm_1") %>%
  dplyr::select(src_subject_id,reshist_addr1_popdensity)

led_l_urban <- read.csv("./linked-external-data/led_l_urban.csv") %>%
  #-- Select baseline visit
  filter(eventname=="baseline_year_1_arm_1") %>%
  dplyr::select(src_subject_id,reshist_addr1_urban_area)

led_l_walk <- read.csv("./linked-external-data/led_l_walk.csv") %>%
  #-- Select baseline visit
  filter(eventname=="baseline_year_1_arm_1") %>%
  dplyr::select(src_subject_id,reshist_addr1_walkindex)

led_l_traffic <- read.csv("./linked-external-data/led_l_traffic.csv") %>%
  #-- Select baseline visit
  filter(eventname=="baseline_year_1_arm_1") %>%
  dplyr::select(src_subject_id,reshist_addr1_traffic_count)

led_l_roadprox <- read.csv("./linked-external-data/led_l_roadprox.csv") %>%
  #-- Select baseline visit
  filter(eventname=="baseline_year_1_arm_1") %>%
  dplyr::select(src_subject_id,reshist_addr1_proxrd)

#-- Neighborhood quality: Social Vulnerability Index (SVI) - overall
led_l_svi <- read.csv("./linked-external-data/led_l_svi.csv") %>%
  #-- Select baseline visit
  filter(eventname=="baseline_year_1_arm_1") %>%
  dplyr::select(src_subject_id,reshist_addr1_svi_tot_20142018)


#-- Neightborhood quality: crime 
led_l_crime <- read.csv("./linked-external-data/led_l_crime.csv") %>%
  #-- Select baseline visit
  filter(eventname=="baseline_year_1_arm_1") %>%
  dplyr::select(src_subject_id,reshist_addr1_p1vlnt)

#-- Neighborhood quality: lead risk
led_l_leadrisk <- read.csv("./linked-external-data/led_l_leadrisk.csv") %>%
  #-- Select baseline visit
  filter(eventname=="baseline_year_1_arm_1") %>%
  dplyr::select(src_subject_id,reshist_addr1_leadrisk)


#-- Natural environment: Natural Space and Satellite, Land Use: Nighttime lights
led_l_urbsat <- read.csv("./linked-external-data/led_l_urbsat.csv") %>%
  #-- Select baseline visit
  filter(eventname=="baseline_year_1_arm_1") %>%
  dplyr::select(src_subject_id,reshist_addr1_urbsat_ntl)


#-- Index of concentration at the extremes
led_l_ice <- read.csv("./linked-external-data/led_l_ice.csv") %>%
  #-- Select baseline visit
  filter(eventname=="baseline_year_1_arm_1") %>%
  dplyr::select(src_subject_id,reshist_addr1_seg_ice_income,
                reshist_addr1_seg_ice_inc_bw)


#-- Laws and biases: ACA Medicaid expansion
led_l_aca <- read.csv("./linked-external-data/led_l_aca.csv") %>%
  #-- Select baseline visit
  filter(eventname=="baseline_year_1_arm_1") %>%
  dplyr::select(src_subject_id, reshist_addr1_ACAexpand)

#-- Demo: hhsize, hhincome, inr, parent education, marital status, family experiences
abcd_p_demo <- read.csv("./abcd-general/abcd_p_demo.csv") %>%
  #-- Select baseline visit
  filter(eventname=="baseline_year_1_arm_1") %>%
  mutate(
      #-- Recode household size < 2 to 2; 
      hhsize = ifelse(demo_roster_v2<2,2,demo_roster_v2),
      #-- Set midpoints to household income ranges
      hhinc_mid = case_when(
        demo_comb_income_v2 == 1 ~ 2500,
        demo_comb_income_v2 == 2 ~ 8500,
        demo_comb_income_v2 == 3 ~ 14000,
        demo_comb_income_v2 == 4 ~ 20500,
        demo_comb_income_v2 == 5 ~ 30000,
        demo_comb_income_v2 == 6 ~ 42000,
        demo_comb_income_v2 == 7 ~ 62500,
        demo_comb_income_v2 == 8 ~ 87500,
        demo_comb_income_v2 == 9 ~ 150000,
        demo_comb_income_v2 == 10 ~ 200000,
        demo_comb_income_v2 %in% c(777,999) ~ NA_real_),
      #-- Income-to-needs ratio based on US Department of Health and Human Services poverty thresholds (https://aspe.hhs.gov/topics/poverty-economic-mobility/poverty-guidelines/prior-hhs-poverty-guidelines-federal-register-references/2018-poverty-guidelines)
      inr = case_when(
        hhsize == 2 ~ hhinc_mid / 16460,
        hhsize == 3 ~ hhinc_mid / 20780,
        hhsize == 4 ~ hhinc_mid / 25100,
        hhsize == 5 ~ hhinc_mid / 29420,
        hhsize == 6 ~ hhinc_mid / 33740,
        hhsize == 7 ~ hhinc_mid / 38060,
        hhsize == 8 ~ hhinc_mid / 42380,
        hhsize >= 9 ~ hhinc_mid / (42380 + 4320*(hhsize-8)),
        is.na(hhsize) & is.na(hhinc_mid) ~ NA_real_
      ),
      #-- Highest level of parental education - if no partner, just primary parent education
      #-- First recode 777 & 999 to NA_real_
      demo_prnt_ed_v2 = case_when(demo_prnt_ed_v2==13 ~ 12,
                                  demo_prnt_ed_v2==14 ~ 12,
                                  demo_prnt_ed_v2==15 ~ 14,
                                  demo_prnt_ed_v2==16 ~ 14,
                                  demo_prnt_ed_v2==17 ~ 14,
                                  demo_prnt_ed_v2==18 ~ 16,
                                  demo_prnt_ed_v2==19 ~ 18,
                                  demo_prnt_ed_v2==20 ~ 20,
                                  demo_prnt_ed_v2==21 ~ 22,
                                  demo_prnt_ed_v2==777 ~ NA_real_,
                                  T~demo_prnt_ed_v2),
      demo_prtnr_ed_v2 = case_when(demo_prtnr_ed_v2==13 ~ 12,
                                   demo_prtnr_ed_v2==14 ~ 12,
                                   demo_prtnr_ed_v2==15 ~ 14,
                                   demo_prtnr_ed_v2==16 ~ 14,
                                   demo_prtnr_ed_v2==17 ~ 14,
                                   demo_prtnr_ed_v2==18 ~ 16,
                                   demo_prtnr_ed_v2==19 ~ 18,
                                   demo_prtnr_ed_v2==20 ~ 20,
                                   demo_prtnr_ed_v2==21 ~ 22,
                                   demo_prtnr_ed_v2 %in% c(777,999) ~ NA_real_,
                                   T~demo_prtnr_ed_v2),
    demo_mar_livingwith = case_when(demo_prnt_marital_v2 %in% c(1,6) ~ 1,
                                    demo_prnt_marital_v2 == 777 ~ NA_real_,
                                    T~0),
    #-- Recode missing partner education with parent ed - because will take max
    demo_prtnr_ed_v2 = case_when(is.na(demo_prtnr_ed_v2) ~ demo_prnt_ed_v2,
                                 T ~ demo_prtnr_ed_v2),
    #-- Take maximum of parent and partner education
    parent_ed = pmax(demo_prnt_ed_v2,demo_prtnr_ed_v2),
    #-- Create race/eth indicators (first recode 2 ppl with missing as "Other")
    race_ethnicity=ifelse(is.na(race_ethnicity),5,race_ethnicity),
    re_nhw=ifelse(race_ethnicity==1,1,0), 
    re_nhb=ifelse(race_ethnicity==2,1,0), 
    re_h=ifelse(race_ethnicity==3,1,0), 
    re_nha=ifelse(race_ethnicity==4,1,0), 
    re_o=ifelse(race_ethnicity==5,1,0),
    #-- Recode missing
    across(c(demo_fam_exp1_v2,demo_fam_exp2_v2,
             demo_fam_exp3_v2, demo_fam_exp4_v2,
             demo_fam_exp5_v2,demo_fam_exp6_v2,demo_fam_exp7_v2),
           ~ dplyr::recode(.,`0`=0,`1`=1,`777`=NA_real_,.default = NaN))
  ) %>%
  dplyr::select(src_subject_id, demo_brthdat_v2, 
                hhsize, hhinc_mid, inr, 
                parent_ed,
                demo_mar_livingwith, 
                demo_fam_exp1_v2,demo_fam_exp2_v2,
                demo_fam_exp3_v2,demo_fam_exp4_v2,
                demo_fam_exp5_v2,demo_fam_exp6_v2,
                demo_fam_exp7_v2, race_ethnicity,
                re_nhw, re_nhb, re_h, re_nha, re_o)

#-- Adversity - Physical abuse - broad psychopathology, KSADS - post-traumatic stress disorder (indiv. Questions)
mh_p_ksads_ptsd <- read.csv("./mental-health/mh_p_ksads_ptsd.csv") %>%
  #-- Select baseline visit
  filter(eventname=="baseline_year_1_arm_1") %>%
  mutate(
    ksads_ptsd_accident = case_when(ksads_ptsd_raw_754_p==1 | ksads_ptsd_raw_755_p==1 ~ 1,
                                    T ~ 0),
    ksads_ptsd_wit_viol = case_when(ksads_ptsd_raw_756_p==1 | ksads_ptsd_raw_757_p==1 | ksads_ptsd_raw_758_p==1 | 
                                      ksads_ptsd_raw_759_p==1 | ksads_ptsd_raw_760_p==1 | ksads_ptsd_raw_766_p==1 ~ 1,
                                    T ~ 0),
    ksads_ptsd_thr_viol = case_when(ksads_ptsd_raw_764_p==1 | ksads_ptsd_raw_765_p==1 ~ 1,
                                    T ~ 0),
    ksads_ptsd_phys_abuse = case_when(ksads_ptsd_raw_761_p==1 | ksads_ptsd_raw_762_p==1 | ksads_ptsd_raw_763_p==1 ~ 1,
                                      T ~ 0),
    ksads_ptsd_sexual_abuse = case_when(ksads_ptsd_raw_767_p==1 | ksads_ptsd_raw_768_p==1 | ksads_ptsd_raw_769_p==1 ~ 1,
                                        T ~ 0)
  ) %>%
  dplyr::select(src_subject_id, ksads_ptsd_accident, 
                ksads_ptsd_wit_viol, ksads_ptsd_thr_viol,
                ksads_ptsd_phys_abuse, ksads_ptsd_sexual_abuse)

#-- Bullying
mh_p_ksads_bg <- read.csv("./mental-health/mh_p_ksads_bg.csv") %>%
  #-- Select baseline visit
  filter(eventname=="baseline_year_1_arm_1") %>%
  mutate(kbi_p_c_bully = case_when(kbi_p_c_bully==2~0,
                                   T~kbi_p_c_bully)) %>%
  dplyr::select(src_subject_id, kbi_p_c_bully)


#-- Adversity - Family history - parent MH & substance use
mh_p_fhx <- read.csv("./mental-health/mh_p_fhx.csv") %>%
  #-- Select baseline visit
  filter(eventname=="baseline_year_1_arm_1") %>%
  dplyr::select(src_subject_id, 
                #-- problem with alcohol / drugs endorsed for mother or father
                famhx_ss_fath_prob_alc_p, famhx_ss_moth_prob_alc_p,
                famhx_ss_moth_prob_dg_p, famhx_ss_fath_prob_dg_p,
                #-- mother/father been to counselor due to emotional / mental problem
                famhx_ss_fath_prob_prf_p, famhx_ss_moth_prob_prf_p,
                #-- mother/father hospitalized due to emotional / mental problem
                famhx_ss_moth_prob_hspd_p, famhx_ss_fath_prob_hspd_p)


#-- Adversity - Household substance use & parental psychopathology - Adult Self-report (parent) ASEBA
mh_p_asr <- read.csv("./mental-health/mh_p_asr.csv") %>%
  #-- Select baseline visit
  filter(eventname=="baseline_year_1_arm_1") %>%
  dplyr::select(src_subject_id, 
                #Household substance use
                asr_q06_p, asr_q90_p, asr_q126_p,
                #Psychopathology
                asr_scr_anxdep_t, asr_scr_withdrawn_t,
                asr_scr_somatic_t, asr_scr_thought_t,
                asr_scr_attention_t, asr_scr_aggressive_t,
                asr_scr_rulebreak_t, asr_scr_intrusive_t
  )


#-- Pubertal development
ph_p_pds <- read.csv("./physical-health/ph_p_pds.csv") %>%
  #-- Select baseline visit
  filter(eventname=="baseline_year_1_arm_1") %>%
  mutate(
    pds_p_ss_cat = case_when(pubertal_sex_p==1 ~ pds_p_ss_male_category,
                             pubertal_sex_p==2 ~ pds_p_ss_female_category)
  ) %>%
  dplyr::select(src_subject_id, pubertal_sex_p, pds_p_ss_cat)


#-- Physical health - Development History Questionnaire
ph_p_dhx <- read.csv("./physical-health/ph_p_dhx.csv") %>%
  #-- Select baseline visit
  filter(eventname=="baseline_year_1_arm_1")  %>%
  mutate(
    birth_weight = birth_weight_lbs+birth_weight_oz/16,
    devhx_4_p = case_when(
      devhx_4_p==332 ~ 32,
      devhx_4_p==389 ~ 39,
      T ~ devhx_4_p
    ),
    birth_complications = case_when(
      devhx_14a3_p==1 | devhx_14b3_p==1 | devhx_14c3_p==1 | devhx_14d3_p==1 | devhx_14e3_p==1 | devhx_14f3_p==1 | devhx_14g3_p==1 | devhx_14h3_p==1 ~ 1,
      T ~ 0
    ),
    devhx_16_p = case_when(devhx_16_p==9990 ~ NA_real_,
                           T~devhx_16_p),
    late_motor_dev = case_when(devhx_20_p ==4 | devhx_20_p ==5 ~ 1,
                               T ~ 0),
    late_speech_dev = case_when(devhx_21_p ==4 | devhx_21_p ==5 ~ 1,
                                T ~ 0),
    across(c(devhx_8_prescript_med,devhx_9_prescript_med,devhx_6_p,devhx_10,devhx_12a_p,devhx_22_3_p),
           ~ dplyr::recode(.,`999`=NA_real_,`1`=1,`0`=0, `-1`=NA_real_, .default = NaN))) %>%
  dplyr::select(src_subject_id, birth_weight, devhx_3_p, devhx_4_p,
                devhx_6_p, devhx_10, 
                birth_complications,devhx_12a_p,
                devhx_16_p,devhx_17_p,devhx_18_p,late_motor_dev,
                late_speech_dev,devhx_22_3_p)


#-- Physical activity - self-reported
ph_y_yrb <- read.csv("./physical-health/ph_y_yrb.csv") %>%
  #-- Select baseline visit
  filter(eventname=="baseline_year_1_arm_1") %>%
  dplyr::select(src_subject_id, physical_activity1_y,physical_activity5_y) 


#-- Sports and organized activities
ph_p_saiq <- read.csv("./physical-health/ph_p_saiq.csv") %>%
  mutate(
    sai_p_perf_arts = case_when(sai_p_activities___0==1 | sai_p_activities___23==1 | sai_p_activities___25==1 ~ 1,
                                T~0),
    sai_p_team_sports = case_when(sai_p_activities___1==1 | sai_p_activities___2==1 |
                                    sai_p_activities___4==1 | sai_p_activities___5==1 |
                                    sai_p_activities___7==1 | sai_p_activities___8==1 |
                                    sai_p_activities___11==1 | sai_p_activities___12==1 |
                                    sai_p_activities___15==1 | sai_p_activities___17==1 | 
                                    sai_p_activities___21==1 ~ 1,
                                  T~0),
    sai_p_non_team_sports = case_when(sai_p_activities___3==1 | sai_p_activities___6==1 |
                                        sai_p_activities___9==1 | sai_p_activities___10==1 |
                                        sai_p_activities___7==1 | sai_p_activities___8==1 |
                                        sai_p_activities___13==1 | sai_p_activities___14==1 |
                                        sai_p_activities___16==1 | sai_p_activities___18==1 | 
                                        sai_p_activities___19==1 | sai_p_activities___20==1 | sai_p_activities___22==1 ~ 1,
                                      T~0),
    sai_p_other = case_when(sai_p_activities___26==1 | sai_p_activities___27 ==1 | sai_p_activities___28==1 | sai_p_activities___29==1 ~ 1,
                            T~0),
    sai_p_any = case_when(sai_p_perf_arts==1 | sai_p_team_sports ==1 | sai_p_non_team_sports==1 | sai_p_other==1 ~ 1,
                          T~0)
  ) %>%
  #-- Select baseline visit
  filter(eventname=="baseline_year_1_arm_1") %>%
  dplyr::select(src_subject_id, 
                sai_p_activities___24, sai_p_perf_arts,
                sai_p_team_sports,sai_p_non_team_sports,
                sai_p_other,sai_p_any) 


#-- Sleep
ph_p_sds <- read.csv("./physical-health/ph_p_sds.csv") %>%
  #-- Select baseline visit
  filter(eventname=="baseline_year_1_arm_1") %>%
  dplyr::select(src_subject_id,sds_p_ss_dims,sds_p_ss_sbd,
                sds_p_ss_da,sds_p_ss_swtd,sds_p_ss_does,
                sds_p_ss_shy,sds_p_ss_total)


#-- Screen time
nt_p_stq <- read.csv("./novel-technologies/nt_p_stq.csv") %>%
  #-- Select baseline visit
  filter(eventname=="baseline_year_1_arm_1") %>%
  mutate(screentime1_p_hours = screentime1_p_hours + screentime1_p_minutes/60,
         screentime2_p_hours = screentime2_p_hours + screentime2_p_minutes/60) %>%
  dplyr::select(src_subject_id,screentime1_p_hours,
                screentime2_p_hours)


#-- Temperament/Personality: Prosocial behavior - parent reported
ce_p_psb <- read.csv("./culture-environment/ce_p_psb.csv") %>%
  #-- Select baseline visit
  filter(eventname=="baseline_year_1_arm_1") %>%
  dplyr::select(src_subject_id, psb_p_ss_mean)


#-- Temperament/Personality: Prosocial behavior - self reported
ce_y_psb <- read.csv("./culture-environment/ce_y_psb.csv") %>%
  #-- Select baseline visit
  filter(eventname=="baseline_year_1_arm_1") %>%
  dplyr::select(src_subject_id, psb_y_ss_mean)

#-- Tempoerament/personality: BIS/BAS
mh_y_bisbas <- read.csv("./mental-health/mh_y_bisbas.csv") %>%
  #-- Select baseline visit
  filter(eventname=="baseline_year_1_arm_1") %>%
  dplyr::select(src_subject_id,bis_y_ss_bas_drive,bis_y_ss_bas_fs,
                bis_y_ss_bas_rr,bis_y_ss_bis_sum)

#-- Tempoerament/personality: BIS/BAS
mh_y_upps <- read.csv("./mental-health/mh_y_upps.csv") %>%
  #-- Select baseline visit
  filter(eventname=="baseline_year_1_arm_1") %>%
  dplyr::select(src_subject_id,upps_y_ss_positive_urgency,
                upps_y_ss_negative_urgency, upps_y_ss_lack_of_planning,
                upps_y_ss_lack_of_perseverance, upps_y_ss_sensation_seeking)


#-- Child psychopathology (parent-reported CBCL)
mh_p_cbcl <- read.csv("./mental-health/mh_p_cbcl.csv") %>%
  #-- Select baseline visit
  filter(eventname=="baseline_year_1_arm_1") %>%
  dplyr::select(src_subject_id,
                cbcl_scr_syn_anxdep_r,cbcl_scr_syn_withdep_r,
                cbcl_scr_syn_somatic_r,cbcl_scr_syn_social_r,
                cbcl_scr_syn_thought_r,cbcl_scr_syn_attention_r,
                cbcl_scr_syn_rulebreak_r,cbcl_scr_syn_aggressive_r)


#-- Number of friends / close friends
mh_y_or <- read.csv("./mental-health/mh_y_or.csv") %>%
  #-- Select baseline visit
  filter(eventname=="baseline_year_1_arm_1") %>%
  mutate(
    num_friends = resiliency5a_y + resiliency6a_y,
    num_close_friends = resiliency5b_y + resiliency6b_y
  ) %>%
  dplyr::select(src_subject_id,num_friends,num_close_friends)


#-- Community risk and protective factors (parent-reported)
su_p_crpf <- read.csv("./substance-use/su_p_crpf.csv") %>%
  #-- Select baseline visit
  filter(eventname=="baseline_year_1_arm_1") %>% 
  mutate(across(su_risk_p_1:su_risk_p_13, ~ dplyr::recode(.,`0`=0,`1`=1,`2`=2,`3`=3,`4`=NA_real_,.default = NaN))) %>%
  mutate(
    parent_su_comm_risk = rowSums(select(., starts_with("su_risk_p")), na.rm = TRUE)
  ) %>%
  dplyr::select(src_subject_id,parent_su_comm_risk)

#-- Rules around substance use (parent-reported)
su_p_pr <- read.csv("./substance-use/su_p_pr.csv") %>%
  #-- Select baseline visit
  filter(eventname=="baseline_year_1_arm_1") %>%
  mutate(across(c(parent_rules_q1, parent_rules_q4, parent_rules_q7, parent_rules_q_10),
                ~ dplyr::recode(.,`6`=0,`1`=1,`2`=2,`3`=3,`4`=4, `5`=5,.default = NaN))) %>%
  mutate(
    parent_subst_rules = rowSums(select(., c("parent_rules_q1", "parent_rules_q4", "parent_rules_q7", "parent_rules_q_10")), na.rm = TRUE)
  ) %>%
  dplyr::select(src_subject_id,parent_subst_rules)

#-- Peer Deviance (youth-reported)
su_y_peerdevia <- read.csv("./substance-use/su_y_peerdevia.csv") %>%
  #-- Select baseline visit
  filter(eventname=="baseline_year_1_arm_1") %>%
  mutate(
    peer_deviance = rowSums(select(., c("peer_deviance_1_4bbe5d", "peer_deviance_2_dd1457", "peer_deviance_3_e1ec2e",
                                        "peer_deviance_4_b6c588", "peer_deviance_5_bffa44", "peer_deviance_6_69562e", 
                                        "peer_deviance_7_beb683", "peer_deviance_8_35702e", "peer_deviance_9_6dd4ef")), na.rm = TRUE),
    some_peer_deviance = ifelse(peer_deviance>0,1,0)
  ) %>%
  dplyr::select(src_subject_id, peer_deviance, some_peer_deviance)


#-- Intention to use substances (alcohol/tobacco/marijuana) (youth-reported)
su_y_path_intuse <- read.csv("./substance-use/su_y_path_intuse.csv") %>%
  #-- Select baseline visit
  filter(eventname=="baseline_year_1_arm_1") %>%
  mutate(across(c(path_alc_youth1:path_alc_youth9),
                ~ dplyr::recode(.,`6`=NA_real_,`5`=NA_real_,`1`=3,`2`=2,`3`=1,`4`=0, .default = NaN))) %>%
  mutate(
    curious_subst = rowSums(select(., c("path_alc_youth1", "path_alc_youth2", "path_alc_youth3")), na.rm = TRUE),
    tobacco_intension = rowSums(select(., c("path_alc_youth1","path_alc_youth4","path_alc_youth7")), na.rm=TRUE),
    alc_intension = rowSums(select(., c("path_alc_youth2","path_alc_youth5","path_alc_youth8")), na.rm=TRUE),
    mj_intension = rowSums(select(., c("path_alc_youth3","path_alc_youth6","path_alc_youth9")), na.rm=TRUE)
  ) %>%
  dplyr::select(src_subject_id, curious_subst, tobacco_intension, alc_intension, mj_intension)



#-------------------------------------------------------------------------#
#-- Assemble the data set
data_list <- list()
data_list <- append(data_list, list(ph_p_meds_wide,ce_y_srpf,ce_p_nsc,
                                    ce_p_fes,ce_y_fes,ce_y_crpbi,
                                    ce_y_pm, ce_p_psb, ce_y_psb,led_l_coi,
                                    led_l_denspop, led_l_urban, led_l_walk,
                                    led_l_traffic, led_l_roadprox, led_l_svi,
                                    led_l_crime, led_l_leadrisk, led_l_urbsat,
                                    led_l_ice,led_l_aca, abcd_p_demo, mh_p_ksads_ptsd,
                                    mh_p_ksads_bg,mh_p_fhx, mh_p_asr, ph_p_pds,
                                    ph_p_dhx,ph_y_yrb,ph_p_saiq,
                                    ph_p_sds,nt_p_stq,mh_y_bisbas,mh_y_upps,
                                    mh_p_cbcl,mh_y_or,su_p_crpf,
                                    su_p_pr,su_y_peerdevia, su_y_path_intuse,
                                    ph_y_anthro))
length(data_list)

data_all <- Reduce(function(x, y) merge(x, y, by="src_subject_id", all=TRUE), data_list)


#-- Save the analytic data set
#save(data_all, file="/Users/Kat/Dropbox/HSPH_T32/Potential projects/ABCD melatonin/temp_data/data_all.Rda")
load("/Users/Kat/Dropbox/HSPH_T32/Potential projects/ABCD melatonin/temp_data/data_all.Rda")


#------------------------------------#
#-- Impute missing predictor data  --#
#--------------------------------------------------------------------#
#-- Imputation code source: https://rpubs.com/kaz_yos/mice-exclude --#
table(complete.cases(data_all))
propmiss <- function(dataframe) lapply(dataframe,function(x) data.frame(nmiss=sum(is.na(x)), n=length(x), propmiss=sum(is.na(x))/length(x)))
misses <-propmiss(data_all) # bmi missing for 23% of obs
misses <- do.call(rbind.data.frame, misses) 
write.csv(misses, file="/Users/Kat/Dropbox/HSPH_T32/Potential projects/ABCD melatonin/temp_data/misses.csv")


## Configure parallelization
## Detect core count
nCores <- min(parallel::detectCores(), 8)
## Used by parallel::mclapply() as default
options(mc.cores = nCores)
## Used by doParallel as default
options(cores = nCores)
## Register doParallel as the parallel backend with foreach
## http://stackoverflow.com/questions/28989855/the-difference-between-domc-and-doparallel-in-r
doParallel::registerDoParallel(cores = nCores)
## Report multicore use
cat("### Using", foreach::getDoParWorkers(), "cores\n")
cat("### Using", foreach::getDoParName(), "as backend\n")

## Extract all variable names in dataset
allVars <- names(data_all)
allVars

## names of variables with missingness
missVars <- names(data_all)[colSums(is.na(data_all)) > 0]
missVars

pred <- quickpred(data_all, mincor = 0.2)
meth <- make.method(data = data_all)

# Don't impute mel use variables (or other med/vit use variables) - only impute predictors
pred[,c("med_use_0", "med_use_1", "med_use_2", "med_use_3", "med_use_4", 
        "vit_use_0", "vit_use_1", "vit_use_2", "vit_use_3", "vit_use_4",
        "mel_use_0", "mel_use_1", "mel_use_2", "mel_use_3", "mel_use_4"  )] <- 0 
meth[c("med_use_0", "med_use_1", "med_use_2", "med_use_3", "med_use_4", 
       "vit_use_0", "vit_use_1", "vit_use_2", "vit_use_3", "vit_use_4",
       "mel_use_0", "mel_use_1", "mel_use_2", "mel_use_3", "mel_use_4")] <- ""

Mice <- mice(data = data_all, m = 30, predictorMatrix = pred, method = meth, maxit = 5, seed=123)

#-- Save Imputed Data!
#save(Mice, file="/Users/Kat/Dropbox/HSPH_T32/Potential projects/ABCD melatonin/temp_data/d_imp_mel_use_UPDATE_Mar5.Rdata")
load("/Users/Kat/Dropbox/HSPH_T32/Potential projects/ABCD melatonin/temp_data/d_imp_mel_use_UPDATE_Mar5.Rdata")

dim(complete(Mice))

summary(complete(Mice)$mel_use_0) 
summary(data_all$mel_use_0)

names(complete(Mice))
mean(complete(Mice)$demo_brthdat_v2)
sqrt(var(complete(Mice)$demo_brthdat_v2))
