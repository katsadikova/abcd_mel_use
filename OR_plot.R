library(tidyverse)
library(ggplot2)
library(ggtext)


tab1a <- read.csv("/Users/Kat/Dropbox/HSPH_T32/Potential projects/ABCD melatonin/results/mel_use_tab1a_imp_update_Jan16_simple.csv") %>%
  filter(Short_Label != ("Race/ethnicity: *"))


p <- ggplot(tab1a, aes(x = OR, y = factor(Short_Label, levels = c('Nationally-normed health and environment domain score (COI), per SD',
                                                                  'Parent behavior inventory: acceptance subscale score, per SD',
                                                                  'Child attention problems raw score (CBCL), per SD',
                                                                  'Parent somatic complains (ASR), per SD *',
                                                                  'Index of Concentration at the Extremes (Income + Race) *',
                                                                  'Child thought problems raw score (CBCL), per SD',
                                                                  'Adult violent crimes in the county, per SD *',
                                                                  'Other vs non-Hispanic white',
                                                                  'Non-Hispanic Black vs non-Hispanic white',
                                                                  'Asian vs non-Hispanic white',
                                                                  'Hispanic vs non-Hispanic white',
                                                                  'Parent been to a doctor or counselor due to emotional/mental problem *',
                                                                  'Problems with initiating and maintaining sleep score, per SD *')))) + 
    geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") + 
    geom_errorbarh(aes(xmax = highCI, xmin = lowCI), size = 1, height = 
                     .5, color = "gray50") +
    geom_point(size = 5, color = "orange") +
    scale_x_continuous(breaks = seq(0, 2.6, 0.2), labels = seq(0, 2.6, 0.2)) +
    theme_bw()+
    theme(axis.text=element_text(size=16),
          panel.grid.minor = element_blank()) +
    ylab("") +
    xlab("Odds ratio") +
    annotate(geom = "text", y =1, x = 2, 
             label = "AUC=0.735", size = 5, hjust = 0) + 
    ggtitle("<span style='font-size: 22pt;'>Odds ratios for the top 10 most important predictors of melatonin use</font>")+
  theme(plot.title = element_markdown())
p
