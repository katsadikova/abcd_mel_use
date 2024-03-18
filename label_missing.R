library(tidyverse)

misses <- read.csv("/Users/Kat/Dropbox/HSPH_T32/Potential projects/ABCD melatonin/temp_data/misses.csv")
labels <- read.csv("/Users/Kat/Dropbox/HSPH_T32/Potential projects/ABCD melatonin/temp_data/Predictors_to_merge.csv")

names(misses) <- c("Variable","nmiss","n","propmiss")
names(labels)

misses_labeled <- merge(misses,labels,by="Variable") %>%
  arrange(Domain,Short_Domain)

write.csv(misses_labeled, file="/Users/Kat/Dropbox/HSPH_T32/Potential projects/ABCD melatonin/temp_data/misses_labeled.csv")
