# Check the conservation of the positive selected genes occured in the branch site model related to the crabtree effect and heat-tolerance

# 2020.6.4
# Hongzhong Lu

#load package
#library(tidyverse)
library(ggplot2)
library(dplyr)
library(stringr)
library(readr)
library(readxl)
library(hongR)


# input more parameters related to ortholog
ortholog_occurence <- read_tsv("data/ortholog_occurence_num_all.tsv")
# update the result
# only choose the positively selected genes related to heat tolerance?
trait_heat_result <- read_csv("/Users/luho/Documents/branch_site_heat/heat_result_all_update2.csv")
# filter based on the select species and select clade
interest_OG <- filter(trait_heat_result, Select_species>=10 & select_clade_all >=3)
# input the data from second calculation
trait_heat_result2 <- read_csv("/Users/luho/Documents/branch_site_heat_2nd/heat_result_all_update2.csv")
trait_heat_result2$Select_species <- as.numeric(trait_heat_result2$Select_species)
# filter based on the select species and select clade
interest_OG2 <- trait_heat_result2[trait_heat_result2$Select_species >= 10 & trait_heat_result2$select_clade_all >=3,]
interest_OG2 <- interest_OG2[!is.na(interest_OG2$OG),]
# calculate the common
common_OGs <- intersect(interest_OG$OG, interest_OG2$OG)
# OG for heat-tolerance
OG_input_heat <- trait_heat_result[,c("OG")]
OG_input_heat$type <- NA
OG_input_heat$type[OG_input_heat$OG %in% common_OGs] <- "Top selected"
OG_input_heat$type[!(OG_input_heat$OG %in% common_OGs)] <- "Reference"
OG_input_heat$species_num <- getSingleReactionFormula(ortholog_occurence$species_num,ortholog_occurence$ID,OG_input_heat$OG)
OG_input_heat$species_num <- as.numeric(OG_input_heat$species_num)

# plot
ggplot(OG_input_heat, aes(x=type, y=species_num, fill=type)) + 
  stat_boxplot(geom ='errorbar', width = 0.25) + # add caps
  geom_boxplot() +
  xlab('Classification') + ylab('Species number') +
  #theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=20, family="Arial"),
        axis.title=element_text(size=24,family="Arial"),
        legend.text = element_text(size=10, family="Arial")) +
  ggtitle('') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) #+
# output size 5 x 5


# wilcon.test
g1 <- OG_input_heat$species_num[OG_input_heat$type=="Top selected"]
g2 <- OG_input_heat$species_num[OG_input_heat$type=="Reference"]
wilcox.test(g1,g2, alternative = "two.sided")
t.test(g1,g2)
