# These scripts are used to produce map related to evolution project
# 2020.4.15
# Hongzhong Lu

#load package
#library(tidyverse)
library(ggplot2)
library(dplyr)
library(stringr)
library(readr)
library(ggplot2)
library(readxl)
library(hongR)

# input the data
yeast_species_classification <- read_csv("data/yeast_species_classification.csv")

yeast_fermentation_data <- read_excel("data/Yeast_fermentation_data.xlsx")


yeast_fermentation_data$crabtree <- getSingleReactionFormula(yeast_species_classification$crabtree_effect,yeast_species_classification$old_species_id,yeast_fermentation_data$Model_name)
# if the ethanol production is zero, then the related species will be regarded as carabtree positive
yeast_fermentation_data$crabtree[yeast_fermentation_data$`Yield: EtOH/glc(g/g)` <=0.001] <- "No"

# some newly added
yeast_fermentation_data$crabtree[yeast_fermentation_data$Strain =="K. phaffii - control\r\n"] <- "No"
yeast_fermentation_data$crabtree[74:77] <- "Yes"


yeast_fermentation_filter <- yeast_fermentation_data %>% filter(.,crabtree !="NA")
yeast_fermentation_filter$crabtree <- as.factor(yeast_fermentation_filter$crabtree)
# remove the duplicated species
yeast_fermentation_filter <- yeast_fermentation_filter[!duplicated(yeast_fermentation_filter$Model_name),]
yeast_fermentation_filter$`Yield: Biomass/Glc`


# statistical analysis
g1 <- yeast_fermentation_filter$`Yield: Biomass/Glc`[yeast_fermentation_filter$crabtree=="Yes"]
g2 <- yeast_fermentation_filter$`Yield: Biomass/Glc`[yeast_fermentation_filter$crabtree=="No"]
t.test(g1, g2) 
wilcox.test(g1,g2, alternative = "two.sided")

m1 <- yeast_fermentation_filter$`Growth rate*: (1/h)`[yeast_fermentation_filter$crabtree=="Yes"]
m2 <- yeast_fermentation_filter$`Growth rate*: (1/h)`[yeast_fermentation_filter$crabtree=="No"]
t.test(m1, m2) # gowth p value 0.054
wilcox.test(m1,m2, alternative = "two.sided")

n1 <- yeast_fermentation_filter$`Glc/Biomass/h((mmol/gDW,h))`[yeast_fermentation_filter$crabtree=="Yes"]
n2 <- yeast_fermentation_filter$`Glc/Biomass/h((mmol/gDW,h))`[yeast_fermentation_filter$crabtree=="No"]
t.test(n1, n2) 
wilcox.test(n1,n2, alternative = "two.sided")

e1 <- yeast_fermentation_filter$`EtOH/Biomass/h((mmol/gDW,h))`[yeast_fermentation_filter$crabtree=="Yes"]
e2 <- yeast_fermentation_filter$`EtOH/Biomass/h((mmol/gDW,h))`[yeast_fermentation_filter$crabtree=="No"]
t.test(e1,e2) 
wilcox.test(e1,e2, alternative = "two.sided")






# plot
yeast_fermentation_filter %>%
  ggplot(aes(x=crabtree,y=`Growth rate*: (1/h)`, fill=crabtree)) +
  stat_boxplot(geom ='errorbar', width = 0.25) + # add caps
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = c(0.85, 0.2)) +
  theme(axis.text=element_text(size=20, family="Arial"),
        axis.title=element_text(size=24,family="Arial"),
        legend.text = element_text(size=10, family="Arial")) +
  ggtitle('') +
  theme(legend.position = "none") +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) +
  xlab("Crabtree") + ylab("Growth rate (1/h)")
# output size 5 x 5 



yeast_fermentation_filter %>%
  ggplot(aes(x=crabtree,y=`Glc/Biomass/h((mmol/gDW,h))`, fill=crabtree)) +
  stat_boxplot(geom ='errorbar', width = 0.25) + # add caps
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = c(0.85, 0.2)) +
  theme(axis.text=element_text(size=20, family="Arial"),
        axis.title=element_text(size=24,family="Arial"),
        legend.text = element_text(size=10, family="Arial")) +
  ggtitle('') +
  theme(legend.position = "none") +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) +
  xlab("Crabtree") + ylab("qS (mM/gDW.h)")



yeast_fermentation_filter %>%
  ggplot(aes(x=crabtree,y=`Yield: EtOH/glc(g/g)`, fill=crabtree)) +
  stat_boxplot(geom ='errorbar', width = 0.25) + # add caps
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = c(0.85, 0.2)) +
  theme(axis.text=element_text(size=20, family="Arial"),
        axis.title=element_text(size=24,family="Arial"),
        legend.text = element_text(size=10, family="Arial")) +
  ggtitle('') +
  theme(legend.position = "none") +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) +
  xlab("Crabtree") + ylab("Yethanol/S (g/g)")


yeast_fermentation_filter %>%
  ggplot(aes(x=crabtree,y=`Yield: Biomass/Glc`, fill=crabtree)) +
  stat_boxplot(geom ='errorbar', width = 0.25) + # add caps
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = c(0.85, 0.2)) +
  theme(axis.text=element_text(size=20, family="Arial"),
        axis.title=element_text(size=24,family="Arial"),
        legend.text = element_text(size=18, family="Arial")) +
  ggtitle('') +
  theme(legend.position = "none") +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) +
  xlab("Crabtree") + ylab("YX/S(C mol/C mol)")












# not used sub figure
yeast_fermentation_filter %>%
  ggplot(aes(x=crabtree,y=`EtOH/Biomass/h((mmol/gDW,h))`, fill=crabtree)) +
  stat_boxplot(geom ='errorbar', width = 0.25) + # add caps
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = c(0.85, 0.2)) +
  theme(axis.text=element_text(size=20, family="Arial"),
        axis.title=element_text(size=24,family="Arial"),
        legend.text = element_text(size=10, family="Arial")) +
  ggtitle('') +
  theme(legend.position = "none") +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) +
  xlab("Crabtree") + ylab("qEthanol (mM/gDW.h)")
