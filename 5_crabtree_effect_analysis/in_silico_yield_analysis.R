# These scripts are used to produce PCA plot for yeast species with different traits
# it seems difficult to classify the crabree positive species with the negative species based on the yield data
# but the yield of ATP is significantly low in crabree positive speceis.
# 2020.9.16
# Hongzhong Lu

#load package
#library(tidyverse)
library(ggplot2)
library(dplyr)
library(stringr)
library(readr)
library(readxl)
library(hongR)
library(ggfortify)

yield_prediction <- read_excel("data/yield_prediction_with_332_yeast_GEMs.xlsx")

# change the data into a dataframe
yeast_species_classification <- read_csv("data/yeast_species_classification.csv")
yeast_species_refine <- yeast_species_classification[!is.na(yeast_species_classification$crabtree_effect),]
yeast_species_refine <- yeast_species_refine[, c("old_species_id","crabtree_effect")]

# only keep the crabtree positive and negative species
yield_prediction_filter0 <- yield_prediction[yield_prediction$species %in% yeast_species_refine$old_species_id,]
yield_prediction_filter0$trait <- getSingleReactionFormula(yeast_species_refine$crabtree_effect, yeast_species_refine$old_species_id, yield_prediction_filter0$species)

# PCA plot again
df2 <- yield_prediction_filter0[, c(2:48)]

autoplot(prcomp(df2), data = yield_prediction_filter0, colour = 'trait') +
  theme(axis.text=element_text(size=20, family="Arial"),
        axis.title=element_text(size=20, family="Arial") ) +
  ggtitle('') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) +
  theme(legend.text=element_text(size=15),
        legend.title =element_text(size=15))


# specially using ATP
ggplot(yield_prediction_filter0, aes(x = trait, y = ATP, fill = trait)) + geom_boxplot() +
  xlab("") + ylab("ATP yield") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.position = "none") +
  theme(
    axis.text = element_text(size = 16, family = "Arial"),
    axis.title = element_text(size = 16, family = "Arial"),
    legend.text = element_text(size = 10, family = "Arial")
  ) +
  ggtitle("") +
  theme(panel.background = element_rect(fill = "white", color = "black", size = 1))

# statistical analysis
# p-value = 1.991e-10
t.test(yield_prediction_filter0$ATP[yield_prediction_filter0$trait=="Yes"], yield_prediction_filter0$ATP[yield_prediction_filter0$trait=="No"])


