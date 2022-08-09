# Similarity analysis based on model structure, trait and evolution distance
library(ape)
library(hongR)
library(ggplot2)
library(readxl)
library(stringr)
library(readr)
library(Biostrings)
library(tidyverse)


# check the model similarity and genomics evolutionary distance

# input the evolutionary distance
similarity_file <- 'data_for_GEMs_similarity/model_similirity_from_kegg.txt'
similarity <- read.table(similarity_file, header = TRUE, stringsAsFactors = FALSE)
unique(similarity$s1)


species_id <- read_excel("data_for_tree/343taxa_speicies-name_clade-name_color-code.xlsx")

similarity$s1 <- str_replace_all(similarity$s1, ".txt", "")
similarity$s2 <- str_replace_all(similarity$s2, ".txt", "")


similarity$s1[similarity$s1=="metschnikowia_matae_maris"] <- "Metschnikowia_matae_maris"
similarity$s2[similarity$s2=="metschnikowia_matae_maris"] <- "Metschnikowia_matae_maris"

similarity$s1_update <- getSingleReactionFormula(species_id$speceis_names_fig2,species_id$old_speceis_names,similarity$s1)
similarity$s2_update <- getSingleReactionFormula(species_id$speceis_names_fig2,species_id$old_speceis_names,similarity$s2)
similarity$source <- "RAVEN_kegg"

# classify the clade information
similarity$clade1 <- getSingleReactionFormula(species_id$`Major clade`,species_id$speceis_names_fig2,similarity$s1_update)
similarity$clade2 <- getSingleReactionFormula(species_id$`Major clade`,species_id$speceis_names_fig2,similarity$s2_update)

similarity00 <- similarity[similarity$clade1==similarity$clade2,]

similarity01 <- similarity00[, c("similirity","clade1")]

similarity02 <-  similarity[similarity$clade1 != similarity$clade2,]
similarity02 <- similarity02[, c("similirity","clade1")]
similarity02$clade1 <- "a_inter_clade"

similarity_combine <- rbind(similarity01, similarity02)
# plot density plot

# plot
ggplot(similarity_combine ,aes(x=clade1, y=similirity)) + 
  stat_boxplot(geom ='errorbar', width = 0.1) + # add caps
  geom_boxplot(fill = "#FF6666") +
  xlab('') + ylab('Similarity') +
  #theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=12, family="Arial"),
        axis.title=element_text(size=16,family="Arial"),
        legend.text = element_text(size=10, family="Arial")) +
  ggtitle('') +
  geom_hline(yintercept=0.845, linetype="dashed", 
             color = "blue", size=1) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) # +
  #coord_flip()


##################### statistical analysis
m1 <- similarity_combine$similirity[similarity_combine$clade1=="a_inter_clade"]
m2 <- similarity_combine$similirity[similarity_combine$clade1=="CUG-Ala"]
t.test(m1, m2)
wilcox.test(m1,m2, alternative = "two.sided")
##################### statistical analysis

