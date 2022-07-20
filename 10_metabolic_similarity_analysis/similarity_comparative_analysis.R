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
evolution_distance <- read.table('data/taxa_pairwise_dist.txt')


# input the the model similarity
similarity_file <- 'data_for_GEMs_similarity/model_similirity_from_biocyc.txt'





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
s_kegg <- similarity[,c("similirity", "source")]



######## biocyc
similarity_file <- 'data_for_GEMs_similarity/model_similirity_from_biocyc.txt'
similarity <- read.table(similarity_file, header = TRUE, stringsAsFactors = FALSE)
unique(similarity$s1)


species_id <- read_excel("data_for_tree/343taxa_speicies-name_clade-name_color-code.xlsx")

similarity$s1 <- str_replace_all(similarity$s1, ".txt", "")
similarity$s2 <- str_replace_all(similarity$s2, ".txt", "")


similarity$s1[similarity$s1=="metschnikowia_matae_maris"] <- "Metschnikowia_matae_maris"
similarity$s2[similarity$s2=="metschnikowia_matae_maris"] <- "Metschnikowia_matae_maris"

similarity$s1_update <- getSingleReactionFormula(species_id$speceis_names_fig2,species_id$old_speceis_names,similarity$s1)
similarity$s2_update <- getSingleReactionFormula(species_id$speceis_names_fig2,species_id$old_speceis_names,similarity$s2)

similarity$source <- "RAVEN_biocyc"
s_biocyc <- similarity[,c("similirity", "source")]



######## manual
similarity_file <- 'data_for_GEMs_similarity/manual_model_similirity.txt'
similarity <- read.table(similarity_file, header = TRUE, stringsAsFactors = FALSE)
unique(similarity$s1)


species_id <- read_excel("data_for_tree/343taxa_speicies-name_clade-name_color-code.xlsx")

similarity$s1 <- str_replace_all(similarity$s1, ".txt", "")
similarity$s2 <- str_replace_all(similarity$s2, ".txt", "")


similarity$s1[similarity$s1=="metschnikowia_matae_maris"] <- "Metschnikowia_matae_maris"
similarity$s2[similarity$s2=="metschnikowia_matae_maris"] <- "Metschnikowia_matae_maris"

similarity$s1_update <- getSingleReactionFormula(species_id$speceis_names_fig2,species_id$old_speceis_names,similarity$s1)
similarity$s2_update <- getSingleReactionFormula(species_id$speceis_names_fig2,species_id$old_speceis_names,similarity$s2)

similarity$source <- "semi_auto_GEMs"
s_manual <- similarity[,c("similirity", "source")]





######### combine dataframe
df3 <- rbind(s_kegg, s_biocyc, s_manual)
# plot density plot
ggplot(df3 , aes(x = similirity, fill = source)) + geom_density(alpha = 0.5) +
  xlim(0.4, 1) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.text=element_text(size=20, family="Arial"),
        axis.title=element_text(size=24, family="Arial")) +
  theme(legend.title = element_text(size=14), legend.text = element_text(size=14)) +
  theme(legend.position = c(0.25, 0.8)) +
  xlab("Similarity") + ylab("Density") +
  theme(legend.text=element_text(size=14, family="Arial"))
ggsave(out <- paste('result/','model_similarity_comparison_from_different_source','.eps', sep = ""), width=5, height=4, dpi=600)



