# conservation score analysis for machine learning

library(ggplot2)
library(readr)
library(dplyr)
library(hongR)
library(stringr)


# input the conservation score data
conservation_score <- read_csv("data/conservation_score_sce_for_ML.csv")
conservation_score <- filter(conservation_score, !is.na(conservation_score))
conservation_score$OG <- str_replace_all(conservation_score$OG, "_aa_aligned.fasta", "")
ortholog <- read_tsv("data/ortholog_occurence_num_all.tsv")
conservation_score$species_num <- getSingleReactionFormula(ortholog$species_num, ortholog$ID, conservation_score$OG)
conservation_score$species_num <- as.numeric(conservation_score$species_num)
conservation_score_filter <- conservation_score[conservation_score$species_num > 2, ]

ggplot(conservation_score_filter, aes(conservation_score)) +
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "blue", bins = 50) +
  xlim(0, 1) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  geom_density(colour = "black", alpha = 0.3, fill = "grey50") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20, face = "bold")) +
  xlab("average_conservation_score")

write.table(conservation_score_filter, "data/conservation_score_sce_for_ML_filter_essential_gene.txt")




conservation_score <- read_csv("data/conservation_score_complex.csv")
conservation_score <- filter(conservation_score, !is.na(conservation_score))
conservation_score$OG <- str_replace_all(conservation_score$OG, "_aa_aligned.fasta", "")
ortholog <- read_tsv("data/ortholog_occurence_num_all.tsv")
conservation_score$species_num <- getSingleReactionFormula(ortholog$species_num, ortholog$ID, conservation_score$OG)
conservation_score$species_num <- as.numeric(conservation_score$species_num)
conservation_score_filter <- conservation_score[conservation_score$species_num > 2, ]

ggplot(conservation_score_filter, aes(conservation_score)) +
        geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "blue", bins = 50) +
        xlim(0, 1) +
        theme(panel.background = element_rect(fill = "white", colour = "black")) +
        geom_density(colour = "black", alpha = 0.3, fill = "grey50") +
        theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20, face = "bold")) +
        xlab("average_conservation_score")

write.table(conservation_score_filter, "data/conservation_score_complex_filter.txt")

# save other information for ML
#ortholog_complex <- read_csv("data/ortholog_complex.csv")
#ortholog_filter_for_complex <- ortholog[ortholog$ID %in% ortholog_complex$OG_ID, ]
#write.table(ortholog_filter_for_complex, "data/conservation_score_complex_filter.txt")




