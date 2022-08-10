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




# loop
similarity$distance <- NA
for (i in 1:nrow(similarity)){
  print(i)
  row0 <- similarity[i,]
  x1 <- row0$s1_update[1]
  x2 <- row0$s2_update[1]
  # select the distance value based on x1 x2
  evolution_distance0 <- evolution_distance[x1, x2]
  print(evolution_distance0)
  if(length(evolution_distance0) > 0) {
    similarity$distance[i] <- evolution_distance0
  } 
  else{
    similarity$distance[i] < NA
    
  }
}

similarity1 <- similarity[!is.na(similarity$distance),]
unique(similarity1$s2)



setdiff(unique(similarity$s1_update), unique(similarity1$s1_update))
setdiff(unique(similarity$s2_update), unique(similarity1$s2_update))
similarity_check <- similarity[similarity$s1_update=="NA",]



# check: Alloascoidea_hylecoeti   metschnikowia_matae_maris

ggplot(similarity1, aes(x=distance, y=similirity) ) +
  geom_hex(bins = 120) +
  scale_fill_continuous(type = "viridis") +
  geom_smooth(method = "gam") +
  theme_bw() +
  xlab("Evolutionary distance") +
  ylab("Model similarity") +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1),
      plot.margin = margin(1, 1, 1, 1, "cm")) +
  theme(axis.text=element_text(size=12, family="Arial"), axis.title=element_text(size=16, family="Arial")) +
  theme(legend.title = element_text(size=16), legend.text = element_text(size=16))
ggsave(out <- paste('result/','manual_model_similarity_evolution_distance','.eps', sep = ""), width=5, height=4, dpi=600)

cor.test(similarity1$distance, similarity1$similirity)





# check the model similarity and trait similarity
# input the the model similarity

# input the evolutionary distance
trait_similarity <- read.table('data_for_GEMs_similarity/trait_similirity.txt', header = TRUE, stringsAsFactors = FALSE)

library(reshape2)
trait_similarity0 <- acast(trait_similarity, s1~s2, value.var="similirity")


# filter
strain_set <- colnames(trait_similarity0)
similarity1 <- similarity1[similarity1$s1 %in% strain_set,]
similarity1 <- similarity1[similarity1$s2 %in% strain_set,]


# loop
similarity1$trait_similarity <- NA
for (i in 1:nrow(similarity1)){
  print(i)

  row0 <- similarity1[i,]
  x1 <- row0$s1[1]
  x2 <- row0$s2[1]

  
  # Note that print(b) fails since b doesn't exist
  skip_to_next <- FALSE
  tryCatch(evolution_distance0 <- trait_similarity0[x1, x2], error = function(e) { skip_to_next <<- TRUE})
  
  if(skip_to_next) {
    evolution_distance0 <- trait_similarity0[x2, x1]
    } 
  
  if (is.na(evolution_distance0)){
    evolution_distance0 <- trait_similarity0[x2, x1]
  }

  if(length(evolution_distance0) > 0) {
    similarity1$trait_similarity[i] <- evolution_distance0
  } 
  else{
    similarity1$trait_similarity[i] < NA
    
  }
}

similarity2 <- similarity1[!is.na(similarity1$trait_similarity),]


cor.test(similarity2$similirity, similarity2$trait_similarity)
ggplot(similarity2, aes(x=similirity, y=trait_similarity) ) +
  geom_hex(bins = 120) +
  scale_fill_continuous(type = "viridis") +
  geom_smooth(method = "gam") +
  theme_bw() +
  xlab("Model similarity") +
  ylab("Trait similarity") +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1),
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  theme(axis.text=element_text(size=12, family="Arial"), axis.title=element_text(size=16, family="Arial")) +
  theme(legend.title = element_text(size=16), legend.text = element_text(size=16))
ggsave(out <- paste('result/','manual_model_VS_trait_similarity','.eps', sep = ""), width=5, height=4, dpi=600)



cor.test(similarity2$distance, similarity2$trait_similarity)
ggplot(similarity2, aes(x=distance, y=trait_similarity) ) +
  geom_hex(bins = 120) +
  scale_fill_continuous(type = "viridis") +
  geom_smooth(method = "gam") +
  theme_bw() +
  xlab("Evolutionary distance") +
  ylab("Trait similarity") +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1),
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  theme(axis.text=element_text(size=12, family="Arial"), axis.title=element_text(size=16, family="Arial")) +
  theme(legend.title = element_text(size=16), legend.text = element_text(size=16))
ggsave(out <- paste('result/','evolution_distance_VS_trait_similarity','.eps', sep = ""), width=5, height=4, dpi=600)











# check the genome similarity
genome_similarity <- read.table('data_for_GEMs_similarity/genome_similirity_based_on_ortholog_gene_existence.txt', header = TRUE, stringsAsFactors = FALSE)

genome_similarity$s1 <- str_replace_all(genome_similarity$s1, ".tsv", "")
genome_similarity$s2 <- str_replace_all(genome_similarity$s2, ".tsv", "")

library(reshape2)
genome_similarity0 <- acast(genome_similarity, s1~s2, value.var="similirity")


# filter
strain_set <- colnames(genome_similarity0)
similarity1 <- similarity1[similarity1$s1 %in% strain_set,]
similarity1 <- similarity1[similarity1$s2 %in% strain_set,]


# loop
similarity1$genome_similarity <- NA
for (i in 1:nrow(similarity1)){
  print(i)
  
  row0 <- similarity1[i,]
  x1 <- row0$s1[1]
  x2 <- row0$s2[1]
  
  
  # Note that print(b) fails since b doesn't exist
  skip_to_next <- FALSE
  tryCatch(evolution_distance0 <- genome_similarity0[x1, x2], error = function(e) { skip_to_next <<- TRUE})
  
  if(skip_to_next) {
    evolution_distance0 <- genome_similarity0[x2, x1]
  } 
  
  if (is.na(evolution_distance0)){
    evolution_distance0 <- genome_similarity0[x2, x1]
  }
  
  if(length(evolution_distance0) > 0) {
    similarity1$genome_similarity[i] <- evolution_distance0
  } 
  else{
    similarity1$genome_similarity[i] < NA
    
  }
}







ggplot(similarity1, aes(x=genome_similarity, y=similirity) ) +
  geom_hex(bins = 120) +
  scale_fill_continuous(type = "viridis") +
  geom_smooth(method = "gam") +
  theme_bw() +
  xlab("Genotype similarity") +
  ylab("Model similarity") +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1),
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  theme(axis.text=element_text(size=12, family="Arial"), axis.title=element_text(size=16, family="Arial")) +
  theme(legend.title = element_text(size=16), legend.text = element_text(size=16))
ggsave(out <- paste('result/','Genome_VS_model_similarity_semi_auto_GEMs','.eps', sep = ""), width=5, height=4, dpi=600)


cor.test(similarity1$genome_similarity, similarity1$similirity)



similarity1$EVO_SCORE <- similarity1$genome_similarity / similarity1$distance
similarity1$EVO_SCORE <- similarity1$genome_similarity *(12-similarity1$distance)

ggplot(similarity1, aes(x=EVO_SCORE, y=similirity) ) +
  geom_hex(bins = 120) +
  scale_fill_continuous(type = "viridis") +
  geom_smooth(method = "gam") +
  theme_bw() +
  xlab("EVO_SCORE") +
  ylab("Model similarity") +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1),
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  theme(axis.text=element_text(size=12, family="Arial"), axis.title=element_text(size=16, family="Arial")) +
  theme(legend.title = element_text(size=16), legend.text = element_text(size=16))
ggsave(out <- paste('result/','EVO_SCORE_VS_model_similarity_semi_auto_GEMs','.eps', sep = ""), width=5, height=4, dpi=600)


cor.test(similarity1$EVO_SCORE, similarity1$similirity)


