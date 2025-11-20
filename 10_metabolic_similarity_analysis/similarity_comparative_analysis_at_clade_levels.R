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




# add other information
# check the model similarity and trait similarity
# input the the model similarity

# input the evolutionary distance
trait_similarity <- read.table('data_for_GEMs_similarity/trait_similirity.txt', header = TRUE, stringsAsFactors = FALSE)

library(reshape2)
trait_similarity0 <- acast(trait_similarity, s1~s2, value.var="similirity")


# filter
strain_set <- colnames(trait_similarity0)
similarity1 <- similarity[similarity$s1 %in% strain_set,]
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


# check 
"Lipomycetaceae"
df_check <- similarity2[similarity2$clade1=="Lipomycetaceae", ]
df_check <- df_check[df_check$clade2=="Lipomycetaceae", ]
plot(df_check$similirity, df_check$trait_similarity)



# classify the clade information
similarity2$clade1 <- getSingleReactionFormula(species_id$`Major clade`,species_id$speceis_names_fig2,similarity2$s1_update)
similarity2$clade2 <- getSingleReactionFormula(species_id$`Major clade`,species_id$speceis_names_fig2,similarity2$s2_update)

similarity00 <- similarity2[similarity2$clade1==similarity2$clade2,]

similarity01 <- similarity00[, c("similirity","trait_similarity","clade1")]

similarity02 <-  similarity2[similarity2$clade1 != similarity2$clade2,]
similarity02 <- similarity02[, c("similirity","trait_similarity","clade1")]
similarity02$clade1 <- "a_inter_clade"

similarity_combine <- rbind(similarity01, similarity02)





# plot density plot

# plot
similarity_combine = similarity_combine[similarity_combine$clade1 != "Sporopachydermia clade", ]


ggplot(similarity_combine ,aes(x=clade1, y=trait_similarity)) + 
  stat_boxplot(geom ='errorbar', width = 0.1) + # add caps
  geom_boxplot(fill = "#FF6666") +
  xlab('') + ylab('Trait similarity') +
  #theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=12, family="Arial"),
        axis.title=element_text(size=16,family="Arial"),
        legend.text = element_text(size=10, family="Arial")) +
  ggtitle('') +
  geom_hline(yintercept=0.56, linetype="dashed", 
             color = "blue", size=1) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) # +



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
  geom_hline(yintercept=0.85, linetype="dashed", 
             color = "blue", size=1) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) # +



unique_clade <- unique(similarity_combine$clade1)
unique_clade  <- unique_clade[unique_clade != "Sporopachydermia clade"]
cor_all <- vector()
for (xx in unique_clade){
  #xx <- "Lipomycetaceae"
  df <- similarity_combine[similarity_combine$clade1==xx,]
  ss <- cor.test(df$similirity, df$trait_similarity)
  yy <- unlist(ss['estimate'])
  cor_all <- c(cor_all,yy)
  
}

new_result <- data.frame(clade=unique_clade, cor=cor_all, stringsAsFactors = FALSE)

ggplot(data=new_result, aes(x=clade, y=cor)) +
  geom_bar(stat="identity", colour="black", fill="white",width = 0.6) +
  xlab('') + ylab('Pearson coefficient') +
  #theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=12, family="Arial"),
        axis.title=element_text(size=16,family="Arial"),
        legend.text = element_text(size=10, family="Arial")) +
  ggtitle('') +
  geom_hline(yintercept=0.18, linetype="dashed", 
             color = "blue", size=1) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) # +



##################### statistical analysis
m1 <- similarity_combine$similirity[similarity_combine$clade1=="a_inter_clade"]
m2 <- similarity_combine$similirity[similarity_combine$clade1=="CUG-Ala"]
t.test(m1, m2)
wilcox.test(m1,m2, alternative = "two.sided")
##################### statistical analysis




