# These scripts are used to produce map related to evolution project
# 2020.4.15
# Hongzhong Lu

#load package
#library(tidyverse)
library(ggplot2)
library(dplyr)
library(stringr)
library(readr)
library(readxl)
library(hongR)



sce_gene_summary <- read_excel("data/sce_gene_summary.xlsx")
sce_gene_summary$sce_gene <- str_replace_all(sce_gene_summary$sce_gene,"\\[", "") %>% str_replace_all(.,"Saccharomyces_cerevisiae@","") %>%
  str_replace_all(.,"\\]", "") %>% str_replace_all(.,"'", "")

# input more parameters related to ortholog
ortholog_occurence <- read_tsv("data/ortholog_occurence_num_all.tsv")
sce_gene_summary0 <- left_join(sce_gene_summary,ortholog_occurence,by=c("OrthologID"="ID") )




# input the core metabolic gene
core_metabolic_gene <- read_csv("data/core_metabolic_gene.csv")

# input the glucose transporter
glucose_transporter <- "YDL245C or YDL247W or YDR342C or YDR343C or YDR345C or YDR536W or YEL069C or YFL011W or YHR092C or YHR094C or YHR096C or YJL214W or YJL219W or YJR158W or YJR160C or YLR081W or YMR011W or YNR072W or YOL156W or YDR387C"
transporter <- unlist(str_split(glucose_transporter," or "))
interest_gene <- c(core_metabolic_gene$locus, transporter)



# find the OGs for this core metabolic gene
sce_og_mapping <- splitAndCombine(sce_gene_summary$sce_gene, sce_gene_summary$OrthologID,",")
sce_og_mapping$v1 <- str_trim(sce_og_mapping$v1, side = "both") # must note the detail, otherwise the result is wrong
core_sce_og <- sce_og_mapping[sce_og_mapping$v1 %in% interest_gene,]


# then find the species specific genes under each OGs

species_OG_all <- list()
id_mapping_dir <- "/Users/luho/Documents/GitHub/Multi_Scale_Evolution_Analysis/pan_genome/result/id_mapping/"
species_name <- list.files(id_mapping_dir)
for (i in species_name){
  print(i)
  input0 <- paste(id_mapping_dir, i, sep = "")
  species0 <- read_tsv(input0)
  species0$ortholog_id <- str_trim(species0$ortholog_id, side = "both")
  core_sce_og$v2 <- str_trim(core_sce_og$v2, side = "both")
  species0_OG <- species0[species0$ortholog_id %in% core_sce_og$v2,]
  # count the copy
  species0_OG_sum <- as.data.frame(table(species0_OG$ortholog_id), stringsAsFactors = FALSE)
  species_name <- str_replace_all(i, ".tsv", "")
  species_OG_all[[species_name]] <-  species0_OG_sum
  
}


# change the data into a dataframe
yeast_species_classification <- read_csv("data/yeast_species_classification.csv")
yeast_species_refine <- yeast_species_classification[!is.na(yeast_species_classification$crabtree_effect),]
yeast_species_refine <- yeast_species_refine[, c("old_species_id","crabtree_effect")]

df_og <- data.frame(OG= unique(core_sce_og$v2))

for(i in yeast_species_refine$old_species_id){
  print(i)
  OG_inf <- species_OG_all[[i]]
  df_og[, i] <- getSingleReactionFormula(OG_inf$Freq, OG_inf$Var1, df_og$OG)
  df_og[, i] <- as.numeric(df_og[, i])
  df_og[, i][is.na(df_og[, i])] <- 0
}


# refine the data format
rownames(df_og) <- df_og$OG
df_og <- df_og[, 2:103]
df_og0 <- as.data.frame(t(df_og))

df_og1 <- df_og0
species <- rownames(df_og1)
df_og1$species <- NA
df_og1$species <-rownames(df_og1)
df_og1$trait <- getSingleReactionFormula(yeast_species_refine$crabtree_effect, yeast_species_refine$old_species_id, df_og1$species)
df_og1$trait <- str_replace_all(df_og1$trait, "No","Negative")
df_og1$trait <- str_replace_all(df_og1$trait, "Yes","Positive")


df1 <- df_og1[, 1:135]

# PCA plot
library(ggfortify)
autoplot(prcomp(df1), data = df_og1, colour = 'trait') +
  theme(axis.text=element_text(size=20, family="Arial"),
        axis.title=element_text(size=20, family="Arial") ) +
  ggtitle('') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) +
  theme(legend.text=element_text(size=15),
        legend.title =element_text(size=15))
#ggsave(out <- paste('result/','PCA-gene','.eps', sep = ""), width=12, height=8, dpi=300)
s0 <- prcomp(df1)

# pca 3D map
library(pca3d)
df2 <- df1[ , apply(df1, 2, var) != 0]
pca <- prcomp(df2, scale.=TRUE)
pca3d(pca, group=df_og1$trait, legend="bottomleft")
snapshotPCA3d(file="ellipses.png") # save the 3D PCA map


# calculate the average gene copy numbers
# choose genes are significantly different in copies nummber between the crabtree positive and negative species
# re-do the PCA analysis for two types of strains.
df_og_crabtree <- df_og1[df_og1$trait=="Positive",]
df_og_not_crabtree <- df_og1[df_og1$trait=="Negative",]
df_og_compare <- data.frame(OG= unique(core_sce_og$v2))
df_og_compare$ave_crabtree <- NA
df_og_compare$ave_non_crabtree <- NA
df_og_compare$p_value <- NA

for (i in 1:length(df_og_compare$OG)){
  print(i)
  og_choose <- df_og_compare$OG[i]
  crabtree0 <- mean(df_og_crabtree[, c(og_choose)])
  non_crabtree0 <- mean(df_og_not_crabtree[, c(og_choose)])
  ss <- t.test(df_og_crabtree[, c(og_choose)], df_og_not_crabtree[, c(og_choose)])
  p_value0 <- ss[["p.value"]]
  df_og_compare$ave_crabtree[i] <- crabtree0
  df_og_compare$ave_non_crabtree[i] <- non_crabtree0
  df_og_compare$p_value[i] <-  p_value0
  
}


df_og_compare_sinificant <- df_og_compare[df_og_compare$p_value <0.05,]

df2 <- df1[, df_og_compare_sinificant$OG]

# plot again
autoplot(prcomp(df2), data = df_og1, colour = 'trait') +
  theme(axis.text=element_text(size=20, family="Arial"),
        axis.title=element_text(size=20, family="Arial") ) +
  ggtitle('') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) +
  theme(legend.text=element_text(size=10),
        legend.title =element_text(size=10)) +
  theme(legend.position = c(0.85, 0.8))
ggsave(out <- paste('result/','classification_between_crabtree_nagative_and_positive_based_gene_copy_number','.eps', sep = ""), width=5, height=5, dpi=600)


