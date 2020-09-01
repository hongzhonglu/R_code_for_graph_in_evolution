# note------------
# here we try to analyze the functions of expanded or contracted genes
# The expansion or contraction of gene families along a specific lineage can be due to chance or can be the result of natural selection.
# 2020-06-02
# Hongzhong Lu


library(ggplot2)
library(readr)
library(dplyr)
library(stringr)
library(hongR)
library(readxl)
library(rjson)


# now the result only contains the FUBAR analysis result. In the future, the fel analysis results will be included.
# classification of OGs with selected sites based on KO or function annotation.
og_panID_mapping <- read_tsv("data/representatives.tsv")

# get the KO id
panID_KO <- read.table("kegg/ko_annotation_of_all_OG.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
panID_KO0 <- left_join(og_panID_mapping, panID_KO, by = c( "representative" = "query") )


# get the pathway
KO_pathway <- read.table("kegg/ko_pathway.txt", sep = "\t", stringsAsFactors = FALSE)
pathway <- read.table("kegg/pathway_list_kegg.txt", sep = "\t", stringsAsFactors = FALSE)
KO_pathway <- filter(KO_pathway,str_detect(KO_pathway$V1, "map"))
KO_pathway$V2 <- str_replace_all(KO_pathway$V2,"ko:","")
KO_pathway0 <- left_join(KO_pathway, pathway, by=c("V1"="V1"))

# third comparsion
# expanded gene in species levels
gene_family3 <- fromJSON(file="cafe_data/expansion_all_species_analysis.json")
# change it into a dataframe
species_all <- vector()
expanded_all <- vector()
extracted_all <- vector()
for (i in 1:length(gene_family3)){
  print(gene_family3[[i]]$organism)
  species_all <- c(species_all, gene_family3[[i]]$organism)
  expanded_all <- c(expanded_all, length(gene_family3[[i]]$expansion))
  extracted_all <- c(extracted_all, length(gene_family3[[i]]$contraction))
}

species_gene_families <- data.frame(species=species_all, expanded=expanded_all, extracted = extracted_all, stringsAsFactors = FALSE)
species_gene_families$ratio_of_expanded_to_extracted <- species_gene_families$expanded/species_gene_families$extracted
# save the data, it is possilbe to further analyze the relation between the gene expansion and proteome size and reaction size.
write.table(species_gene_families, "result/gene_family_expansion_extraction_for_332_species.txt", row.names = FALSE, sep = "\t")
plot(density(species_gene_families$expanded))


# change the complex list data into simple list data 
species_num <- length(gene_family3)
species_expanded <- list()
species_extracted <- list()
for (i in 1: species_num){
  print(i)
  species_name <- gene_family3[[i]]$organism
  species_expanded[[species_name]] <- gene_family3[[i]]$expansion
  species_extracted[[species_name]] <- gene_family3[[i]]$contraction
}
# then we next extract the species with crabtree effect
yeast_species_classification <- read_csv("data/yeast_species_classification.csv")
yeast_species_classification <- yeast_species_classification[!is.na(yeast_species_classification$crabtree_effect),]
species_crabtree <- yeast_species_classification$old_species_id[yeast_species_classification$crabtree_effect=="Yes"]
species_no_crabtree <- yeast_species_classification$old_species_id[yeast_species_classification$crabtree_effect=="No"]

# put the species from crabtree and non crabtree together
og_crabtree <- vector()
for (i in species_crabtree){
  print(i)
  og_choose <- species_expanded[[i]]
  if(length(og_choose) >=1){
    og_crabtree <- c(og_crabtree, og_choose)
    
  }
}

og_crabtree_unique <- unique(og_crabtree)
og_crabtree_df <- data.frame(og=og_crabtree_unique, stringsAsFactors = FALSE)
# then check whether og exist in each species with crabtree effect

species_crabtree <- species_crabtree[1:(length(species_crabtree)-1)]
for (i in species_crabtree){
  print(i)
  #i <- "Saccharomyces_cerevisiae"
  og_choose <- species_expanded[[i]]
  og_crabtree_df[[i]] <- NA
  og_crabtree_df[[i]] <- as.numeric(og_crabtree_df$og %in%  og_choose )
}
og_crabtree_df$species_count <- apply(og_crabtree_df[, c(2:26)],1,sum)
og_crabtree_df_filter1 <- og_crabtree_df[og_crabtree_df$species_count>=25,]
og_crabtree_df_filter1 <- og_crabtree_df_filter1[,c("og", "species_count")]



# non no_crabtree together
og_no_crabtree <- vector()
for (i in species_no_crabtree){
  print(i)
  og_choose <- species_expanded[[i]]
  if(length(og_choose) >=1){
    og_no_crabtree <- c(og_no_crabtree, og_choose)
    
  }
}

og_no_crabtree_unique <- unique(og_no_crabtree)
og_no_crabtree_df <- data.frame(og=og_no_crabtree_unique, stringsAsFactors = FALSE)
# then check whether og exist in each species with no_crabtree effect
species_no_crabtree <- species_no_crabtree[1:(length(species_no_crabtree)-1)]
for (i in species_no_crabtree){
  print(i)
  #i <- "Saccharomyces_cerevisiae"
  og_choose <- species_expanded[[i]]
  og_no_crabtree_df[[i]] <- NA
  og_no_crabtree_df[[i]] <- as.numeric(og_no_crabtree_df$og %in%  og_choose )
}
og_no_crabtree_df$species_count <- apply(og_no_crabtree_df[, c(2:76)],1,sum)
plot(density(og_no_crabtree_df$species_count))
og_no_crabtree_df_filter1 <- og_no_crabtree_df[og_no_crabtree_df$species_count>=60,]
og_no_crabtree_df_filter1 <- og_no_crabtree_df_filter1[,c("og", "species_count")]



# only choose expanded OGs exist in crabtree effective species but not in crabtree nagative species
og_crabtree_df_filter2 <- og_crabtree_df_filter1[!(og_crabtree_df_filter1$og %in% og_no_crabtree_df_filter1$og), ]

og_crabtree_df_filter2$KO <- getSingleReactionFormula(panID_KO0$ko,panID_KO0$ortho_id, og_crabtree_df_filter2$og)
# be careful about it as one KO can belong to multiple pathways
og_crabtree_df_filter2$pathway <- getSingleReactionFormula(KO_pathway0$V2.y,KO_pathway0$V2.x, og_crabtree_df_filter2$KO)
og_crabtree_df_filter2$panID <- getSingleReactionFormula(panID_KO0$representative,panID_KO0$ortho_id, og_crabtree_df_filter2$og)
#extracted sce gene
sce_gene <- og_crabtree_df_filter2$panID[str_detect(og_crabtree_df_filter2$panID,"Saccharomyces_cerevisiae")]
sce_gene <- str_replace(sce_gene,"Saccharomyces_cerevisiae@","")
print(paste0(sce_gene, collapse=","))