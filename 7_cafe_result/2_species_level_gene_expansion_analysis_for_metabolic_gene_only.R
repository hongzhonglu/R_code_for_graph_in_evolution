# note------------
# here we try to analyze the functions of expanded or contracted genes
# The expansion or contraction of gene families along a specific lineage can be due to chance or can be the result of natural selection.
# 2020-09-08
# Hongzhong Lu
# No interesting results could be found from expansion gene analysis of crabtree positive species

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
gene_family3 <- fromJSON(file="cafe_data/CAFE_all_species_analysis.json")
# input the metabolic gene from draft yeast GEMs by RAVEN2
ortholog_id_biocyc <- read_csv("cafe_data/ortholog_id_biocyc.csv")
ortholog_id_kegg <- read_csv("cafe_data/ortholog_id_kegg.csv")
length(union(ortholog_id_biocyc$ID, ortholog_id_kegg$ID))



# change it into a dataframe
species_all <- vector()
expanded_all <- vector()
extracted_all <- vector()
for (i in 1:length(gene_family3)){
  print(gene_family3[[i]]$organism)
  species_all <- c(species_all, gene_family3[[i]]$organism)
  expanded_all <- c(expanded_all, length(which(gene_family3[[i]]$expansion %in% ortholog_id_biocyc$ID)))
  extracted_all <- c(extracted_all, length(which(gene_family3[[i]]$contraction %in% ortholog_id_biocyc$ID)))
}

species_gene_families <- data.frame(species=species_all, expanded=expanded_all, extracted = extracted_all, stringsAsFactors = FALSE)
species_gene_families$ratio_of_expanded_to_extracted <- species_gene_families$expanded/species_gene_families$extracted
# save the data, it is possilbe to further analyze the relation between the gene expansion and proteome size and reaction size.
write.table(species_gene_families, "result/M_gene_family_expansion_extraction_for_332_species.txt", row.names = FALSE, sep = "\t")






# plot the result
plot(density(species_gene_families$expanded))
plot(density(species_gene_families$extracted))

B1 <- select(species_gene_families, species, expanded)
B2 <- select(species_gene_families, species, extracted)
colnames(B1) <- c("species","Num")
colnames(B2) <- c("species","Num")
B1$group <- "expanded"
B2$group <- "contracted"
merge_df <- rbind.data.frame(B1, B2)

ggplot(merge_df, aes(x = Num, fill = group)) + geom_histogram(alpha = 0.5, bins = 100) +
  xlim(0, 500) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = c(0.8, 0.8))


ggplot(merge_df, aes(x = Num, fill = group)) + 
  geom_density(alpha = 0.5, bins = 10) +
  xlim(0, 700) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 24)) +
  theme(legend.position = c(0.7, 0.8),legend.title=element_text(size=20), 
        legend.text=element_text(size=20))
#output 8 x 6

t.test(species_gene_families$expanded, species_gene_families$extracted)
wilcox.test(species_gene_families$expanded, species_gene_families$extracted, alternative = "two.sided")
