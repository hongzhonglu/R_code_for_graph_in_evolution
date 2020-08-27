# note------------
# here we try to analyze the functions of genes with positive selected sites based on pathways.
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




# input the data file
# 332 species as a whole
gene_family <- fromJSON(file="cafe_data/CAFE_whole_tree_analysis.json")
gene_expand_family <- data.frame(og=gene_family$expansion, stringsAsFactors = FALSE)
gene_expand_family$KO <- getSingleReactionFormula(panID_KO0$ko,panID_KO0$ortho_id,gene_expand_family$og)
gene_expand_family$pathway <- getSingleReactionFormula(KO_pathway0$V2.y,KO_pathway0$V2.x,gene_expand_family$KO)
#gene_expand_family <- gene_expand_family[gene_expand_family$pathway!="NA",]
gene_expand_family$panID <- getSingleReactionFormula(panID_KO0$representative,panID_KO0$ortho_id,gene_expand_family$og)
# as most og could find sce OGs, then we can conduct the enrichment analysis based on sce
sce_gene <- gene_expand_family$panID[str_detect(gene_expand_family$panID,"Saccharomyces_cerevisiae")]
sce_gene <- str_replace(sce_gene,"Saccharomyces_cerevisiae@","")
print(paste0(sce_gene, collapse=","))


gene_contract_family <- data.frame(og=gene_family$contraction, stringsAsFactors = FALSE)
gene_contract_family$KO <- getSingleReactionFormula(panID_KO0$ko,panID_KO0$ortho_id,gene_contract_family$og)
gene_contract_family$pathway <- getSingleReactionFormula(KO_pathway0$V2.y,KO_pathway0$V2.x,gene_contract_family$KO)



# second comparsion
# whole genome branch  vs non whole genome branch
gene_family <- fromJSON(file="cafe_data/CAFE_WGD_node_analysis.json")
gene_expand_family <- data.frame(og=gene_family[[2]]$expansion, stringsAsFactors = FALSE)
gene_expand_family$KO <- getSingleReactionFormula(panID_KO0$ko,panID_KO0$ortho_id,gene_expand_family$og)
gene_expand_family$pathway <- getSingleReactionFormula(KO_pathway0$V2.y,KO_pathway0$V2.x,gene_expand_family$KO)


gene_contract_family <- data.frame(og=gene_family[[2]]$contraction, stringsAsFactors = FALSE)
gene_contract_family$KO <- getSingleReactionFormula(panID_KO0$ko,panID_KO0$ortho_id,gene_contract_family$og)
gene_contract_family$pathway <- getSingleReactionFormula(KO_pathway0$V2.y,KO_pathway0$V2.x,gene_contract_family$KO)



# expanded gene in branch levels
gene_family <- fromJSON(file="cafe_data/CAFE_branch_clade_analysis.json")
