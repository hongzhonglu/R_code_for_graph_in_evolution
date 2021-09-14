# These scripts are used to produce map related to evolution project
# This result is based on the two parallel calculation
# heat-tolerance as the main trait
# Compare the model prediction and selection analysis
# Compare the proteomics and selection analysis
# This result is based on the one calculation
# 2021.07.23 updated on the new MAC
# Hongzhong Lu


#load package
#library(tidyverse)
library(ggplot2)
library(dplyr)
library(stringr)
library(readr)
library(readxl)
library(hongR)

# input the data from the first calculation
trait_heat_result <- read_csv("/Users/xluhon/Documents/branch_site_heat/heat_result_all_update2.csv")

# filter based on the select species and select clade
interest_OG <- filter(trait_heat_result, Select_species>=10 & select_clade_all >=3)


# input the data from second calculation
trait_heat_result2 <- read_csv("/Users/xluhon/Documents/branch_site_heat_2nd/heat_result_all_update2.csv")
trait_heat_result2$Select_species <- as.numeric(trait_heat_result2$Select_species)
# filter based on the select species and select clade
interest_OG2 <- trait_heat_result2[trait_heat_result2$Select_species >= 10 & trait_heat_result2$select_clade_all >=3,]
interest_OG2 <- interest_OG2[!is.na(interest_OG2$OG),]


# calculate the common
common_OGs <- intersect(interest_OG$OG, interest_OG2$OG)
filter_result <- trait_heat_result[trait_heat_result$OG %in% common_OGs, ]
write.table(filter_result, "result/commom_select_gene_for_heat_in_two_independent_calculation.txt", row.names = FALSE, sep = "\t")


# plot a heatmap based on the species number under different combination of clade and species with positive selection
selection_analysis <- trait_heat_result[, c("OG","select_clade_all", "Select_species")]
selection_analysis <- filter(selection_analysis, Select_species>=6 & select_clade_all >=2)
# change the data format
selection_analysis$combine <- paste(selection_analysis$select_clade_all, selection_analysis$Select_species,sep="&")
selection_analysis_table <- as.data.frame(table(selection_analysis$combine),stringsAsFactors = FALSE)
# build a matrix
ss <- matrix(, nrow = 3, ncol = 12)
for(i in 1:3){
  for(j in 1:12){
  i0 <- i + 1
  j0 <- j + 5
  search_index <- paste(i0, j0, sep = "&")
  print(search_index)
  value <- selection_analysis_table$Freq[which(selection_analysis_table$Var1==search_index)]
  if(length(value)){
    ss[i,j] <- value} 
  else{
    ss[i,j] <- 0
  }
  }
}
rownames(ss) <- paste(2:4,"_clades", sep = "")
colnames(ss) <- paste(6:17,"_species", sep = "")


# heatmap for Figure 6B
library(superheat)
require(scales)
library(pheatmap)
library(RColorBrewer)
library(grid)

ss_df <- as.data.frame(ss)
ss_df <- ss_df[order(row.names(ss_df),decreasing = TRUE), ]
#ss_df <- ss_df[, order(colnames(ss_df),decreasing = FALSE)]
new_col_name <- paste(17:6, "_species", sep = "")
ss_df <- ss_df[, new_col_name]
pheatmap(ss_df,
         method = c("pearson"),
         clustering_method = "complete",
         treeheight_row = 40,
         treeheight_col = 40,
         cluster_row = FALSE,
         cluster_col = FALSE,
         show_rownames = T,
         show_colnames = T,
         legend = T,
         fontsize = 14,
         color = colorRampPalette(c("white", "SandyBrown", "firebrick3"))(100))

# when output as pdf, the best size 5.25 x 2.26 inches
write.csv(ss_df, "result/dataset_for_heat_map_of_positive_selected_gene_related_to_heat_tolerance.csv")






figure5e_all_enzymes_2020_04 <- read_csv("data/branch_site_heat/figure5e_all_enzymes_2020_04.csv")
fccs_at_40 <- read_csv("data/branch_site_heat/fccs_at_40_100_models_figure5a_2020_04.csv")

trait_heat_result$fcc <- getSingleReactionFormula(fccs_at_40$mean, fccs_at_40$X1,trait_heat_result$Entry)

trait_heat_result$vote_percentage <-  getSingleReactionFormula(figure5e_all_enzymes_2020_04$vote_percentage, figure5e_all_enzymes_2020_04$locus,trait_heat_result$loucs)

trait_heat_result$fcc <- as.numeric(trait_heat_result$fcc)
trait_heat_result$vote_percentage <- as.numeric(trait_heat_result$vote_percentage)
trait_heat_result$number_select <- as.numeric(trait_heat_result$number_select)

write.csv(trait_heat_result, "data/branch_site_heat/trait_heat_result_integrated.csv")


# combine the result from fel and absrel
# the fel_result is based on the site model
fel_result_summary <- read_csv("data/branch_site_heat/fel_result_summary.csv")
fel_result_summary$OG <- str_replace_all(fel_result_summary$OG,".csv","")
heat_result_all_update2 <- read_excel("data/branch_site_heat/heat_result_all_update2.xlsx", sheet = "top_200_OGs")
heat_result_all_update2 <- left_join(heat_result_all_update2, fel_result_summary, by = c("OG" = "OG"))
write.csv(heat_result_all_update2, "data/branch_site_heat/combine_fel_absrel_heat.csv")





# interesction of model and evolution analysis
library(VennDiagram)
gene_model <- figure5e_all_enzymes_2020_04$locus
gene_evolution <- trait_heat_result$sce[trait_heat_result$OG %in% common_OGs]
gene_evolution <- gene_evolution[!is.na(gene_evolution)]
# get the intersection
common_gene <- intersect(gene_model, gene_evolution)

paste0(gene_evolution, collapse = ",")
paste0(common_gene, collapse = ",")
paste0(gene_model, collapse = ",")
#plot the graph
VennDiagram::venn.diagram(x= list(Model = gene_model , Evolution = gene_evolution), 
                          filename = "result/Compare_model_prediction_and_evolution_analysis.png", height = 1000, width = 1000,resolution =300, imagetype="png", col="transparent",
                          fill=c("blue","red"),alpha = 0.50, cex=0.45, cat.cex=0.45)






# omics analysis under 38 C
up_regulated_gene <- read_excel("data/branch_site_heat/Multi-omics from our own lab.xlsx", 
                                           sheet = "up_regulation")
up_regulated_gene0 <- up_regulated_gene$up_regulated

down_regulated_gene <- read_excel("data/branch_site_heat/Multi-omics from our own lab.xlsx", 
                                sheet = "down_regulation")
down_regulated_gene0 <- down_regulated_gene$down_regulated

VennDiagram::venn.diagram(x= list(Evolution = gene_evolution, Up_regulate=up_regulated_gene0, Down_regulated=down_regulated_gene0), 
                          filename = "result/Combine_evolution_analaysis_with_omics_data.png", height = 1000, width = 1000,resolution =300, imagetype="png", col="transparent",
                          fill=c("blue","red","green"),alpha = 0.50, cex=0.45, cat.cex=0.45)








## compare the species number in the filtered OGs and the remaining reference OGs
OG_df <- trait_heat_result[,c("OG")]
OG_df$Type <- NA
OG_df$Type[OG_df$OG %in% common_OGs] <- "Top selected"
OG_df$Type[!(OG_df$OG %in% common_OGs)] <- "Reference"
# input the species number for each OG
ortholog_occurance <- read_tsv("data/ortholog_occurence_num_all.tsv")
OG_df$Species_num <- getSingleReactionFormula(ortholog_occurance$species_num,ortholog_occurance$ID,OG_df$OG)
OG_df$Species_num <- as.numeric(OG_df$Species_num)
# plot a box graph
ggplot(data=OG_df, aes(x=Type, y=Species_num)) +
  geom_boxplot(fill="steelblue") +
  #theme(axis.text.x = element_text(angle = 0, hjust = 1)) + 
  theme(panel.background = element_rect(fill = "white", color="black", size = 1),
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  theme(axis.text=element_text(size=16, family="Arial"),
        axis.title=element_text(size=20, family="Arial"),
        legend.text = element_text(size=20, family="Arial"))


##################################################
# combine the result from experimental data
##################################################
gene_in_vivo <- read_excel("data/Ogataea_polymorpha_gene_for_heat_tolerance.xlsx")
sce_ortholog <- gene_in_vivo$`Systematic name of Saccharomyces cerevisiae ortholog`


# find the common gene between evolution and experimental data
common_evolution_and_validation <- intersect(sce_ortholog, gene_evolution)








# Note: the followed script is not used!
# analyze the evolution rate difference between these verified genes
Ogataea_polymorpha_dN_dS <- read_csv("data/Ogataea_polymorpha_dN_dS.csv")
Ogataea_polymorpha_dN_dS$OG <- str_replace_all(Ogataea_polymorpha_dN_dS$OG,"_yn00.csv","")


sce_gene_summary <- read_excel("data/sce_gene_summary.xlsx")
sce_gene_summary$sce_gene <- str_replace_all(sce_gene_summary$sce_gene, "Saccharomyces_cerevisiae@","") %>%
  str_replace_all(.,"\\'", "") %>%  str_replace_all(.,"\\[", "") %>% str_replace_all(.,"\\]", "") 
sce_gene_summary0 <- sce_gene_summary[, c("OrthologID", "sce_gene")]
sce_gene_summary1 <- splitAndCombine(sce_gene_summary0$sce_gene,sce_gene_summary0$OrthologID,",")
colnames(sce_gene_summary1) <- c("sce_gene","OGid")
sce_gene_summary1$sce_gene <- str_trim(sce_gene_summary1$sce_gene, side = "both")







