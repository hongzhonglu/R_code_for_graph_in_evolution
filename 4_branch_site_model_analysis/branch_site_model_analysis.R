# These scripts are used to produce map related to evolution project
# 2020.04.15
# Hongzhong Lu


#load package
#library(tidyverse)
library(ggplot2)
library(dplyr)
library(stringr)
library(readr)
library(readxl)
library(hongR)

trait_heat_result <- read_excel("data/branch_site_heat/enzyme_need_analysis_for_heat_with_result.xlsx")
figure5e_all_enzymes_2020_04 <- read_csv("data/branch_site_heat/figure5e_all_enzymes_2020_04.csv")
fccs_at_40 <- read_csv("data/branch_site_heat/fccs_at_40_100_models_figure5a_2020_04.csv")

trait_heat_result$fcc <- getSingleReactionFormula(fccs_at_40$mean, fccs_at_40$X1,trait_heat_result$Entry)

trait_heat_result$vote_percentage <-  getSingleReactionFormula(figure5e_all_enzymes_2020_04$vote_percentage, figure5e_all_enzymes_2020_04$locus,trait_heat_result$loucs)

trait_heat_result$fcc <- as.numeric(trait_heat_result$fcc)
trait_heat_result$vote_percentage <- as.numeric(trait_heat_result$vote_percentage)
trait_heat_result$number_select <- as.numeric(trait_heat_result$number_select)

write.csv(trait_heat_result, "data/branch_site_heat/trait_heat_result_integrated.csv")

# combine the result from fel and absrel
fel_result_summary <- read_csv("data/branch_site_heat/fel_result_summary.csv")
fel_result_summary$OG <- str_replace_all(fel_result_summary$OG,".csv","")
heat_result_all_update2 <- read_excel("data/branch_site_heat/heat_result_all_update2.xlsx", sheet = "top_200_OGs")
heat_result_all_update2 <- left_join(heat_result_all_update2, fel_result_summary, by = c("OG" = "OG"))
write.csv(heat_result_all_update2, "data/branch_site_heat/combine_fel_absrel_heat.csv")

# interesction of model and evolution analysis
library(VennDiagram)
gene_model <- figure5e_all_enzymes_2020_04$locus
gene_evolution <- heat_result_all_update2$sce[!is.na(heat_result_all_update2$sce)]
# get the intersection
common_gene <- intersect(gene_model, gene_evolution)
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




