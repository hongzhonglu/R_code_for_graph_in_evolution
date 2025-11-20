#---------------------------------
# analyze the relation between genotype and phenotype for maltose
# gene for maltose degration is under OG1179
#---------------------------------
library(ape)
library(hongR)
library(ggplot2)
library(readxl)
library(stringr)
library(readr)
library(Biostrings)
library(tidyverse)


substrate_all <- read_excel("data_substrate_metabolism/Supplementary Table 1.xlsx", sheet = "Sheet1")

maltose1 <- substrate_all
maltose0 <- filter(maltose1, Maltose=="1" | Maltose=="0" )




# analyze the existence of OG in each yeast species
# input each OG one by one
og_dir <- "/Users/xluhon/Documents/OG1179/protein_refine/"
OG <- list.files(og_dir)
for (i in 1:length(OG)) {
  print(OG[i])
  OG0 <- OG[i]
  OG_in <- paste(og_dir, OG0, sep = "")

  # firstly check the existence of files
  s <- readDNAStringSet(OG_in)
  # extract species
  species_inf <- s@ranges@NAMES
  # extract species id
  species_inf0 <- sub("@.*", "", species_inf)
  species_inf1 <- data.frame(table(species_inf0), stringsAsFactors = FALSE)
  maltose0[OG0] <- 0
  maltose0[OG0] <- getSingleReactionFormula(species_inf1$Freq, species_inf1$species_inf0, maltose0$`Original names`)
  maltose0[, OG0] <- as.numeric(unlist(maltose0[, OG0]))
}

write.table(maltose0, "data_substrate_metabolism/connect_gene_copy_num_with_trait.txt", sep = "\t", row.names = FALSE)


# connect the site with trait
OG1179_conserved_site <- read_excel("data_substrate_metabolism/OG1179_conserved_site.xlsx", sheet = "Sheet2")
OG1179_conserved_site$Maltose_auto <- getSingleReactionFormula(maltose1$Maltose, maltose1$`Original names`, OG1179_conserved_site$Species)
write.table(OG1179_conserved_site, "data_substrate_metabolism/connect_site_conservation_with_trait.txt", sep = "\t", row.names = FALSE)

