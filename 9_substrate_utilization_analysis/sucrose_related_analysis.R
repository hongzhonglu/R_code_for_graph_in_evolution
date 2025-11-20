#---------------------------------
# analyze the relation between genotype and phenotype for sucrose
#---------------------------------
library(ape)
library(hongR)
library(ggplot2)
library(readxl)
library(stringr)
library(readr)
library("Biostrings")
library(tidyverse)

# input interested OG
OG_K01193 <- read_excel("data_substrate_metabolism/OG_K01193.xlsx")

# get the OG id based on pan id in the above table
og_pan <- read_tsv("data/representatives.tsv")

OG_K01193$OG <- getSingleReactionFormula(og_pan$ortho_id,og_pan$representative,OG_K01193$query)


# input the species data
yeast_species <- read_excel("data_for_tree/genome_summary_332_yeasts_heat_Ethanol_updated_02_20.xlsx", sheet = "332 budding yeasts")
yeast_species <- yeast_species[, c(1,2,3,4,5)]


# analyze the existence of OG in each yeast species
# input each OG one by one
og_dir <- "/Users/xluhon/Documents/linux_project/ortholog_343/cds/"
for (i in 1: nrow(OG_K01193)) {
  
  print(OG_K01193$OG[i])
  OG0 <- OG_K01193$OG[i]
  OG_in <- paste(og_dir, OG0, "_code.fasta",sep = "")
  
  if (file.exists(OG_in)){
  
  
    # firstly check the existence of files
    s <- readDNAStringSet(OG_in)
    # extract species
    species_inf <- s@ranges@NAMES
    # extract species id
    species_inf0 <- sub("@.*", "", species_inf)
    species_inf1 <- data.frame(table(species_inf0), stringsAsFactors = FALSE)
    yeast_species[OG0] <- 0
    yeast_species[OG0] <- getSingleReactionFormula(species_inf1$Freq,species_inf1$species_inf0,yeast_species$old_species_id)
    yeast_species[,  OG0] <- as.numeric( unlist(yeast_species[, OG0]))
  }
  
}


# input the trait data
trait_yeast <- read_excel("data_substrate_metabolism/trait_experiment_data_cell.xlsx")

yeast_species$Sucrose <- getSingleReactionFormula(trait_yeast$Sucrose, trait_yeast$updated_species, yeast_species$`Species name`)
yeast_species$Raffinose <- getSingleReactionFormula(trait_yeast$Raffinose, trait_yeast$updated_species, yeast_species$`Species name`)



# input the blast result
blastnewbasedAll_yeasts2 <- read_excel("data_substrate_metabolism/blastnewbasedAll_yeasts2.xlsx")
blastnewbasedAll_yeasts20 <- filter(blastnewbasedAll_yeasts2, pident >= 30)
blastnewbasedAll_yeasts20$species_ID <- sub("@.*", "", blastnewbasedAll_yeasts20$hitID)
species_all <- data.frame(table(blastnewbasedAll_yeasts20$species_ID), stringsAsFactors = FALSE)






