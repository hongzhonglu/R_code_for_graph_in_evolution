library(ggplot2)
library(readr)
library(dplyr)
library(stringr)
library(hongR)
library(readxl)
library (vioplot) # for violpot
library(caroline) # for violpot

###################################################################################################
# input the clade information
library(readxl)
clade0 <- read_excel("data/343taxa_speicies-name_clade-name_color-code.xlsx")
clade0 <- clade0[,"Major_clade"]
summary <- as.data.frame(table(clade0), stringsAsFactors = FALSE)


###############################################################################################################
clade_choose <- summary[summary$Freq >=9 , ]
clade_choose0 <- clade_choose$clade0

# this is used to import the original data
# dn_ds_clade <- list()
# for (i in clade_choose0){
#  print(i)
#  in_dir <- paste("/home/luhongzhong/ortholog_", i, "/result_paml_parse/gene_dn_ds.csv", sep="")
#  dn_ds_clade[[i]] <- read.csv(in_dir, stringsAsFactors = FALSE)
#  dn_ds_clade[[i]]$clade <- i
#}

# to save the data in the current directory
# save(dn_ds_clade, file="data/dn_ds_clade.RData")
load("data/dn_ds_clade.RData")


dn_ds_merge <- dn_ds_clade[["CUG-Ala"]]
clade_choose1 <- clade_choose0[clade_choose0 != "CUG-Ala"]
for (i in clade_choose1){
  s2 <- dn_ds_clade[[i]]
  dn_ds_merge <- rbind.data.frame(dn_ds_merge, s2)
}

# analysis the merged dn_ds
dn_ds_merge <- filter(dn_ds_merge, !is.na(dN_dS))
dn_ds_merge$OG <- str_replace_all(dn_ds_merge$OG, ".out_yn00", "")
ggplot(dn_ds_merge,aes(x=dN_dS, fill=clade)) + geom_density(alpha=0.3) +
  xlim(0,0.75) +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20, face = "bold")) +
  theme(panel.background = element_rect(fill = "white", colour = "black"))
