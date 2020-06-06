# These scripts are used to produce map related to evolution project
# 2019.2.20
# Hongzhong Lu

#load package
#library(tidyverse)
library(ggplot2)
library(dplyr)
library(stringr)
library(readr)
library(readxl)
library(hongR)

# input the data 
# some general function used for evolution analysis
getConsevedSite <- function(OGID, site_type = "conserved", max_gap_ratio = 0.3) {

  # this function is used to output the interesting sites using for the evolution analysis
  # based on the protein 3D structures
  # now this function is only suitable for result from FUBAR
  # also only used for panID from sce.


  # input
  # OGID
  # output
  # coordinate dataframe of the selected sites


  # input the evolution information based on the OGID
  # it should be careful that result from different methods will be also different in the column name
  evolution_dir <- "/Users/luho/Google Drive/R application and code/protein 3D structure QC and QA/Evolution analysis/Evolution_data/result_conserved_site/"
  column_choose <- c("gap_ratio", "site_ratio", "conserved_sites", "refsite")
  # read and furthe reduce the data size
  s0 <- paste(evolution_dir, OGID, sep = "")

  if (!file.exists(s0)) {
    print("NO evolution data file can be found!!")
    return(NA)
  }

  # check whether OGID exist??
  if (length(OGID) == 0) {
    print("NO Ortholog ID can be found!!")
    return(NA)
  } else {
    evolution_df <- read_csv(s0)
    all_col <- colnames(evolution_df)
    if (!("refsite" %in% all_col)) {
      return(NA)
    } else {
      evolution_df <- evolution_df %>% select(., column_choose)
      # filter out gaps based on the ref-site
      evolution_df1 <- filter(evolution_df, refsite != "-")
      # add the coordinates information of the reference sequence
      coordinate_ref <- 1:nrow(evolution_df1)
      evolution_df1$coordinate_ref <- coordinate_ref
      # filter out gaps based on the gap ratio
      gap_percentage <- max_gap_ratio # be careful in this step! shoud observe how this parameter affect the output
      evolution_df2 <- filter(evolution_df1, gap_ratio <= gap_percentage)
      # change the data formula
      if (site_type == "conserved") {
        positive_site <- evolution_df2[evolution_df2$conserved_sites != "{}", ]
        if (nrow(positive_site) >= 1) {
          interest_site <- data.frame(Pos = positive_site$coordinate_ref, Ref = positive_site$refsite, Alt = "X", stringsAsFactors = FALSE)
          return(interest_site)
        } else {
          print("No select sites! Please check.")
          return(NA)
        }
      }
    }
  }
} # This function is not used here


# all OG with conserved site data from sce
#data_dir <- "/Users/luho/Google Drive/R application and code/protein 3D structure QC and QA/Evolution analysis/Evolution_data/result_conserved_site/"
#OG_sce <- list.files(data_dir)

#conserve_site <- vector(length = length(OG_sce))
#t <- 0
#for (i in OG_sce){
#  print(i)
#  t <- t+1
#  s0 <- getConsevedSite(OGID=i, site_type = "conserved", max_gap_ratio = 0.3)
#  site_num <- nrow(s0)
#  if(length(site_num)){
#    conserve_site[t] <- site_num } else{
#    conserve_site[t] <- NA
#  }
#}

#conserve_site_inf <- data.frame(OG=OG_sce, conserve_site=conserve_site, stringsAsFactors = FALSE)
#conserve_site_inf$OG <- str_replace_all(conserve_site_inf$OG,"_conserved_site.csv","")



#save(conserve_site_inf, file = "data/conserve_site_inf.RData")


# the dn_ds_all is calculated based on the site model, see the code from the evolution analysis
load("data/conserve_site_inf.RData")



# id mapping
id_mapping <- read_tsv("data/representatives.tsv")
conserve_site_inf$locus <- getSingleReactionFormula(id_mapping$representative,id_mapping$ortho_id,conserve_site_inf$OG)
conserve_site_inf$locus <- str_replace_all(conserve_site_inf$locus, "Saccharomyces_cerevisiae@","")

# protein length
sgd_gene_annotation <- read_tsv("data/yeast_protein_sequence_SGD.tsv")

conserve_site_inf$protein_length <- getSingleReactionFormula(sgd_gene_annotation$`Sequence length`,sgd_gene_annotation$`Gene Systematic name`,conserve_site_inf$locus)
conserve_site_inf$protein_length <- as.numeric(conserve_site_inf$protein_length)
conserve_site_inf$ratio_conserved <- conserve_site_inf$conserve_site / conserve_site_inf$protein_length
conserve_site_inf0 <- conserve_site_inf[!is.na(conserve_site_inf$ratio_conserved),]
plot(density(conserve_site_inf0$ratio_conserved))
plot(density(conserve_site_inf0$conserve_site))

ggplot(conserve_site_inf0, aes(x=ratio_conserved)) +
  geom_density(fill="lightblue", alpha=1) +
  xlab("Ratio of conserved sites") + 
  ylab("Density") +
  theme(axis.text=element_text(size=12,face="bold", family="Arial"),
        axis.title=element_text(size=16,face="bold", family="Arial"),
        legend.text = element_text(size=20,face="bold", family="Arial")) +
  ggtitle('') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1))



# possible enrichment analysis
conserve_site_inf1 <- conserve_site_inf0[order(conserve_site_inf0$ratio_conserved),]
group1 <- paste0(conserve_site_inf1$locus[1:100],collapse = " ")
group2 <- paste0(conserve_site_inf1$locus[(length(conserve_site_inf1$locus)-99):length(conserve_site_inf1$locus)],collapse = " ")

