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
dn_ds_function <- list()
all_file <- list.files("data/function_site_dn_ds")
for (x in all_file){
  infile <- paste("data/function_site_dn_ds/", x, sep = "")
  site <- read_csv(infile)
  x1 <- str_replace(x, ".csv", "")
  print(x1)
  dn_ds_function[[x1]] <- as.vector(unlist(site[, x1]))
}


dn_ds_df = data.frame(dn_ds=c(), type=c())
for (x in all_file){
  infile <- paste("data/function_site_dn_ds/", x, sep = "")
  site <- read_csv(infile)
  x1 <- str_replace(x, ".csv", "")
  print(x1)
  dn_ds_0 <- site[, x1]
  dn_ds_0$type <- x1
  colnames(dn_ds_0) <- c("dn_ds","type")
  dn_ds_df <- rbind.data.frame(dn_ds_df, dn_ds_0)
}


# input the site classification data
site_classification <- read_excel("data/site classification.xlsx")
site_type <- unique(dn_ds_df$type)
site_type_df <- data.frame(site=site_type, class=NA, stringsAsFactors = FALSE)

site_type_df$class <- getSingleReactionFormula(site_classification$class, site_classification$Note1,site_type_df$site)
site_type_df$class[site_type_df$site=="all_site"] <- "All"
site_type_df$class[site_type_df$site=="interface"] <- "3D structure"
site_type_df$combine <- paste(site_type_df$class, site_type_df$site, sep = "_")
site_type_df$combine <- str_replace_all(site_type_df$combine, "All_","1_") %>%
  str_replace_all(., "Sites_","2_") %>% str_replace_all(., "Molecule processing_","3_") %>%
  str_replace_all(., "Regions_","4_") %>% str_replace_all(., "Amino acids modifications_","5_") %>%
  str_replace_all(., "Secondary structure_","6_") %>% str_replace_all(., "3D structure","7_")
site_type_df$combine <- str_to_lower(site_type_df$combine)
# update the site name in dn_ds_df
dn_ds_df0 <- left_join(dn_ds_df,site_type_df, c("type" = "site"))

# plot
dn_ds_df0 %>%
  ggplot(aes(x=combine,y=dn_ds, fill=combine)) +
  stat_boxplot(geom ='errorbar', width = 0.25) + # add caps
  geom_boxplot(outlier.colour = "gray75", outlier.shape = 1) +
  ylim(0, 1) +
  theme_bw() +                                       # control background and border
  theme(axis.line = element_line(colour = "black"),  # control background and border
        panel.grid.major = element_blank(),          # control background and border   
        panel.grid.minor = element_blank(),          # control background and border
        panel.border = element_blank(),              # control background and border
        panel.background = element_blank()) +        # control background and border
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = c(0.85, 0.2)) +
  theme(axis.text=element_text(size=16, family="Arial"),
        axis.title=element_text(size=20,family="Arial"),
        legend.text = element_text(size=10, family="Arial")) +
  theme(legend.position = "none") +
  xlab("") + ylab("Site-wise dN/dS")
ggsave(out <- paste('result/','dN_dS_distribution_based_function','.svg', sep = ""), width=15, height=5, dpi=300)

paste0(unique(dn_ds_df0$combine),collapse = ";")

# t.test
t.test(dn_ds_function[["Active site"]], dn_ds_function[["all_site"]])
t.test(dn_ds_function[["Site"]], dn_ds_function[["Peptide"]])

  #t.test(dn_ds_function[["Site"]], dn_ds_function[["all_site"]])
  #t.test(dn_ds_function[["Binding site"]], dn_ds_function[["all_site"]])
# wilcon.test
wilcox.test(dn_ds_function[["Active site"]], dn_ds_function[["all_site"]], alternative = "two.sided")
wilcox.test(dn_ds_function[["Site"]], dn_ds_function[["Peptide"]], alternative = "two.sided")