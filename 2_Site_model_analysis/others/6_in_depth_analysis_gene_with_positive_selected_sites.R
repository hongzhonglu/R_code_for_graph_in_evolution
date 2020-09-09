# note------------
# This is old result based on MACSE alignment, which will be removed laterly.
# 2020-06-02
# Hongzhong Lu


library(ggplot2)
library(readr)
library(dplyr)
library(stringr)
library(hongR)
library(readxl)
library (vioplot) # for violpot
library(caroline) # for violpot


# now the result only contains the FUBAR analysis result. In the future, the fel analysis results will be included.
###################################################################################################
# positive selected genes initial analysis
# read all OG information
gene_dn_ds_all_new <- read.table("result/gene_dn_ds_all_new.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# positive selected gene analysis across 343 yeast species
fubar_02_08 <- read_csv("data/fubar_02_08.csv")
select_gene_dn_ds <- merge.data.frame(gene_dn_ds_all_new, fubar_02_08, by.x = "OG", by.y="OG")
select_gene_dn_ds$Gene_number <- NA
select_gene_dn_ds$Gene_number[select_gene_dn_ds$positive_num >= 3 ] <- "d over_3_sites"
select_gene_dn_ds$Gene_number[select_gene_dn_ds$positive_num ==2 ] <- "c 2_sites"
select_gene_dn_ds$Gene_number[select_gene_dn_ds$positive_num ==1 ] <- "b 1_sites"
select_gene_dn_ds$Gene_number[select_gene_dn_ds$positive_num ==0 ] <- "a no_sites"
select_gene_dn_ds$Gene_number <- as.factor(select_gene_dn_ds$Gene_number)
# plot
ggplot(select_gene_dn_ds,aes(x=Gene_number, y=dN_dS, fill=Gene_number)) + geom_boxplot() +
  ylim(0,1) +
  #theme(panel.background = element_rect(fill = "white", colour = "black")) + # generate the whole border
  theme_bw() +                                       # control background and border
  theme(axis.line = element_line(colour = "black"),  # control background and border
        panel.grid.major = element_blank(),          # control background and border   
        panel.grid.minor = element_blank(),          # control background and border
        panel.border = element_blank(),              # control background and border
        panel.background = element_blank()) +        # control background and border

  theme(legend.position = c(0.2, 0.8)) +
  theme(axis.text=element_text(size=12,face="bold", family="Arial"),
        axis.title=element_text(size=24,face="bold", family="Arial"),
        legend.text = element_text(size = 13, family = "Arial")) +  ggtitle('') +
        theme(legend.position = "none")



SG1 <-  filter(select_gene_dn_ds, Gene_number == "a no_sites")
SG2 <-  filter(select_gene_dn_ds, Gene_number == "b 1_sites")
SG3 <-  filter(select_gene_dn_ds, Gene_number == "c 2_sites")
SG4 <-  filter(select_gene_dn_ds, Gene_number == "d over_3_sites")

# comparsion values from different groups
t.test(SG4$dN_dS,SG1$dN_dS)
t.test(SG4$dN_dS,SG2$dN_dS)
t.test(SG4$dN_dS,SG3$dN_dS)

t.test(SG2$dN_dS,SG1$dN_dS)
t.test(SG3$dN_dS,SG1$dN_dS)

# calculate the number of OGs with different selection number
# plot
ggplot(select_gene_dn_ds, aes(Gene_number)) +
  geom_bar(fill = "seagreen4") +
  theme(legend.position = c(0.2, 0.8)) +
  #theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  
  theme_bw() +                                       # control background and border
  theme(axis.line = element_line(colour = "black"),  # control background and border
        panel.grid.major = element_blank(),          # control background and border   
        panel.grid.minor = element_blank(),          # control background and border
        panel.border = element_blank(),              # control background and border
        panel.background = element_blank()) +        # control background and border
  
  theme(axis.text=element_text(size=12,face="bold", family="Arial"),
        axis.title=element_text(size=24,face="bold", family="Arial"),
        legend.text = element_text(size = 13, family = "Arial"))