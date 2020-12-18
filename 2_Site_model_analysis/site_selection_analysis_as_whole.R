# Plot the distribution of site level dN/dS
# 2020.06.07
# Hongzhong Lu

#load package
#library(tidyverse)
library(ggplot2)
library(dplyr)
library(stringr)
library(readr)
library(readxl)
library(hongR)

# some general function used for evolution analysis
filterSiteBasedGap <- function(OGID, max_gap_ratio = 0.3) {

  # this function is used to output the interesting sites using for the evolution analysis
  # based on the protein 3D structures
  # now this function is only suitable for result from FUBAR
  # also only used for panID from sce.


  # input
  # panID0, panID of each OG group
  # output
  # coordinate dataframe of the selected sites
  evolution_dir <- "/Users/luho/Google Drive/R application and code/protein 3D structure QC and QA/Evolution analysis/Evolution_data/result_site/"
  #column_choose <- c("gap_ratio", "site_ratio", "alpha", "beta", "dN_dS", "Prob[alpha>beta]", "Prob[alpha<beta]")
  column_choose <- c("gap_ratio", "dN_dS", "Prob[alpha>beta]", "Prob[alpha<beta]")
  # read and furthe reduce the data size
  s0 <- paste(evolution_dir, OGID, sep = "")

  evolution_df <- read_csv(s0) %>% select(., column_choose)
  # filter out gaps based on the ref-site
  evolution_df1 <- evolution_df
  # filter out gaps based on the gap ratio
  gap_percentage <- max_gap_ratio # be careful in this step! shoud observe how this parameter affect the output
  evolution_df2 <- filter(evolution_df1, gap_ratio <= gap_percentage)
  return(evolution_df2)
} # not used function here!

# all OG with conserved site data from sce
#data_dir <- "/Users/luho/Google Drive/R application and code/protein 3D structure QC and QA/Evolution analysis/Evolution_data/result_site/"
#OG_sce <- list.files(data_dir)

# test 
#dn_ds_all <- vector()

#for (i in OG_sce){
#  print(i)
#  og0 <- i
#  result0 <- filterSiteBasedGap(OGID=og0, max_gap_ratio = 0.3)
#  dn_ds <- result0$dN_dS
#  dn_ds_all <- c(dn_ds_all, dn_ds)
#}

#save(dn_ds_all, file = "data/dn_ds_all.RData")






# the dn_ds_all is calculated based on the site model, see the code from the evolution analysis
load("data/dn_ds_all.RData")

# density plot is not good
plot(density(dn_ds_all), xlim=c(0,5))

# bar plot
g <- c()
g[1] <- length(dn_ds_all[dn_ds_all<=0.2])
g[2] <- length(dn_ds_all[dn_ds_all > 0.2 & dn_ds_all <= 0.4])
g[3] <- length(dn_ds_all[dn_ds_all > 0.4 & dn_ds_all <= 0.6])
g[4] <- length(dn_ds_all[dn_ds_all > 0.6 & dn_ds_all <= 0.8])
g[5] <- length(dn_ds_all[dn_ds_all > 0.8 & dn_ds_all <= 1])
g[6] <- length(dn_ds_all[dn_ds_all > 1 & dn_ds_all <= 1.2])
g[7] <- length(dn_ds_all[dn_ds_all > 1.2 & dn_ds_all <= 1.4])
g[8] <- length(dn_ds_all[dn_ds_all > 1.4 & dn_ds_all <= 1.6])
g[9] <- length(dn_ds_all[dn_ds_all > 1.6 & dn_ds_all <= 1.8])
g[10] <- length(dn_ds_all[dn_ds_all > 1.8 & dn_ds_all <= 2])
g[11] <- length(dn_ds_all[dn_ds_all > 2 & dn_ds_all <= 2.2])
g[12] <- length(dn_ds_all[dn_ds_all > 2.2 & dn_ds_all <= 2.4])
g[13] <- length(dn_ds_all[dn_ds_all > 2.4 & dn_ds_all <= 2.6])
g[14]<- length(dn_ds_all[dn_ds_all > 2.6 & dn_ds_all <= 2.8])
g[15] <- length(dn_ds_all[dn_ds_all > 2.8 & dn_ds_all <= 3])
g[16]<- length(dn_ds_all[dn_ds_all > 3])
g_name <- paste("g", 1:16, sep = "")
g_name2 <- c("0-0.2", "0.2-0.4","0.4-0.6","0.6-0.8","0.8-1.0","1.0-1.2","1.2-1.4","1.4-1.6","1.6-1.8","1.8-2.0","2.0-2.2","2.2-2.4",
             "2.4-2.6","2.6-2.8","2.8-3.0","3.0-")
df <- data.frame(group= g_name, num = g, range0 = g_name2)
df$num_log <- log10(df$num)
df$order0 <- 1:16

df <- df[order(df$order0),]
df$range0 <-factor(df$range0, levels=df$range0)


ggplot(data=df, aes(x=range0, y=num_log)) +
  geom_bar(stat="identity", color="blue",fill="blue", alpha=0.3) +
  xlab("dN/dS") + 
  ylab("log10 (Count)") + 
  theme(panel.background = element_rect(fill = "white", color="black", size = 1),
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  theme(axis.text=element_text(size=16, family="Arial"),
        axis.title=element_text(size=20, family="Arial"),
        legend.text = element_text(size=20, family="Arial")) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave(out <- paste('result/','all_site_level_dN_dS_distribution','.svg', sep = ""), width=12, height=6, dpi=600)




ggplot(data=df, aes(x=range0, y=num)) +
  geom_bar(stat="identity", fill="steelblue") +
  xlab("dN/dS") + 
  ylab("Count") + 
  theme(panel.background = element_rect(fill = "white", color="black", size = 1),
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  theme(axis.text=element_text(size=12, family="Arial"),
        axis.title=element_text(size=16, family="Arial"),
        legend.text = element_text(size=20, family="Arial")) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))


##################### not used ################
library(plotrix)
x <- c(1:5, 6.9, 7)
y <- 2^x
from <- 33
to <- 110
plot(x, y, type="b", xlab="index", ylab="value")
gap.plot(x, y, gap=c(from,to), type="b", xlab="index", ylab="value")
axis.break(2, from, breakcol="snow", style="gap")
axis.break(2, from*(1+0.02), breakcol="black", style="slash")
axis.break(4, from*(1+0.02), breakcol="black", style="slash")
axis(2, at=from)


twogrp <- df$num
gap.barplot(twogrp,gap=c(5500,6000),xlab="Index",
            ylab="Group values",main="Barplot with gap")