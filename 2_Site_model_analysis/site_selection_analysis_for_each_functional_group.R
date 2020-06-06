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


ggplot(dn_ds_df, aes(dn_ds, fill = type, colour = type)) +
  geom_density(alpha = 0.1) +
  xlim(0, 2)

dn_ds_df %>%
  ggplot(aes(x=type,y=dn_ds, fill=type, alpha=0.5)) +
    geom_boxplot(alpha = 0.1) +
    ylim(0, 1) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(legend.position = c(0.85, 0.2)) +
    theme(axis.text=element_text(size=12, family="Arial"),
        axis.title=element_text(size=12,family="Arial"),
        legend.text = element_text(size=10, family="Arial")) +
    ggtitle('') +
    theme(legend.position = "none") +
    theme(panel.background = element_rect(fill = "white", color="black", size = 1))


t.test(dn_ds_function[["Site"]], dn_ds_function[["all_site"]])

