# note------------
# here we try to analyze the functions of genes with positive selected sites based on pathways.
# 2020-06-02
# Hongzhong Lu


library(ggplot2)
library(readr)
library(dplyr)
library(stringr)
library(hongR)
library(readxl)

# align by macse for pruned Ortholog group
fubar_macse_prune <- read_csv("data/site_model_data/fubar_02_08.csv")
fel_macse_prune <- read_csv("data/site_model_data/fel_macse_prune_3400.csv")


# align by macse for pruned Ortholog group
fubar_guidance_prune <- read_csv("data/site_model_data/fubar_guidance_prune.csv")
fel_guidance_prune <- read_csv("data/site_model_data/fel_guidance_prune.csv")


# align by macse for unpruned Ortholog group
fubar_guidance_unprune <- read_csv("data/site_model_data/fubar_guidance_unprune.csv")


# filter based on the sites number
site_num <- 1
fubar_macse_prune0 <- filter(fubar_macse_prune, positive_num >= site_num)
fel_macse_prune0 <- filter (fel_macse_prune, site_pvalue_0.1 >= site_num)
#fel_macse_prune0 <- filter (fel_macse_prune, site_pvalue_0.05 >= site_num)
fubar_guidance_prune0 <- filter(fubar_guidance_prune, positive_num >= site_num)
fel_guidance_prune0 <- filter (fel_guidance_prune, site_pvalue_0.1 >= site_num)
#fel_guidance_prune0 <- filter (fel_guidance_prune, site_pvalue_0.05 >= site_num)
fubar_guidance_unprune0 <- filter(fubar_guidance_unprune, positive_num >= site_num)

# insect analysis
length(intersect(fubar_guidance_prune0$OG, fubar_guidance_unprune0$OG))
length(intersect(fubar_guidance_prune0$OG, fubar_macse_prune0$OG))



# prepare the plot
type <- c("fubar_macse_prune", "fubar_guidance_prune","fubar_guidance_unprune","fel_macse_prune", 
          "fel_guidance_prune")
number <- c(length(fubar_macse_prune0$OG), length(fubar_guidance_prune0$OG), length(fubar_guidance_unprune0$OG),
            length(fel_macse_prune0$OG), length(fel_guidance_prune0$OG))
all_result <- data.frame(source=type, num=number, stringsAsFactors = FALSE)
# remove the fel_macse_prune as it is only analyze the result from 3400 OGs
all_result0 <- filter(all_result, source !="fel_macse_prune")

ggplot(data=all_result0, aes(x=source, y=num)) +
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = c(0.85, 0.2)) +
  theme(axis.text=element_text(size=12, family="Arial"),
        axis.title=element_text(size=12,family="Arial"),
        legend.text = element_text(size=10, family="Arial")) +
  ggtitle('') +
  theme(legend.position = "none") #+
  #theme(panel.background = element_rect(fill = "white", color="black", size = 1)


ggplot(data=all_result, aes(x=source, y=num)) +
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = c(0.85, 0.2)) +
  theme(axis.text=element_text(size=12, family="Arial"),
        axis.title=element_text(size=12,family="Arial"),
        legend.text = element_text(size=10, family="Arial")) +
  ggtitle('') +
  theme(legend.position = "none") #+
#theme(panel.background = element_rect(fill = "white", color="black", size = 1)



