library(ggplot2)
library(readr)
library(dplyr)
library(stringr)
library(hongR)

# dN_dS for 1011 sce sequence project
# The result is not good
gene_dn_ds <-  read_csv("data/gene_dn_ds_complex.csv") %>% filter(., !is.na(dN_dS)) %>%
  filter(., dN_dS < 3 & dN_dS > 0)
# plot the column density plot
ggplot(gene_dn_ds, aes(dN_dS)) +
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "blue",  bins = 50) +
  xlim(0, 3) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  geom_density(colour = "black", alpha = 0.3, fill = "grey50") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20, face = "bold")) +
  geom_vline(
    xintercept = 1, linetype = "dotted",
    color = "red", size = 1
  ) +
  geom_vline(
    xintercept = 3, linetype = "dotted",
    color = "red", size = 1
  )
