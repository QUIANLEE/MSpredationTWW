# Load necessary packages

library(tidyverse)
library(ggvenn)
library(RColorBrewer)
library(dplyr)
library(tibble)
library(MASS)
library(reshape2)
library(Hmisc)
library(igraph)
library(ggplot2)
library(ggpubr)
library(ggsignif)

# Clear workspace
rm(list = ls())

############################# Conventional Soil Analysis ##########################
setwd("C:/Users/User/OneDrive/Desktop/qian/R/CYC-filter0.01/18dna-cs")

# Load metadata and OTU table
metadata18 <- read.table("group.txt", sep = "\t", header = TRUE)
otu18 <- read.table("featuretable.txt", sep = "\t", header = TRUE, check.names = FALSE, row.names = 1) %>% 
  t() %>% as.data.frame() %>% rownames_to_column("sample")

# Merge and reshape table
mergedtable18 <- merge(otu18, metadata18, by = "sample") %>%
  select(-sample, -group) %>% relocate(Replicate, .before = everything()) %>% 
  t() %>% as.data.frame()

colnames(mergedtable18) <- as.character(mergedtable18[1, ])
mergedtable18 <- mergedtable18[-1, ]

mergedtableall <- mergedtable18 %>% mutate_all(as.numeric) %>% t() %>% as.data.frame()

# Remove columns with all 0 or 1
mergedtableall <- mergedtableall[, colSums(mergedtableall == 0) != nrow(mergedtableall)]
mergedtableall <- mergedtableall[, colSums(mergedtableall == 1) != nrow(mergedtableall)] %>% t() %>% as.data.frame()

# Long format transformation
mergedtableall <- rownames_to_column(mergedtableall, var = "ASV")
mergedtableall_longer <- pivot_longer(mergedtableall, -ASV)
mergedtableall_filter <- filter(mergedtableall_longer, value != 0) %>% rename(Replicate = name)

mergedtable <- merge(mergedtableall_filter, metadata18, by = "Replicate") %>%
  select(-sample, -Replicate)

# Create Venn lists
venn_list_cs <- list(
  CSn_0 = filter(mergedtable, group == "CSn-0")$ASV,
  CSn_14 = filter(mergedtable, group == "CSn-14")$ASV,
  CSn_28 = filter(mergedtable, group == "CSn-28")$ASV,
  CSc_14 = filter(mergedtable, group == "CSc-14")$ASV,
  CSc_28 = filter(mergedtable, group == "CSc-28")$ASV
)

# Generate Venn plots
vennn_CS_14 <- ggvenn(
  venn_list_cs, columns = c("CSn_14", "CSc_14"), show_percentage = FALSE, 
  fill_color = c("#F9C00F", "#9AC9DB"), text_size = 8, auto_scale = TRUE, set_name_size = 0
) +
  ggtitle("Conventional soil") +
  theme(plot.title = element_text(size = 30, hjust = 0.5, family = "Arial")) +
  annotate("text", x = 0.6, y = 0.8, label = "CHX", size = 8, family = "Arial")

vennn_CS_28 <- ggvenn(
  venn_list_cs, columns = c("CSn_28", "CSc_28"), show_percentage = FALSE,
  fill_color = c("#F9C00F", "#9AC9DB"), text_size = 8, auto_scale = TRUE, set_name_size = 0
) +
  ggtitle("Conventional soil") +
  theme(plot.title = element_text(size = 30, hjust = 0.5, family = "Arial")) +
  annotate("text", x = 0.6, y = 0.8, label = "CHX", size = 8, family = "Arial")


############################# Organic Soil Analysis ##########################
setwd("C:/Users/User/OneDrive/Desktop/qian/R/CYC-filter0.01/18dna-os")

metadata18 <- read.table("group.txt", sep = "\t", header = TRUE) %>% select(-starts_with("X"))
otu18 <- read.table("featuretable.txt", sep = "\t", header = TRUE, check.names = FALSE, row.names = 1) %>% 
  t() %>% as.data.frame() %>% rownames_to_column("sample")

mergedtable18 <- merge(otu18, metadata18, by = "sample") %>% 
  select(-sample, -group) %>% relocate(Replicate, .before = everything()) %>%
  t() %>% as.data.frame()

colnames(mergedtable18) <- as.character(mergedtable18[1, ])
mergedtable18 <- mergedtable18[-1, ]
mergedtableall <- mergedtable18 %>% mutate_all(as.numeric) %>% t() %>% as.data.frame()

# Remove all-0 or all-1 columns
mergedtableall <- mergedtableall[, colSums(mergedtableall == 0) != nrow(mergedtableall)]
mergedtableall <- mergedtableall[, colSums(mergedtableall == 1) != nrow(mergedtableall)] %>% t() %>% as.data.frame()

mergedtableall <- rownames_to_column(mergedtableall, var = "ASV")
mergedtableall_longer <- pivot_longer(mergedtableall, -ASV)
mergedtableall_filter <- filter(mergedtableall_longer, value != 0) %>% rename(Replicate = name)

mergedtable <- merge(mergedtableall_filter, metadata18, by = "Replicate") %>% select(-sample, -Replicate)

# Create Venn lists
venn_list_os <- list(
  OSn_0 = filter(mergedtable, group == "OSn-0")$ASV,
  OSn_14 = filter(mergedtable, group == "OSn-14")$ASV,
  OSn_28 = filter(mergedtable, group == "OSn-28")$ASV,
  OSc_14 = filter(mergedtable, group == "OSc-14")$ASV,
  OSc_28 = filter(mergedtable, group == "OSc-28")$ASV
)

vennn_OS_14 <- ggvenn(
  venn_list_os, columns = c("OSn_14", "OSc_14"), show_percentage = FALSE,
  fill_color = c("#F9C00F", "#9AC9DB"), text_size = 8, auto_scale = TRUE, set_name_size = 0
) +
  ggtitle("Organic soil") +
  theme(plot.title = element_text(size = 30, hjust = 0.5, family = "Arial")) +
  annotate("text", x = 0.75, y = 0.5, label = "CHX", size = 8, family = "Arial")

vennn_OS_28 <- ggvenn(
  venn_list_os, columns = c("OSn_28", "OSc_28"), show_percentage = FALSE,
  fill_color = c("#F9C00F", "#9AC9DB"), text_size = 8, auto_scale = TRUE, set_name_size = 0
) +
  ggtitle("Organic soil") +
  theme(plot.title = element_text(size = 30, hjust = 0.5, family = "Arial")) +
  annotate("text", x = 0.65, y = 0.6, label = "CHX", size = 8, family = "Arial")

# Combine all Venn diagrams into a single figure
combined_plot <- ggarrange(
  vennn_OS_14, vennn_CS_14, 
  vennn_OS_28 + theme(plot.title = element_blank()),
  vennn_CS_28 + theme(plot.title = element_blank()),
  ncol = 2, nrow = 2, common.legend = FALSE
)

# Save combined plot
ggsave("C:/Users/User/OneDrive/Desktop/qian/R/paper1_figs/merged_venn.jpg", plot = combined_plot, device = "jpeg", dpi = 600, width = 32, height = 25, units = "cm")
ggsave("C:/Users/User/OneDrive/Desktop/qian/R/paper1_figs/merged_venn.tiff", plot = combined_plot, device = "tiff", dpi = 600, width = 32, height = 25, units = "cm")
