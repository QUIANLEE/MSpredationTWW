library(pheatmap)
library(phyloseq)
library(tidyverse)

# Clear environment
rm(list = ls())

# Set working directory and load metadata
setwd("C:/Users/User/OneDrive/Desktop/qian/R/CYC-filter0.01/18dna-cs")
metadata18 <- read.table("group(+cyc)1.txt", sep = "\t", header = TRUE)

# Load and process OTU table (18S)
otu18 <- read.table("level-8.txt", sep = "\t", header = TRUE, check.names = FALSE, row.names = 1) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("sample")

# Merge OTU and metadata tables, clean columns
mergedtable18 <- merge(otu18, metadata18, by = "sample") %>% 
  dplyr::select(-sample, -group) %>% 
  relocate(Replicate, .before = everything())

# Add time labels based on experimental groupings
metadata18 <- metadata18 %>%
  mutate(time = case_when(
    grepl("CS-Day0\\(CHX-free\\)", group) ~ "Day0",
    grepl("CS-Day14\\(with CHX\\)", group) ~ "Day14(CHX)",
    grepl("CS-Day28\\(with CHX\\)", group) ~ "Day28(CHX)",
    TRUE ~ NA_character_
  )) %>% 
  dplyr::select(-group)

# Load 16S metadata and OTU table (keep only Escherichia-Shigella)
setwd("C:/Users/User/OneDrive/Desktop/qian/R/CYC-with E.coil+unknow/16dna-cs")
metadata16 <- read.table("group(+cyc)1.txt", sep = "\t", header = TRUE)
otu16 <- read.table("genus.txt", sep = "\t", header = TRUE, check.names = FALSE, row.names = 1) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("sample")

# Merge and clean 16S data
mergedtable16 <- merge(otu16, metadata16, by = "sample") %>% 
  dplyr::select(-sample, -group) %>% 
  dplyr::select(Replicate, `g__Escherichia-Shigella`)

# Merge 18S and 16S tables
mergedtableall <- merge(mergedtable18, mergedtable16, by = "Replicate") %>% 
  arrange(Replicate)

rownames(mergedtableall) <- mergedtableall$Replicate
mergedtableall <- mergedtableall[, -1]

# Scale and clean merged data
mergedtableall <- mergedtableall %>% 
  scale() %>% 
  t() %>% 
  as.data.frame()

mergedtableall[] <- lapply(mergedtableall, as.numeric)
mergedtableall[is.na(mergedtableall)] <- 0

# Clean taxon names
mergedtableall <- rownames_to_column(mergedtableall, var = "Row_Name")
mergedtableall$Row_Name <- mergedtableall$Row_Name %>% 
  gsub("d8__", "", .) %>% 
  gsub("st__un_d[2-7]__", "un_", .) %>% 
  gsub("g__Escherichia-Shigella", "E.coli", .)

# Load protist-E.coli correlation data
setwd("C:/Users/User/OneDrive/Desktop/qian/R/network/18cs+cyc-ecoli")
network <- read.csv("ecoli_related_taxa_trophic.csv", sep = ",", header = TRUE, row.names = 1)

# Filter for E.coli and bacterivorous protists
network_filtered <- network %>% 
  filter(Trophic %in% c("Consumer", "E.coli")) %>% 
  rename(Row_Name = Rank8)

# Filter merged table based on relevant taxa
mergedtableall <- mergedtableall %>% 
  inner_join(network_filtered %>% select(Row_Name), by = "Row_Name") %>% 
  distinct()

mergedtableall <- column_to_rownames(mergedtableall, var = "Row_Name")

# Format metadata for heatmap annotation
metadata18$group <- factor(metadata18$time, levels = c("Day0", "Day14(CHX)", "Day28(CHX)"))
rownames(metadata18) <- metadata18$Replicate
metadata18 <- metadata18[, -3]
metadata18 <- metadata18 %>% 
  dplyr::select(group) %>% 
  rename(Time = group)

# Separate E.coli from protists
group_ecoli <- mergedtableall[rownames(mergedtableall) == "E.coli", , drop = FALSE]
group_protists <- mergedtableall[rownames(mergedtableall) != "E.coli", , drop = FALSE]

# Cluster protists
clustered_protists <- group_protists[hclust(dist(group_protists))$order, ]

# Combine E.coli with clustered protists
combined_data <- rbind(group_ecoli, clustered_protists)

# Annotate rows
annotation_row <- data.frame(Taxon = rep(c("E.coli", "Bacterivorous protists"), c(1, nrow(clustered_protists))))
rownames(annotation_row) <- rownames(combined_data)

# Plot heatmap
heatmap_plot <- pheatmap(combined_data,
                         annotation_col = metadata18,
                         annotation_row = annotation_row,
                         cutree_cols = 3,
                         cellwidth = 20,
                         cellheight = 20,
                         cluster_rows = FALSE,
                         cluster_cols = FALSE,
                         gaps_col = c(5, 10),
                         gaps_row = 1)

print(heatmap_plot)

# Save heatmap
ggsave("18cs+cyc-heatmap-network-L8.jpg", plot = heatmap_plot, device = "jpeg", dpi = 600, width = 22, height = 14, units = "cm")
ggsave("18cs+cyc-heatmap-network-L8.tiff", plot = heatmap_plot, device = "tiff", dpi = 600, width = 22, height = 14, units = "cm")
