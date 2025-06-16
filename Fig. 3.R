library(phyloseq)
library(tidyverse)
library(fantaxtic)
library(magrittr)
library(ggnested)
library(knitr)
library(gridExtra)
library(MicrobiotaProcess)
library(ggplot2)
library(dplyr)
library(ggpubr)
rm(list = ls())
windowsFonts(Arial = windowsFont("Arial"))
############Conventional soil-Taxonomic profiles#######
setwd(dir = 'C:/Users/User/OneDrive/Desktop/qian/R/CYC-filter0.01/18dna-cs')

otu <- "featuretable.qza"
rep <- "rep-seqs.qza"
tree <- "rooted-tree.qza"
tax <- "taxonomy.qza"
sample <- "sample-metadata-dnacs.tsv"
ps_dada2 <- import_qiime2(otuqza=otu, taxaqza=tax,refseqqza=rep,
                          mapfilename=sample,treeio=tree)
ps <- ps_dada2
# Merge samples by a metadata variable
ps_grouped <- merge_samples(ps, "group")

# Step 2: Extract metadata and update the grouped metadata manually
# Get sample data (metadata) from the original phyloseq object
metadata <- data.frame(sample_data(ps))

# Group the metadata to summarize information like "type"
metadata_grouped <- metadata %>%
  group_by(group) %>%
  summarize(type = first(type)) %>%
  rename(Sample = group)  # Or use another suitable method to summarize "type"

# Step 2: Transform the grouped phyloseq object to relative abundance
ps_grouped_rel <- transform_sample_counts(ps_grouped, function(x) x / sum(x))

# Step 3: Convert phyloseq object to data frame
df <- psmelt(ps_grouped_rel)

# Step 4: Filter species with relative abundance < 1% and group them as "Other"
df_filtered <- df %>%
  group_by(Sample) %>%
  mutate(across(c(OTU, Kingdom, Phylum, Class, Order, Family, Genus, Species, Rank8, Rank9), ~ ifelse(Abundance < 0, "Others", .))) %>%
  group_by(Sample, OTU, Kingdom, Phylum, Class, Order, Family, Genus, Species, Rank8, Rank9) %>%
  summarize(Abundance = sum(Abundance), .groups = 'drop') %>%
  distinct()%>%
  left_join(metadata_grouped, by = "Sample")#ADD the column"type" into the df

#removing the prefix 'd4__' from each value in the column
df_filtered$Phylum <- gsub("d2__", "", df_filtered$Phylum)
df_filtered$Class <- gsub("d3__", "", df_filtered$Class)
df_filtered$Class <- gsub("c__un_d2__TSAR", "Unclassified_TSAR", df_filtered$Class)
df_filtered$Phylum <- gsub("TSAR:plas", "TSAR", df_filtered$Phylum)
df_filtered$Kingdom <- gsub("d1__Eukaryota:plas", "d1__Eukaryota", df_filtered$Kingdom)
df_filtered$Class <- gsub("Stramenopiles:plas", "Stramenopiles", df_filtered$Class)
df_filtered$Phylum <- as.factor(df_filtered$Phylum) 
unique(df_filtered$Phylum)

# Step 5: Convert the "Sample" column to a factor and set levels in the desired order
df_filtered$Sample <- gsub("CSn-14", "14", df_filtered$Sample)
df_filtered$Sample <- gsub("CSn-28", "28", df_filtered$Sample)
df_filtered$Sample <- gsub("CSn-0", "0", df_filtered$Sample)
df_filtered$Sample <- gsub("CSc-14", "14\n(CHX)", df_filtered$Sample)
df_filtered$Sample <- gsub("CSc-28", "28\n(CHX)", df_filtered$Sample)
df_filtered$Sample <- factor(df_filtered$Sample, levels = c("0", "14", "28","14\n(CHX)","28\n(CHX)"))

# Fix factor levels for the type variable
df_filtered$type <- factor(df_filtered$type, levels = c("Day0", "CHX-free", "with CHX"))
# Reorder levels for the Phylum factor to change legend order
df_filtered$Phylum <- factor(df_filtered$Phylum, levels = c("Amoebozoa", "Cryptista","Eukaryota_X","Excavata","Obazoa", "TSAR:plas","TSAR","Haptista","Others"))

###Convert data frame back to phyloseq object

###1-Create OTU Table
# Pivot the data frame to wide format to get OTU table format
otu_table_df <- df_filtered %>%
  select(OTU, Sample, Abundance) %>%
  pivot_wider(names_from = Sample, values_from = Abundance, values_fill = list(Abundance = 0))

# Convert to matrix and then to otu_table
otu_mat <- as.matrix(otu_table_df[, -1])  # Remove the OTU column for conversion
rownames(otu_mat) <- otu_table_df$OTU
otu_table_ps <- otu_table(otu_mat, taxa_are_rows = TRUE)

###2-Create Taxonomy Table
# Extract unique taxonomic information for each OTU
tax_table_df <- df_filtered %>%
  select(OTU, Kingdom, Phylum, Class, Order, Family, Genus,Species,Rank8,Rank9) %>%
  distinct() #dplyr package, is used to remove duplicate rows from a data frame or tibble based on one or more columns.

# Convert to matrix and then to tax_table
tax_mat <- as.matrix(tax_table_df[, -1])  # Remove OTU column
rownames(tax_mat) <- tax_table_df$OTU
tax_table_ps <- phyloseq::tax_table(tax_mat)

###3-Create Sample Metadata (Sample Data):
# Extract unique sample data
sample_data_df <- df_filtered %>%
  select(Sample, type) %>%
  distinct() %>%
  column_to_rownames("Sample")

# Convert to sample_data
sample_data_ps <- sample_data(sample_data_df)

###4- Create phyloseq object from components
ps_new <- phyloseq(otu_table_ps, tax_table_ps, sample_data_ps)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Now run nested_top_taxa again
top_nested <- nested_top_taxa(ps_new,
                              top_tax_level = "Phylum",
                              nested_tax_level = "Class",
                              grouping = "type",
                              n_top_taxa = 8,
                              n_nested_taxa = 10)

# Plot the relative abundances at two levels.
#pdf("18dnacs-L2+3-0.01.pdf",width = 7)
bar_cs <- plot_nested_bar(ps_obj = top_nested$ps_obj, 
                          top_level = "Phylum",
                          nested_level = "Class",
                          palette = c("Amoebozoa"= "#FDE18E","Cryptista" = "#40588F","Excavata" = "#9AC9DB", 
                                      "Haptista" = "#7F7F7F", "Obazoa" = "#EED5B7",
                                      "TSAR" = "#FFE8E5","TSAR:plas" = "#FFA500",
                                      "Eukaryota_X" = "#74DEB8","Others" = "#FFFFCC")
) +
  facet_grid(. ~ type, scales = "free", space = "free_x") +
  labs(title = "Conventional soil", x = "Time(Days)", y = "Relative Abundance", fill = "Phylum") +
  
  # Theme adjustments
  theme(
    legend.title = element_blank(),
    axis.title.x = element_blank(),  # x-axis title formatting
    axis.title.y = element_text(size = 20, color = "black", family = "Arial"),  # y-axis title formatting
    axis.text.x = element_text(size = 20, color = "black", family = "Arial", angle = 0, hjust = 0.5),  # Rotate x-axis labels
    axis.text.y = element_text(size = 20, color = "black", family = "Arial"), #y-axis text formatting
    plot.title = element_text(size = 30, hjust = 0.5, vjust = 2, family = "Arial"),  # ✅ Center and bold the top title#
    strip.text = element_blank()  # ✅ Remove facet titles
  ) +
  geom_bar(stat = "identity", position = "fill", width = 0.2) +
  guides(fill = guide_legend(ncol = 1))  # Set legend to one column

# Print the plot
print(bar_cs)

############Organic soil-Taxonomic profiles#######
setwd(dir = 'C:/Users/User/OneDrive/Desktop/qian/R/CYC-filter0.01/18dna-os')

otu <- "featuretable.qza"
rep <- "rep-seqs.qza"
tree <- "rooted-tree.qza"
tax <- "taxonomy.qza"
sample <- "sample-metadata-dnaos.tsv"
ps <- import_qiime2(otuqza=otu, taxaqza=tax, refseqqza=rep,
                    mapfilename=sample,treeio=tree)

# Merge samples by a metadata variable
ps_grouped <- merge_samples(ps, "group")

#  Extract metadata and update the grouped metadata manually
# Get sample data (metadata) from the original phyloseq object
metadata <- data.frame(sample_data(ps))

# Group the metadata to summarize information like "type"
metadata_grouped <- metadata %>%
  group_by(group) %>%
  summarize(type = first(type)) %>%
  rename(Sample = group)  # Or use another suitable method to summarize "type"

#  Transform the grouped phyloseq object to relative abundance
ps_grouped_rel <- transform_sample_counts(ps_grouped, function(x) x / sum(x))

#  Convert phyloseq object to data frame
df <- psmelt(ps_grouped_rel)

#  Filter species with relative abundance < 1% and group them as "Other"
df_filtered <- df %>%
  group_by(Sample) %>%
  mutate(across(c(OTU, Kingdom, Phylum, Class, Order, Family, Genus, Species, Rank8, Rank9), ~ ifelse(Abundance < 0, "Others", .))) %>%
  group_by(Sample, OTU, Kingdom, Phylum, Class, Order, Family, Genus, Species, Rank8, Rank9) %>%
  summarize(Abundance = sum(Abundance), .groups = 'drop') %>%
  distinct()%>%
  left_join(metadata_grouped, by = "Sample")#ADD the column"type" into the df

#removing the prefix 'd4__' from each value in the column
df_filtered$Phylum <- gsub("d2__", "", df_filtered$Phylum)
df_filtered$Class <- gsub("d3__", "", df_filtered$Class)
df_filtered$Class <- gsub("c__un_d2__TSAR", "Unclassified_TSAR", df_filtered$Class)

#  Convert the "Sample" column to a factor and set levels in the desired order
df_filtered$Sample <- gsub("OSn-14", "14", df_filtered$Sample)
df_filtered$Sample <- gsub("OSn-28", "28", df_filtered$Sample)
df_filtered$Sample <- gsub("OSn-0", "0", df_filtered$Sample)
df_filtered$Sample <- gsub("OSc-14", "14\n(CHX)", df_filtered$Sample)
df_filtered$Sample <- gsub("OSc-28", "28\n(CHX)", df_filtered$Sample)
df_filtered$Sample <- factor(df_filtered$Sample, levels = c("0", "14", "28","14\n(CHX)","28\n(CHX)"))

# Fix factor levels for the type variable
df_filtered$type <- factor(df_filtered$type, levels = c("Day0", "CHX-free", "with CHX"))
# Reorder levels for the Phylum factor to change legend order
df_filtered$Phylum <- factor(df_filtered$Phylum, levels = c("Amoebozoa", "Cryptista","Eukaryota_X","Excavata","Obazoa", "TSAR:plas","TSAR","Haptista","Others"))

###1-Create OTU Table
# Pivot the data frame to wide format to get OTU table format
otu_table_df <- df_filtered %>%
  select(OTU, Sample, Abundance) %>%
  pivot_wider(names_from = Sample, values_from = Abundance, values_fill = list(Abundance = 0))

# Convert to matrix and then to otu_table
otu_mat <- as.matrix(otu_table_df[, -1])  # Remove the OTU column for conversion
rownames(otu_mat) <- otu_table_df$OTU
otu_table_ps <- otu_table(otu_mat, taxa_are_rows = TRUE)

###2-Create Taxonomy Table
# Extract unique taxonomic information for each OTU
tax_table_df <- df_filtered %>%
  select(OTU, Kingdom, Phylum, Class, Order, Family, Genus,Species,Rank8,Rank9) %>%
  distinct() #dplyr package, is used to remove duplicate rows from a data frame or tibble based on one or more columns.

# Convert to matrix and then to tax_table
tax_mat <- as.matrix(tax_table_df[, -1])  # Remove OTU column
rownames(tax_mat) <- tax_table_df$OTU
tax_table_ps <- phyloseq::tax_table(tax_mat)

###3-Create Sample Metadata (Sample Data):
# Extract unique sample data
sample_data_df <- df_filtered %>%
  select(Sample, type) %>%
  distinct() %>%
  column_to_rownames("Sample")

# Convert to sample_data
sample_data_ps <- sample_data(sample_data_df)

###4- Create phyloseq object from components
ps_new <- phyloseq(otu_table_ps, tax_table_ps, sample_data_ps)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Now run nested_top_taxa again
top_nested <- nested_top_taxa(ps_new,
                              top_tax_level = "Phylum",
                              nested_tax_level = "Class",
                              grouping = "type",
                              top_merged_label = "Others",
                              n_top_taxa = 8,
                              n_nested_taxa = 10)

# Plot the relative abundances at two levels.
#pdf("18dnaos-L2+3-0.01.pdf",width = 7)
bar_os <- plot_nested_bar(ps_obj = top_nested$ps_obj, 
                          top_level = "Phylum",
                          nested_level = "Class",
                          palette = c("Amoebozoa"= "#FDE18E","Cryptista" = "#40588F","Excavata" = "#9AC9DB", 
                                      "Haptista" = "#7F7F7F", "Obazoa" = "#EED5B7",
                                      "TSAR" = "#FFE8E5","TSAR:plas" = "#FFA500",
                                      "Eukaryota_X" = "#74DEB8","Others" = "#FFFFCC")
) +
  facet_grid(. ~ type, scales = "free", space = "free_x")+
  labs(title = "Organic soil", x = "Time(Days)", y = "Relative Abundance", fill = "Phylum") +
  
  # Theme adjustments
  theme(
    legend.title = element_blank(),
    axis.title.x = element_blank(),  # x-axis title formatting
    axis.title.y = element_text(size = 20, color = "black", family = "Arial"),  # y-axis title formatting
    axis.text.x = element_text(size = 20, color = "black", family = "Arial", angle = 0, hjust = 0.5),  # Rotate x-axis labels
    axis.text.y = element_text(size = 20, color = "black", family = "Arial"), #y-axis text formatting
    plot.title = element_text(size = 30, hjust = 0.5, vjust = 2, family = "Arial"),  # ✅ Center and bold the top title#
    strip.text = element_blank()  # ✅ Remove facet titles
  ) +
  geom_bar(stat = "identity", position = "fill", width = 0.2) +
  guides(fill = guide_legend(ncol = 1))  # Set legend to one column

# Print the plot
print(bar_os)

############Conventional soil-functional profiles#######
setwd(dir = 'C:/Users/User/OneDrive/Desktop/qian/R/CYC-filter0.01/18dna-cs')

otu <- "featuretable.qza"
rep <- "rep-seqs.qza"
tree <- "rooted-tree.qza"
tax <- "taxonomy.qza"
sample <- "sample-metadata-dnacs.tsv"
ps <- import_qiime2(otuqza=otu, taxaqza=tax,refseqqza=rep,
                    mapfilename=sample,treeio=tree)

# Merge samples by a metadata variable
ps_grouped <- merge_samples(ps, "group")

#  Extract metadata and update the grouped metadata manually
# Get sample data (metadata) from the original phyloseq object
metadata <- data.frame(sample_data(ps))

# Group the metadata to summarize information like "type"
metadata_grouped <- metadata %>%
  group_by(group) %>%
  summarize(type = first(type)) %>%
  rename(Sample = group)  # Or use another suitable method to summarize "type"

#  Transform the grouped phyloseq object to relative abundance
ps_grouped_rel <- transform_sample_counts(ps_grouped, function(x) x / sum(x))

# Step 3: Convert phyloseq object to data frame
df <- psmelt(ps_grouped_rel)

#  Filter species with relative abundance < 1% and group them as "Other"
df_filtered <- df %>%
  group_by(Sample) %>%
  mutate(across(c(OTU, Kingdom, Phylum, Class, Order, Family, Genus, Species, Rank8, Rank9), ~ ifelse(Abundance < 0, "Others", .))) %>%
  group_by(Sample, OTU, Kingdom, Phylum, Class, Order, Family, Genus, Species, Rank8, Rank9) %>%
  summarize(Abundance = sum(Abundance), .groups = 'drop') %>%
  distinct()%>%
  left_join(metadata_grouped, by = "Sample")#ADD the column"type" into the df

#removing the prefix 'd4__' from each value in the column
df_filtered$Family <- gsub("d5__", "", df_filtered$Family)
df_filtered$Family <- gsub("f__un_d2__", "Unclassified_", df_filtered$Family)
df_filtered$Family <- gsub("f__un_d3__", "Unclassified_", df_filtered$Family)
df_filtered$Family <- gsub("f__un_d4__", "Unclassified_", df_filtered$Family)
# df_filtered$Family <- gsub("Phaeophyceae:plas", "Phaeophyceae", df_filtered$Family)
df_filtered$Phylum <- gsub("d2__", "", df_filtered$Phylum)
# TSAR:plas specifically refers to plastid-bearing lineages within TSAR.
# It highlights taxa that have retained photosynthetic plastids, often through endosymbiosis.
# This means that TSAR:plas includes organisms with photosynthetic capabilities
df_filtered$Phylum <- gsub("TSAR:plas", "TSAR", df_filtered$Phylum)
df_filtered$Kingdom <- gsub("d1__Eukaryota:plas", "d1__Eukaryota", df_filtered$Kingdom)
df_filtered$Phylum <- as.factor(df_filtered$Phylum) 
unique(df_filtered$Phylum)

STrophic <- c("Bacillariophyceae"="Phototroph",
              "Bicoecea"="Consumer",
              "Centramoebia"="Consumer",
              "Centroplasthelida_XX"="Consumer",
              "Choanoflagellatea"="Consumer",
              "Chrysophyceae"="Phototroph",#"Mixotroph",
              "Coccidiomorphea"="Parasite",
              "Colpodellidea"="Consumer",
              "Colponemidia_X"="Consumer",
              "Cryptophyceae"="Phototroph",#"Mixotroph",
              "Dinophyceae"="Phototroph",#"Mixotroph",
              "Discoba_XX"="Consumer",
              "Echinamoebida"="Consumer",
              "Elardia"="Consumer",
              "Endomyxa"="Consumer",
              "Endomyxa-Phytomyxea"="Parasite",
              "Euglenida"="Consumer",
              "Filosa-Granofilosea"="Consumer",
              "Filosa-Imbricatea"="Consumer",
              "Filosa-Sarcomonadea"="Consumer",
              "Filosa-Thecofilosea"="Consumer",
              "Flabellinia"="Consumer",
              "Globothalamea"="Consumer",
              "Gregarinomorphea"="Parasite",
              "Hemimastigophora_XX"="Consumer",
              "Heterolobosea"="Consumer",
              "Heterotrichea"="Consumer",
              "Hyphochytriomyceta"="Saprotroph",
              "Ichthyophonida"="Parasite",
              "Kinetoplastea"="Heterotroph",
              "Litostomatea"="Consumer",
              "Monothalamids"="Consumer",
              "Myxogastria_Fuscisporidia"="Consumer",
              "Novel-clade-10-12"="Consumer",
              "Nucleariidea"="Consumer",
              "Oligohymenophorea"="Consumer",
              "Opalozoa"="Consumer",
              "Peronosporomycetes"="Parasite",
              "Phaeophyceae:plas"="Phototroph",
              "Phyllopharyngea"="Consumer",
              "Pterocystida"="Consumer",
              "Sagenista"="Consumer",
              "Sainouroidea"="Consumer",
              "Spirotrichea"="Consumer",
              "Variosea"="Consumer",
              "Unclassified_TSAR"="Unknown",
              "Unclassified_Alveolata"="Unknown",
              "Unclassified_Stramenopiles"="Unknown",
              "Unclassified_Bigyra"="Unknown",
              "Unclassified_Centroplasthelida_X"="Consumer",
              "Unclassified_Ciliophora"="Consumer",
              "Unclassified_Gyrista"="Unknown",
              "Apusomonadidae"="Consumer",
              "Mediophyceae"="Phototroph",
              "Unclassified_Euglenozoa" ="Unknown",
              "Unclassified_Foraminifera" ="Unknown"
)

df_filtered <- df_filtered %>%
  mutate(Trophic = STrophic[Family])


# Step 1: Extract rows where Time == "Day0" & Treatment == "-CHX"
rows_to_duplicate <- df_filtered %>%
  filter(Sample == "CSn-0", 
         OTU == "03749dcc36d37f46aabc97c4b02c3765", 
         Abundance == 0, 
         Trophic== "Heterotroph")

# Step 2: Modify the Treatment column to "+CHX"
new_rows <- rows_to_duplicate %>%
  mutate(Phylum = "TSAR", 
         OTU = "03749dcc36d37f46aabc97c4b02c3852",
         Abundance = 0)

# Step 3: Append new rows to the original data
df_filtered <- bind_rows(df_filtered, new_rows)

# Ensure Abundance is numeric and remove any possible NA values
df_filtered <- df_filtered %>%
  mutate(Abundance = as.numeric(Abundance)) %>%
  replace_na(list(Abundance = 0))  # Replace any NA with 0
####
#Reorder levels
df_filtered$type <- factor(df_filtered$type, levels = c("Day0", "CHX-free", "with CHX"))

df_filtered$Sample <- gsub("CSn-14", "14", df_filtered$Sample)
df_filtered$Sample <- gsub("CSn-28", "28", df_filtered$Sample)
df_filtered$Sample <- gsub("CSn-0", "0", df_filtered$Sample)
df_filtered$Sample <- gsub("CSc-14", "14\n(CHX)", df_filtered$Sample)
df_filtered$Sample <- gsub("CSc-28", "28\n(CHX)", df_filtered$Sample)
df_filtered$Sample <- factor(df_filtered$Sample, levels = c("0", "14", "28","14\n(CHX)","28\n(CHX)"))


######Convert data frame back to phyloseq object
###1-Create OTU Table
# Pivot the data frame to wide format to get OTU table format
otu_table_df <- df_filtered %>%
  select(OTU, Sample, Abundance) %>%
  pivot_wider(names_from = Sample, values_from = Abundance, values_fill = list(Abundance = 0))

# Convert to matrix and then to otu_table
otu_mat <- as.matrix(otu_table_df[, -1])  # Remove the OTU column for conversion
rownames(otu_mat) <- otu_table_df$OTU
otu_table_ps <- otu_table(otu_mat, taxa_are_rows = TRUE)

###2-Create Taxonomy Table
# Extract unique taxonomic information for each OTU
tax_table_df <- df_filtered %>%
  select(OTU,Trophic,Kingdom, Phylum, Class, Order, Family, Genus,Species,Rank8,Rank9) %>%
  distinct() #dplyr package, is used to remove duplicate rows from a data frame or tibble based on one or more columns.

# Convert to matrix and then to tax_table
tax_mat <- as.matrix(tax_table_df[, -1])  # Remove OTU column
rownames(tax_mat) <- tax_table_df$OTU
tax_table_ps <- phyloseq::tax_table(tax_mat)

###3-Create Sample Metadata (Sample Data):
# Extract unique sample data
sample_data_df <- df_filtered %>%
  select(Sample, type) %>%
  distinct() %>%
  column_to_rownames("Sample")

# Convert to sample_data
sample_data_ps <- sample_data(sample_data_df)

###4- Create phyloseq object from components
ps_new <- phyloseq(otu_table_ps, tax_table_ps, sample_data_ps)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Now run nested_top_taxa again
top_nested <- nested_top_taxa(ps_new,
                              top_tax_level = "Trophic",
                              nested_tax_level = "Phylum",
                              grouping = "type",
                              n_top_taxa = 10,
                              n_nested_taxa = 10)

# Plot the relative abundances at two levels.
bar_csT <- plot_nested_bar(ps_obj = top_nested$ps_obj, 
                           top_level = "Trophic",
                           nested_level = "Phylum",
                           palette = c(Consumer = "#FF8ACE",
                                       Parasite = "#FDE18E", 
                                       Heterotroph ="#8A7B61",
                                       Saprotroph = "#FFD9D0",
                                       Phototroph = "#74DEB8", 
                                       #Mixotroph = "#D94701",
                                       Unknown="#8ACEFF")
) +
  facet_grid(. ~ type, scales = "free", space = "free_x")+
  labs(title = "Conventional soil", x = "Time(Days)", y = "Relative Abundance", fill = "Phylum") +  # Set x-axis title to an empty space (not completely empty)
  theme(
    legend.title = element_blank(),
    axis.title.x = element_text(size = 20, color = "black", family = "Arial"),  # x-axis title formatting
    axis.title.y = element_text(size = 20, color = "black", family = "Arial"),  # Add space to the right of y-axis title
    axis.text.x = element_text(size = 20, color = "black", family = "Arial", angle = 0, hjust = 0.5),  # Rotate x-axis labels
    axis.text.y = element_text(size = 20, color = "black", family = "Arial"), #y-axis text formatting
    strip.text = element_blank(),
    plot.title = element_blank())+  # ✅ Center and bold the top title#
  guides(
    fill = guide_legend(ncol = 1)  # ✅ Make legend a single column
  )+
  geom_bar(stat = "identity", position = "fill", width = 0.2) 

print(bar_csT)

############Organic soil-functional profiles#######
setwd(dir = 'C:/Users/User/OneDrive/Desktop/qian/R/CYC-filter0.01/18dna-os')
otu <- "featuretable.qza"
rep <- "rep-seqs.qza"
tree <- "rooted-tree.qza"
tax <- "taxonomy.qza"
sample <- "sample-metadata-dnaos.tsv"
ps <- import_qiime2(otuqza=otu, taxaqza=tax,refseqqza=rep,
                    mapfilename=sample,treeio=tree)

# Merge samples by a metadata variable
ps_grouped <- merge_samples(ps, "group")

#   Extract metadata and update the grouped metadata manually
# Get sample data (metadata) from the original phyloseq object
metadata <- data.frame(sample_data(ps))

# Group the metadata to summarize information like "type"
metadata_grouped <- metadata %>%
  group_by(group) %>%
  summarize(type = first(type)) %>%
  rename(Sample = group)  # Or use another suitable method to summarize "type"

#   Transform the grouped phyloseq object to relative abundance
ps_grouped_rel <- transform_sample_counts(ps_grouped, function(x) x / sum(x))

#   Convert phyloseq object to data frame
df <- psmelt(ps_grouped_rel)

#   Filter species with relative abundance < 1% and group them as "Other"
df_filtered <- df %>%
  group_by(Sample) %>%
  mutate(across(c(OTU, Kingdom, Phylum, Class, Order, Family, Genus, Species, Rank8, Rank9), ~ ifelse(Abundance < 0, "Others", .))) %>%
  group_by(Sample, OTU, Kingdom, Phylum, Class, Order, Family, Genus, Species, Rank8, Rank9) %>%
  summarize(Abundance = sum(Abundance), .groups = 'drop') %>%
  distinct()%>%
  left_join(metadata_grouped, by = "Sample")#ADD the column"type" into the df

#removing the prefix 'd4__' from each value in the column
df_filtered$Family <- gsub("d5__", "", df_filtered$Family)
df_filtered$Family <- gsub("f__un_d2__", "Unclassified_", df_filtered$Family)
df_filtered$Family <- gsub("f__un_d3__", "Unclassified_", df_filtered$Family)
df_filtered$Family <- gsub("f__un_d4__", "Unclassified_", df_filtered$Family)
df_filtered$Phylum <- gsub("d2__", "", df_filtered$Phylum)

#write.table (df_filtered,file ="df_filtered.xls", sep ="\t",col.names = NA)

STrophic <- c("Apusomonadidae"="Consumer",
              "Bacillariophyceae"="Phototroph",
              "Bicoecea"="Consumer",
              "Centramoebia"="Consumer",
              "Centroplasthelida_XX"="Consumer",
              "Choanoflagellatea"="Consumer",
              "Chrysophyceae"="Phototroph",#"Mixotroph",
              "Coccidiomorphea"="Parasite",
              "Colpodea"="Consumer",
              "Colpodellidea"="Consumer",
              "Cryptophyceae"="Phototroph",#"Mixotroph",
              "Diplonemea"="Consumer",
              "Discoba_XX"="Consumer",
              "Echinamoebida"="Consumer",
              "Elardia"="Consumer",
              "Endomyxa"="Consumer",
              "Endomyxa-Phytomyxea"="Parasite",
              "Euglenida"="Consumer",
              "Eustigmatophyceae"="Phototroph",
              "Filosa-Granofilosea"="Consumer",
              "Filosa-Imbricatea"="Consumer",
              "Filosa-Sarcomonadea"="Consumer",
              "Filosa-Thecofilosea"="Consumer",
              "Flabellinia"="Consumer",
              "Globothalamea"="Consumer",
              "Gregarinomorphea"="Parasite",
              "Hemimastigophora_XX"="Consumer",
              "Heterolobosea"="Consumer",
              "Heterotrichea"="Consumer",
              "Hyphochytriomyceta"="Saprotroph",
              "Ichthyophonida"="Parasite",
              "Kinetoplastea"="Heterotroph",
              "Litostomatea"="Consumer",
              "Mediophyceae"="Phototroph",
              "Monothalamids"="Consumer",
              "Myxogastria_Fuscisporidia"="Consumer",
              "Myxogastria_Lucisporidia"="Consumer",
              "Oligohymenophorea"="Consumer",
              "Opalozoa"="Consumer",
              "Peronosporomycetes"="Parasite",
              "Phyllopharyngea"="Consumer",
              "Pterocystida"="Consumer",
              "Sagenista"="Consumer",
              "Variosea"="Consumer",
              "Xanthophyceae"="Phototroph",
              "Unclassified_TSAR"="Unknown",
              "Unclassified_Alveolata"="Unknown",
              "Unclassified_Stramenopiles"="Unknown",
              "Unclassified_Apicomplexa"="Parasite",
              "Unclassified_Bigyra"="Unknown",
              "Unclassified_Centroplasthelida_X"="Consumer",
              "Unclassified_Cercozoa"="Heterotroph",
              "Unclassified_Ciliophora"="Consumer",
              "Unclassified_Dinoflagellata"="Unknown",
              "Unclassified_Euglenozoa"="Unknown",
              "Unclassified_Foraminifera"="Unknown",
              "Unclassified_Gyrista"="Unknown")

df_filtered <- df_filtered %>%
  mutate(Trophic = STrophic[Family])
write.table (df_filtered,file ="df_filtered-trophic.xls", sep ="\t",col.names = NA)
#Reorder levels
df_filtered$type <- factor(df_filtered$type, levels = c("Day0", "CHX-free", "with CHX"))
df_filtered$Sample <- gsub("OSn-14", "14", df_filtered$Sample)
df_filtered$Sample <- gsub("OSn-28", "28", df_filtered$Sample)
df_filtered$Sample <- gsub("OSn-0", "0", df_filtered$Sample)
df_filtered$Sample <- gsub("OSc-14", "14\n(CHX)", df_filtered$Sample)
df_filtered$Sample <- gsub("OSc-28", "28\n(CHX)", df_filtered$Sample)
df_filtered$Sample <- factor(df_filtered$Sample, levels = c("0", "14", "28","14\n(CHX)","28\n(CHX)"))


######Convert data frame back to phyloseq object
###1-Create OTU Table
# Pivot the data frame to wide format to get OTU table format
otu_table_df <- df_filtered %>%
  select(OTU, Sample, Abundance) %>%
  pivot_wider(names_from = Sample, values_from = Abundance, values_fill = list(Abundance = 0))

# Convert to matrix and then to otu_table
otu_mat <- as.matrix(otu_table_df[, -1])  # Remove the OTU column for conversion
rownames(otu_mat) <- otu_table_df$OTU
otu_table_ps <- otu_table(otu_mat, taxa_are_rows = TRUE)

###2-Create Taxonomy Table
# Extract unique taxonomic information for each OTU
tax_table_df <- df_filtered %>%
  select(OTU,Trophic,Kingdom, Phylum, Class, Order, Family, Genus,Species,Rank8,Rank9) %>%
  distinct() #dplyr package, is used to remove duplicate rows from a data frame or tibble based on one or more columns.

# Convert to matrix and then to tax_table
tax_mat <- as.matrix(tax_table_df[, -1])  # Remove OTU column
rownames(tax_mat) <- tax_table_df$OTU
tax_table_ps <- phyloseq::tax_table(tax_mat)

###3-Create Sample Metadata (Sample Data):
# Extract unique sample data
sample_data_df <- df_filtered %>%
  select(Sample, type) %>%
  distinct() %>%
  column_to_rownames("Sample")

# Convert to sample_data
sample_data_ps <- sample_data(sample_data_df)

###4- Create phyloseq object from components
ps_new <- phyloseq(otu_table_ps, tax_table_ps, sample_data_ps)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Now run nested_top_taxa again
top_nested <- nested_top_taxa(ps_new,
                              top_tax_level = "Trophic",
                              nested_tax_level = "Phylum",
                              grouping = "type",
                              n_top_taxa = 8,
                              n_nested_taxa = 10)

# Plot the relative abundances at two levels.
bar_osT <- plot_nested_bar(ps_obj = top_nested$ps_obj, 
                              top_level = "Trophic",
                              nested_level = "Phylum",
                              palette = c(Consumer = "#FF8ACE",
                                          Parasite = "#FDE18E", 
                                          Heterotroph = "#8A7B61",
                                          Saprotroph = "#FFD9D0",
                                          Phototroph = "#74DEB8", 
                                          Unknown = "#8ACEFF")
) +
  facet_grid(. ~ type, scales = "free", space = "free_x") +
  labs(title = "Organic soil", x = "Time(Days)",  y = "Relative Abundance", fill = "Phylum") +  # Set x-axis title to an empty space (not completely empty)
  theme(
    legend.title = element_blank(),
    axis.title.x = element_text(size = 20, color = "black", family = "Arial"),  # x-axis title formatting
    axis.title.y = element_text(size = 20, color = "black", family = "Arial"),  # Add space to the right of y-axis title
    axis.text.x = element_text(size = 20, color = "black", family = "Arial", angle = 0, hjust = 0.5),  # Rotate x-axis labels
    axis.text.y = element_text(size = 20, color = "black", family = "Arial"), #y-axis text formatting
    strip.text = element_blank(),
    plot.title = element_blank())+  # Center and bold the top title
  guides(
    fill = guide_legend(ncol = 1)  # ✅ Make legend a single column
  )+
  geom_bar(stat = "identity", position = "fill", width = 0.2) 

print(bar_osT)
############Code to Merge Two Boxplots into One Figure#####################
# Group A & B with a shared legend (from bar_cs)
group_ab <- ggarrange(
  bar_os + theme(legend.position = "none",
                 plot.margin = margin(20, 10, 10, 40),
                 plot.title = element_blank(),
                 axis.title.y = element_text(margin = margin(r = 35))),
  bar_cs + theme(plot.margin = margin(20, 10, 10, 40),
                 plot.title = element_blank(),
                 axis.title.y = element_text(margin = margin(r = 35))),
  # labels = c("A", "B"),
  # font.label = list(size = 25, face = "bold", family = "Arial"),
  ncol = 2, nrow = 1,
  align = "v",
  common.legend = TRUE,
  legend = "right"
)

# Group C & D with a shared legend (from bar_osT)
group_cd <- ggarrange(
  bar_osT + theme(legend.position = "none",
                  plot.margin = margin(20, 10, 10, 40),
                  axis.title.y = element_text(margin = margin(r = 35))),
  bar_csT + theme(plot.margin = margin(20, 10, 10, 40),
                  axis.title.y = element_text(margin = margin(r = 35))),
  # labels = c("C", "D"),
  # font.label = list(size = 25, face = "bold", family = "Arial"),
  ncol = 2, nrow = 1,
  align = "v",
  common.legend = TRUE,
  legend = "right"
)

# Combine the two groups vertically
combined_plot <- ggarrange(
  group_ab,
  group_cd,
  ncol = 1, nrow = 2
)

# Print the combined plot
print(combined_plot)

#legend = "right"  # Legend position: "top", "bottom", "left", or "right"
# Print the combined plot
print(combined_plot)

# Save the combined plot
ggsave("C:/Users/User/OneDrive/Desktop/qian/R/paper1_figs/merged_otu_trophic1.jpg", plot = combined_plot, device = "jpeg", 
       dpi = 600, width = 37, height = 34, units = "cm")

ggsave("C:/Users/User/OneDrive/Desktop/qian/R/paper1_figs/merged_otu_trophic1.tiff", plot = combined_plot, device = "tiff", 
       dpi = 600, width = 37, height = 34, units = "cm")







