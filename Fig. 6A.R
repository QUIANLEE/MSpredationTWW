library(dplyr)
library(tibble)
library(MASS)
library(reshape2)
library(Hmisc)
library(igraph)
library(ggplot2)
library(ggpubr)
library(ggsignif)

# Clear environment
rm(list = ls())

# -------------------- Load and Process 18S Data --------------------
setwd('C:/Users/User/OneDrive/Desktop/qian/R/CYC-filter0.01/18dna-os')
metadata18 <- read.table("group(-cyc).txt", sep = "\t", header = TRUE)
otu18 <- read.table("level-8.txt", sep = "\t", header = TRUE, check.names = FALSE, row.names = 1) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample")

mergedtable18 <- merge(otu18, metadata18, by = "sample") %>%
  dplyr::select(-c(sample, group)) %>%
  relocate(Replicate, .before = everything()) %>%
  t() %>%
  as.data.frame()

colnames(mergedtable18) <- as.character(mergedtable18[1, ])
mergedtable18 <- mergedtable18[-1, ] %>% mutate_all(as.numeric)

# CLR transformation
mergedtable18_pseu <- mergedtable18 + 1
mergedtable18_clr <- microbiome::transform(mergedtable18_pseu, "clr") %>% as.data.frame()

# -------------------- Load and Process 16S (E. coli) Data --------------------
setwd('C:/Users/User/OneDrive/Desktop/qian/R/CYC-with E.coil+unknow/16dna-os')
metadata16 <- read.table("group(-cyc).txt", sep = "\t", header = TRUE)
otu16 <- read.table("genus.txt", sep = "\t", header = TRUE, check.names = FALSE, row.names = 1) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample")

mergedtable16 <- merge(otu16, metadata16, by = "sample") %>%
  dplyr::select(-sample) %>%
  relocate(Replicate, .before = everything()) %>%
  t() %>%
  as.data.frame()

colnames(mergedtable16) <- as.character(mergedtable16[1, ])
mergedtable16 <- mergedtable16[-1, ] %>% mutate_all(as.numeric)

mergedtable16_pseu <- mergedtable16 + 1
mergedtable16_clr <- microbiome::transform(mergedtable16_pseu, "clr") %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  dplyr::select(`g__Escherichia-Shigella`) %>%
  t() %>%
  as.data.frame()

# Standardize taxon labels
mergedtable18_clr <- rownames_to_column(mergedtable18_clr, var = "Row_Name") %>%
  mutate(Row_Name = gsub("d8__|st__un_d[2-7]__", "un_", Row_Name))
mergedtable16_clr <- rownames_to_column(mergedtable16_clr, var = "Row_Name") %>%
  mutate(Row_Name = gsub("g__Escherichia-Shigella", "E.coli", Row_Name))

mergedtable18_clr <- column_to_rownames(mergedtable18_clr, var = "Row_Name")
mergedtable16_clr <- column_to_rownames(mergedtable16_clr, var = "Row_Name")

mergedtable18_clr <- mergedtable18_clr %>% t() %>% as.data.frame() %>% rownames_to_column("sample")
mergedtable16_clr <- mergedtable16_clr %>% t() %>% as.data.frame() %>% rownames_to_column("sample")

# Merge 16S and 18S datasets
mergedtableall <- merge(mergedtable18_clr, mergedtable16_clr, by = "sample") %>%
  t() %>%
  as.data.frame()
colnames(mergedtableall) <- as.character(mergedtableall[1, ])
mergedtableall <- mergedtableall[-1, ] %>% mutate_all(as.numeric) %>% as.matrix()

# -------------------- Filter by prevalence (10%) --------------------
dt <- mergedtableall
binary_dt <- dt
binary_dt[binary_dt > 1] <- 1
prevalent_rows <- which(rowSums(binary_dt)/ncol(binary_dt) > 0.1)
mergedtableall <- mergedtableall[prevalent_rows, ]

# -------------------- Spearman Correlation --------------------
sp.cor <- rcorr(t(mergedtableall), type = "spearman")
r.cor <- sp.cor$r
p.cor <- sp.cor$P

write.csv(r.cor, "rvaluecor_level8.csv")
write.csv(p.cor, "pvaluecor_level8.csv")

# Adjust p-values
p.adj <- p.adjust(p.cor, method = "BH")

# Filter based on p < 0.05 and |r| > 0.6
r.cor[p.cor >= 0.05 | abs(r.cor) < 0.6] <- 0
p.cor[p.cor >= 0.05 | abs(r.cor) < 0.6] <- 0

write.csv(r.cor, "Rvalue_filter_level8.csv")
write.csv(p.cor, "Pvalue_filter_level8.csv")

# Replace NA with 1 for safety
p.cor[is.na(p.cor)] <- 1

# -------------------- E. coli-specific Correlation --------------------
r.ecoli <- r.cor %>% as.data.frame() %>%
  filter(abs(E.coli) > 0) %>%
  t() %>% as.data.frame() %>%
  filter(abs(E.coli) > 0)

p.ecoli <- p.cor %>% as.data.frame() %>%
  filter(abs(E.coli) != 0) %>%
  t() %>% as.data.frame() %>%
  filter(abs(E.coli) != 0)

diag(r.ecoli) <- 0
diag(p.ecoli) <- NA

write.csv(r.ecoli, "R.ecoli_level8_sig.csv")
write.csv(p.ecoli, "P.ecoli_level8_sig.csv")

# -------------------- Network Construction --------------------
r.ecoli <- read.csv("R.ecoli_level8_sig.csv", row.names = 1)
A <- graph_from_adjacency_matrix(as.matrix(r.ecoli), weight = TRUE, mode = "undirected")
A <- simplify(A)
E(A)$width <- abs(E(A)$weight) * 5

# Add node attributes: abundance
abund_data <- rowSums(mergedtableall)
abund_df <- as.data.frame(abund_data)
V(A)$Abundance <- log10(abund_df[V(A)$name, , drop = FALSE]$abund_data)

# Add taxonomy
taxonomy <- read.table("level2-8.txt", sep = "\t", header = TRUE, row.names = 1)
taxonomy <- rownames_to_column(taxonomy, var = "Row_Name")
abund_df <- rownames_to_column(abund_df, var = "Row_Name")
taxonomy <- taxonomy %>% semi_join(abund_df, by = "Row_Name")
taxonomy <- column_to_rownames(taxonomy, var = "Row_Name")
V(A)$G <- as.character(taxonomy[V(A)$name, ])

# Add label and degree
V(A)$label <- V(A)$name
V(A)$degree <- degree(A)

# Save the network
write_graph(A, "os-cyc-L8-ecoli.graphml", format = "graphml")
write_graph(A, "os-cyc-L8-ecoli.gml", format = "gml")
