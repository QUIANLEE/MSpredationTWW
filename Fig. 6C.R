library(dplyr)
library(tibble)
library(MASS)
library(reshape2)
library(Hmisc)
library(igraph)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(microbiome)
library(psych)

# Clear environment and load 18S data (absolute OTU table)
rm(list = ls())
setwd("C:/Users/User/OneDrive/Desktop/qian/R/CYC-filter0.01/18dna-cs")
metadata18 <- read.table("group(-cyc).txt", sep = "\t", header = TRUE)
otu18 <- read.table("level-8.txt", sep = "\t", header = TRUE, check.names = FALSE, row.names = 1) %>% 
  t() %>% as.data.frame() %>% rownames_to_column("sample")
mergedtable18 <- merge(otu18, metadata18, by = "sample") %>%
  dplyr::select(-sample) %>% relocate(Replicate, .before = everything()) %>%
  t() %>% as.data.frame()
colnames(mergedtable18) <- as.character(mergedtable18[1, ])
mergedtable18 <- mergedtable18[-1, ] %>% mutate_all(as.numeric)

# CLR transformation
mergedtable18_pseu <- mergedtable18 + 1
mergedtable18_clr <- microbiome::transform(mergedtable18_pseu, "clr") %>% as.data.frame()

# Load 16S data and retain only Escherichia-Shigella
setwd("C:/Users/User/OneDrive/Desktop/qian/R/CYC-with E.coil+unknow/16dna-cs")
metadata16 <- read.table("group(-cyc).txt", sep = "\t", header = TRUE)
otu16 <- read.table("genus.txt", sep = "\t", header = TRUE, check.names = FALSE, row.names = 1) %>% 
  t() %>% as.data.frame() %>% rownames_to_column("sample")
mergedtable16 <- merge(otu16, metadata16, by = "sample") %>%
  dplyr::select(-sample) %>% relocate(Replicate, .before = everything()) %>%
  t() %>% as.data.frame()
colnames(mergedtable16) <- as.character(mergedtable16[1, ])
mergedtable16 <- mergedtable16[-1, ] %>% mutate_all(as.numeric)

# CLR transformation
mergedtable16_pseu <- mergedtable16 + 1
mergedtable16_clr <- microbiome::transform(mergedtable16_pseu, "clr") %>%
  as.data.frame() %>% t() %>% as.data.frame() %>%
  dplyr::select(`g__Escherichia-Shigella`) %>% t() %>% as.data.frame()

# Rename taxa
mergedtable18_clr <- rownames_to_column(mergedtable18_clr, var = "Row_Name")
mergedtable16_clr <- rownames_to_column(mergedtable16_clr, var = "Row_Name")
mergedtable18_clr$Row_Name <- gsub("d8__|st__un_d[2-7]__", "un_", mergedtable18_clr$Row_Name)
mergedtable16_clr$Row_Name <- gsub("g__Escherichia-Shigella", "E.coli", mergedtable16_clr$Row_Name)

mergedtable18_clr <- column_to_rownames(mergedtable18_clr, var = "Row_Name") %>% t() %>% as.data.frame() %>% rownames_to_column("sample")
mergedtable16_clr <- column_to_rownames(mergedtable16_clr, var = "Row_Name") %>% t() %>% as.data.frame() %>% rownames_to_column("sample")

# Merge 16S and 18S tables
mergedtableall <- merge(mergedtable18_clr, mergedtable16_clr, by = "sample") %>%
  t() %>% as.data.frame()
colnames(mergedtableall) <- as.character(mergedtableall[1, ])
mergedtableall <- mergedtableall[-1, ] %>% mutate_all(as.numeric)

# Filter taxa detected in >10% of samples
dt <- mergedtableall
dt[dt > 1] <- 1
no <- which(rowSums(dt) / ncol(dt) > 0.1)
mergedtableall <- mergedtableall[no, ] %>% as.matrix()

# Spearman correlation
sp.cor <- rcorr(t(mergedtableall), type = "spearman")
r.cor <- sp.cor$r
p.cor <- sp.cor$P

# Filter by correlation thresholds
r.cor[p.cor >= 0.05 | abs(r.cor) < 0.6] <- 0
write.csv(r.cor, 'C:/Users/User/OneDrive/Desktop/qian/R/network/18cs-cyc-ecoli/Rvalue_filter_level8.csv')
p.cor[p.cor >= 0.05 | abs(r.cor) < 0.6] <- 0
write.csv(p.cor, 'C:/Users/User/OneDrive/Desktop/qian/R/network/18cs-cyc-ecoli/Pvalue_filter_level8.csv')
p.cor[is.na(p.cor)] <- 1

# Keep only taxa correlated with E.coli
r.ecoli <- r.cor %>% as.data.frame() %>% filter(abs(E.coli) > 0) %>% t() %>% as.data.frame() %>% filter(abs(E.coli) > 0)
p.ecoli <- p.cor %>% as.data.frame() %>% filter(abs(E.coli) != 0) %>% t() %>% as.data.frame() %>% filter(abs(E.coli) != 0)
diag(r.ecoli) <- 0
write.csv(r.ecoli, 'C:/Users/User/OneDrive/Desktop/qian/R/network/18cs-cyc-ecoli/R.ecoli_level8_sig.csv')
diag(p.ecoli) <- NA
write.csv(p.ecoli, 'C:/Users/User/OneDrive/Desktop/qian/R/network/18cs-cyc-ecoli/P.ecoli_level8_sig.csv')

# Filter abundance table based on taxa in r.ecoli
r.ecoli <- read.csv("C:/Users/User/OneDrive/Desktop/qian/R/network/18cs-cyc-ecoli/R.ecoli_level8_sig.csv", row.names = 1)
mergedtableall <- mergedtableall %>% as.data.frame() %>% mutate(X = row.names(.)) %>%
  filter(X %in% rownames(r.ecoli)) %>% dplyr::select(-X) %>% mutate_all(as.numeric) %>% as.matrix()

# Format correlation matrix
colnames(r.ecoli) <- gsub("\\.", "-", colnames(r.ecoli))
r.ecoli <- as.matrix(r.ecoli)

# Create graph
A <- graph_from_adjacency_matrix(r.ecoli, weight = TRUE, mode = "undirected") %>% simplify()
E(A)$width <- abs(E(A)$weight) * 5

# Add abundance info
data <- rowSums(mergedtableall)
data1 <- as.data.frame(data)
vertex_names <- V(A)$name
df_igraph_size <- data1[vertex_names, , drop = FALSE]
df_igraph_size2 <- log10(df_igraph_size)
V(A)$Abundance <- df_igraph_size2$data

# Add taxonomy info
data2 <- read.table("C:/Users/User/OneDrive/Desktop/qian/R/network/18cs-cyc-ecoli/level2-8.txt", sep = "\t", header = TRUE, check.names = FALSE, row.names = 1)
data1 <- rownames_to_column(data1, var = "Row_Name")
data2 <- rownames_to_column(data2, var = "Row_Name") %>% semi_join(data1, by = "Row_Name")
data2 <- column_to_rownames(data2, var = "Row_Name")
data1 <- column_to_rownames(data1, var = "Row_Name")
df_igraph_col <- data2[V(A)$name, ]
V(A)$G <- as.character(df_igraph_col)

# Add labels and degree
V(A)$label <- V(A)$name
V(A)$degree <- degree(A)

# Plot network
plot(A, vertex.size = 20, vertex.label.cex = 0.8)

# Save network files
write_graph(A, "C:/Users/User/OneDrive/Desktop/qian/R/network/18cs-cyc-ecoli/cs-cyc-L8-ecoli.graphml", format = "graphml")
write_graph(A, "C:/Users/User/OneDrive/Desktop/qian/R/network/18cs-cyc-ecoli/cs-cyc-L8-ecoli.gml", format = "gml")
