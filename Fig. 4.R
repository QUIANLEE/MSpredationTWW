library(DESeq2)
library(ggplot2)
library(ggrepel)
library(readxl) 
library(tidyverse)
library(dplyr)
library(tibble)

########################Conventional soil####################################################
rm(list = ls())
setwd(dir = 'C:/Users/User/OneDrive/Desktop/qian/R/CYC-filter0.01/18dna-cs')

groupdata <- read.table("group(-cyc28 vs. +cyc28).txt", sep = "\t", header = T)

level5 <- read.delim("level-5.txt",check.names = FALSE, row.names = 1)%>%
  as.data.frame(header=NULL,check.names=T) %>%
  select("41","42","43","44","45","46","47","48","9","10")
level5 <- rownames_to_column(level5, var = "Row_Name")

level5$Row_Name <- gsub("d5__", "", level5$Row_Name)
level5$Row_Name <- gsub("d3__", "", level5$Row_Name)
level5$Row_Name <- gsub("f__un_d2__", "d2_", level5$Row_Name)
level5$Row_Name <- gsub("f__un_d3__", "d3_", level5$Row_Name)
level5$Row_Name <- gsub("f__un_d4__", "d4_", level5$Row_Name)
level5$Row_Name <- gsub("f__un_", "Un_", level5$Row_Name)
level5<- column_to_rownames(level5, var = "Row_Name")

# Step 1: Create DESeqDataSet 
dds <- DESeqDataSetFromMatrix(countData = level5,colData = groupdata,design = ~group)

# Step 2: Run differential abundance analysisdds <- DESeq(dds) 
dds <- DESeq(dds)
suppressMessages(dds)

# Step 3: Extract results
res <- results(dds,contrast = c('group','CSn-28','CSc-28'))
summary(res) 

# Step 4: Organize results 
deseq_res <- as.data.frame(res[order(res$padj), ]) 
deseq_res$otu_id <- rownames(deseq_res)
#write.table(deseq_res[c(7, 1:6)], 'C:/Users/User/OneDrive/Desktop/qian/R/CYC-filter0.01/18dna-cs/deseq_pic/DESeq2-level5(+cyc28.vs.-cyc28).txt', row.names = FALSE, sep = '\t', quote = FALSE)

# Step 5: Classify significance
for (i in 1:nrow(deseq_res)) { 
  if 
  (deseq_res[i,'log2FoldChange'] 
   %in% NA | abs(deseq_res[i,'log2FoldChange']) < 0.95) 
    deseq_res[i,'select_change'] <- 'n' else deseq_res[i,'select_change'] <- 'y'   
    if 
    (deseq_res[i,'padj'] 
     %in% NA | abs(deseq_res[i,'padj']) >= 0.05) 
      deseq_res[i,'select_padj'] <- 'n' else deseq_res[i,'select_padj'] <- 'y' 
      deseq_res[i,'select'] <- paste(deseq_res[i,'select_change'], deseq_res[i,'select_padj'], sep 
                                     = '') 
} 
deseq_res$sig <-"Not Significant"
deseq_res$sig[(deseq_res$log2FoldChange >= 0.95)&(deseq_res$padj<0.05)]<-"Higher in CHX-free"
deseq_res$sig[(deseq_res$log2FoldChange < (-0.95))&(deseq_res$padj<0.05)]<-"Higher in with-CHX"
deseq_res$name <-deseq_res$otu_id
deseq_res$name[(deseq_res$sig=="Not Significant")]<-NA

# Step 6: Volcano plot
volcano_plot_pvalue_cs <- ggplot(deseq_res, aes(log2FoldChange, -log10(padj))) +
  geom_point(aes(color = sig, fill = sig), alpha = 0.6, size = 7, shape = 21, stroke = 0.7) +
  scale_color_manual(values = c("Higher in with-CHX" = '#006FBC', "Not Significant" = '#BBBBBC', "Higher in CHX-free" = '#F32327')) +
  scale_fill_manual(values = c("Higher in with-CHX" = '#006FBC', "Not Significant" = '#BBBBBC', "Higher in CHX-free" = '#F32327')) +
  theme(
    panel.grid = element_blank(), 
    panel.background = element_rect(color = 'black', fill = 'transparent'),
    legend.title = element_blank(), 
    legend.key = element_rect(fill = 'transparent', color = 'transparent'), # Remove box outline
    legend.background = element_rect(fill = 'transparent'),
    legend.position = "bottom", # Move legend to the bottom
    plot.title = element_text(size = 20,  family = "Arial", hjust = 0.5),
    axis.title.x = element_text(size = 15, family = "Arial"),
    axis.title.y = element_text(size = 15, family = "Arial"),
    axis.text.x = element_text(size = 15, family = "Arial"),
    axis.text.y = element_text(size = 15, family = "Arial")
  ) +
  geom_vline(xintercept = c(-1, 1), color = 'grey', size = 0.7, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), color = 'grey', size = 0.7, linetype = "dashed") +
  geom_text_repel(aes(label = name), 
                  size = 4,
                  family = "Arial",       # Set font family
                  box.padding = 0.5,      # Space around text
                  max.overlaps = 20) +
  labs(title = "Conventional soil", 
       x = "log2FoldChange", 
       y = "-log10(padj)")

# Display the volcano plot
print(volcano_plot_pvalue_cs)


########################Organic soil####################################################
setwd(dir = 'C:/Users/User/OneDrive/Desktop/qian/R/CYC-filter0.01/18dna-os')
groupdata <- read.table("group(-cyc28 vs. +cyc28).txt", sep = "\t", header = T)

level5 <- read.delim("level-5.txt",check.names = FALSE, row.names = 1)%>%
  as.data.frame(header=NULL,check.names=T) %>%
  select("31","32","33","34","35","36","37","38","39","40")
level5 <- rownames_to_column(level5, var = "Row_Name")

level5$Row_Name <- gsub("d5__", "", level5$Row_Name)
level5$Row_Name <- gsub("d3__", "", level5$Row_Name)
level5$Row_Name <- gsub("f__un_d2__", "d2_", level5$Row_Name)
level5$Row_Name <- gsub("f__un_d3__", "d3_", level5$Row_Name)
level5$Row_Name <- gsub("f__un_d4__", "d4_", level5$Row_Name)
level5$Row_Name <- gsub("f__un_", "Un_", level5$Row_Name)
level5$Row_Name <- gsub("d4_Centroplasthelida_X", "d4_Centroplasthelida", level5$Row_Name)

level5<- column_to_rownames(level5, var = "Row_Name")

# Step 1: Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = level5,colData = groupdata,design = ~group)

# Step 2: Run differential abundance analysis
dds <- DESeq(dds) 
suppressMessages(dds) 

# Step 3: Extract results
res <- results(dds,contrast = c('group','OSn-28','OSc-28'))
summary(res) 

# Step 4: Organize results
deseq_res <- as.data.frame(res[order(res$padj), ]) 
deseq_res$otu_id <- rownames(deseq_res)
#write.table(deseq_res[c(7, 1:6)], 'C:/Users/User/OneDrive/Desktop/qian/R/CYC-filter0.01/18dna-os/deseq_pic/DESeq2-level5(+cyc28.vs.-cyc28).txt', row.names = FALSE, sep = '\t', quote = FALSE)

# Step 5: Classify significance
for (i in 1:nrow(deseq_res)) { 
  if 
  (deseq_res[i,'log2FoldChange'] 
   %in% NA | abs(deseq_res[i,'log2FoldChange']) < 0.95) 
    deseq_res[i,'select_change'] <- 'n' else deseq_res[i,'select_change'] <- 'y'   
    if 
    (deseq_res[i,'padj'] 
     %in% NA | abs(deseq_res[i,'padj']) >= 0.05) 
      deseq_res[i,'select_padj'] <- 'n' else deseq_res[i,'select_padj'] <- 'y' 
      deseq_res[i,'select'] <- paste(deseq_res[i,'select_change'], deseq_res[i,'select_padj'], sep 
                                     = '') 
} 
deseq_res$sig <-"Not Significant"
deseq_res$sig[(deseq_res$log2FoldChange >= 0.95)&(deseq_res$padj<0.05)]<-"Higher in CHX-free"
deseq_res$sig[(deseq_res$log2FoldChange < (-0.95))&(deseq_res$padj<0.05)]<-"Higher in with-CHX"
deseq_res$name <-deseq_res$otu_id
deseq_res$name[(deseq_res$sig=="Not Significant")]<-NA

# Step 6: Volcano plot
volcano_plot_pvalue_os <- ggplot(deseq_res, aes(log2FoldChange, -log10(padj))) +
  geom_point(aes(color = sig, fill = sig), alpha = 0.6, size = 7, shape = 21, stroke = 0.7) +
  scale_color_manual(values = c("Higher in with-CHX" = '#006FBC', "Not Significant" = '#BBBBBC', "Higher in CHX-free" = '#F32327')) +
  scale_fill_manual(values = c("Higher in with-CHX" = '#006FBC', "Not Significant" = '#BBBBBC', "Higher in CHX-free" = '#F32327')) +
  theme(
    panel.grid = element_blank(), 
    panel.background = element_rect(color = 'black', fill = 'transparent'),
    legend.title = element_blank(), 
    legend.key = element_rect(fill = 'transparent', color = 'transparent'), # Remove box outline
    legend.background = element_rect(fill = 'transparent'),
    legend.position = "bottom", # Move legend to the bottom
    plot.title = element_text(size = 20,  family = "Arial", hjust = 0.5),
    axis.title.x = element_text(size = 15, family = "Arial"),
    axis.title.y = element_text(size = 15, family = "Arial"),
    axis.text.x = element_text(size = 15, family = "Arial"),
    axis.text.y = element_text(size = 15, family = "Arial")
  ) +
  geom_vline(xintercept = c(-1, 1), color = 'grey', size = 0.6, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), color = 'grey', size = 0.6, linetype = "dashed") +
  geom_text_repel(aes(label = name), 
                  size = 4,
                  family = "Arial",       # Set font family
                  box.padding = 0.5,      # Space around text
                  max.overlaps = 20) +
  labs(title = "Organic soil", 
       x = "log2FoldChange", 
       y = "-log10(padj)")

# Display the volcano plot
print(volcano_plot_pvalue_os)


############Code to Merge Two Boxplots into One Figure#####################
library(ggpubr)
# Arrange plots side by side
merged_plot <- ggarrange(volcano_plot_pvalue_os + theme(
                           plot.margin = margin(10, 10, 10, 20),  # Increase left margin
                           axis.title.y = element_text(margin = margin(r = 35))  # Add space between title and tick labels
                            ), 
                         volcano_plot_pvalue_cs+ theme(
                            plot.margin = margin(10, 10, 10, 20),  # Increase left margin
                            axis.title.y = element_text(margin = margin(r = 35))  # Add space between title and tick labels
                            ), 
                         labels = c("A", "B"),  # Add labels for identification
                         font.label = list(size = 20), 
                         ncol = 2, nrow = 1, 
                         common.legend = TRUE,
                         legend = "bottom", # Arrange in one row, two columns
                         align = "hv")  # Align horizontally and vertically

# Print the final plot
print(merged_plot)

# Save the merged plot as an image
ggsave("C:/Users/User/OneDrive/Desktop/qian/R/paper1_figs/merged_deseq.jpg", plot = merged_plot, device = "jpeg", 
       dpi = 600, width = 25, height = 12, units = "cm")

ggsave("C:/Users/User/OneDrive/Desktop/qian/R/paper1_figs/merged_deseq.tiff", plot = merged_plot, device = "tiff", 
       dpi = 600, width = 25, height = 12, units = "cm")