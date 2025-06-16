# Load required libraries
library(tidyverse)
library(dplyr)
library(tibble)
library(reshape2)
library(stringr)
library(rstatix)     # For normality tests
library(ggpubr)      # For plot combination and annotation
library(FSA)         # For Dunn's test
library(multcompView)# For significance letters

# ---------------------------
# Process: Conventional Soil
# ---------------------------
rm(list = ls())
setwd("C:/Users/User/OneDrive/Desktop/qian/R/CYC-filter0.01/18dna-cs")

# Import metadata and OTU table
data_meta <- read.table("group.txt", sep = "\t", header = TRUE)
data_otu <- read.table("featuretable-rarefied.txt", sep = "\t", header = TRUE, check.names = FALSE, row.names = 1) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample")

# Merge OTU and metadata
merged_data <- merge(data_otu, data_meta, by = "sample", all.x = TRUE) %>%
  select(-c(sample, group)) %>%
  relocate(Replicate, .before = everything()) %>%
  t() %>%
  as.data.frame()

# Extract and assign column names
colnames(merged_data) <- merged_data[1, ]
merged_data <- merged_data[-1, ]

# Convert all values to numeric and transpose
merged_data <- merged_data %>%
  mutate_all(as.numeric) %>%
  t() %>%
  as.data.frame()

# Filter out ASVs that are all 0 or all 1
merged_data <- merged_data[, colSums(merged_data == 0) != nrow(merged_data)]
merged_data <- merged_data[, colSums(merged_data == 1) != nrow(merged_data)]

# Binarize: Convert all non-zero entries to 1
merged_data[merged_data != 0] <- 1

# Sum ASVs across samples
asv_sum <- colSums(merged_data) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("Sample")

# Convert long format for plotting
asv_long <- asv_sum %>%
  pivot_longer(cols = -Sample, names_to = "Category", values_to = "Value") %>%
  mutate(
    Treatment = case_when(str_starts(Category, "CSc") ~ "+CHX", str_starts(Category, "CSn") ~ "-CHX"),
    Time = case_when(
      str_detect(Category, "-0-") ~ "Day0",
      str_detect(Category, "-14-") ~ "Day14",
      str_detect(Category, "-28-") ~ "Day28"
    ),
    Replicate = str_extract(Category, "[1-5]$")
  ) %>%
  select(Treatment, Time, Replicate, Value)

# Shapiro–Wilk normality test
asv_long %>% group_by(Treatment, Time) %>% summarise(p_value = shapiro.test(Value)$p.value, .groups = "drop")

# Compute summary statistics
summary_stats <- asv_long %>%
  group_by(Treatment, Time) %>%
  summarise(
    Mean = mean(Value),
    SE = sd(Value)/sqrt(n()),
    .groups = "drop"
  )

# Manually add missing group (+CHX Day0) if needed
summary_stats <- bind_rows(summary_stats, summary_stats %>% filter(Time == "Day0") %>% mutate(Treatment = "+CHX"))

# Set factor levels
summary_stats$Time <- factor(summary_stats$Time, levels = c("Day0", "Day14", "Day28"))
summary_stats$Treatment <- factor(summary_stats$Treatment, levels = c("-CHX", "+CHX"))

# Dunn's test for pairwise comparison
asv_long$Group <- interaction(asv_long$Treatment, asv_long$Time)
dunn_result <- dunnTest(Value ~ Group, data = asv_long, method = "bonferroni")
dunn_df <- as.data.frame(dunn_result$res)
colnames(dunn_df) <- c("Comparison", "Z", "p.unadj", "p.adj")

# Convert comparisons into matrix for multcompView
comp_split <- strsplit(dunn_df$Comparison, " - ")
group1 <- sapply(comp_split, `[`, 1)
group2 <- sapply(comp_split, `[`, 2)
groups <- unique(c(group1, group2))
p_mat <- matrix(1, nrow = length(groups), ncol = length(groups), dimnames = list(groups, groups))
for (i in seq_along(group1)) {
  p_mat[group1[i], group2[i]] <- dunn_df$p.adj[i]
  p_mat[group2[i], group1[i]] <- dunn_df$p.adj[i]
}

# Assign letters
letters_df <- multcompLetters(p_mat, threshold = 0.05)$Letters %>%
  enframe(name = "Group", value = "label")

# Merge letters into summary
summary_stats$Group <- interaction(summary_stats$Treatment, summary_stats$Time)
summary_stats <- left_join(summary_stats, letters_df, by = "Group") %>%
  mutate(label_y = Mean + SE + 5)

# Plot bar chart
bar_cs <- ggplot(summary_stats, aes(x = Treatment, y = Mean, fill = Time)) +
  geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.6) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2, position = position_dodge(0.7)) +
  geom_text(aes(y = label_y, label = label), position = position_dodge(0.7), vjust = 0, size = 5, color = "black", family = "Arial") +
  scale_fill_manual(values = c("Day0" = "#737373", "Day14" = "#1f77b4", "Day28" = "#ffcc00")) +
  ylim(0, 200) +
  labs(title = "Conventional Soil", x = "", y = "Observed ASVs") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 15, color = "black", family = "Arial"),
    axis.text.y = element_text(size = 15, color = "black", family = "Arial"),
    axis.title.y = element_text(size = 15, color = "black", family = "Arial"),
    plot.title = element_text(size = 20, hjust = 0.5, family = "Arial"),
    legend.title = element_blank()
  )

# Save output
ggsave("C:/Users/User/OneDrive/Desktop/qian/R/paper1_figs/Observed_ASVs_CS.tiff", bar_cs, dpi = 600, width = 16, height = 12, units = "cm")






# ---------------------------
# Process: Organic Soil
# ---------------------------
setwd(dir = 'C:/Users/User/OneDrive/Desktop/qian/R/CYC-filter0.01/18dna-os')

# Import metadata and OTU table
metadata18 <- read.table("group.txt", sep = "\t", header = T) %>%dplyr::select(-c(X,X.1,X.2,X.3,X.4,X.5,X.6))
otu18<- read.table("featuretable-rarefied.txt", sep = "\t", header = T,check.names = FALSE,row.names = 1) %>% t()%>% as.data.frame()%>%rownames_to_column("sample")#Import data file

# Merge OTU and metadata
mergedtable18 <- merge(otu18, metadata18, by = "sample", all.x = TRUE)%>%dplyr::select(-c(sample,group))%>%relocate(Replicate, .before = everything())#%>% t()%>% as.data.frame()
mergedtable18 <- t(mergedtable18)%>% as.data.frame()
column_names <- as.character(mergedtable18[1, ])
mergedtable18 <- mergedtable18[-1, ]

# Extract and assign column names
colnames(mergedtable18) <- column_names

mergedtableall <- mergedtable18 %>%
  mutate_all(as.numeric)%>% t()%>% as.data.frame() 

# Filter out ASVs that are all 0 
zero_columns <- apply(mergedtableall, 2, function(col) all(col == 0))#筛选全为0的列
mat_zeros <-mergedtableall[, zero_columns]
mergedtableall <- mergedtableall[, colSums(mergedtableall == 0) != nrow(mergedtableall)]#%>% t()
# Filter out ASVs that are all 1
one_columns <- apply(mergedtableall, 2, function(col) all(col == 1))#筛选全为1的列
mat_ones <-mergedtableall[, one_columns]
mergedtableall <- mergedtableall[, colSums(mergedtableall == 1) != nrow(mergedtableall)]%>% t()%>% as.data.frame()

# Binarize: Convert all non-zero entries to 1
mergedtableall[mergedtableall != 0] <- 1

# Calculate column sums
col_sums <- colSums(mergedtableall)

# Append as last row
df <- rbind(mergedtableall, col_sums)

# Optionally rename the last row
rownames(df)[nrow(df)] <- "Sum"

#To keep only the last row of a dataframe and output it in R
ASV_num <- tail(df, 1)

# Convert row names into a column
ASV_num <- ASV_num %>% 
  tibble::rownames_to_column("Sample")

ASV_num_long <- ASV_num %>%
  pivot_longer(cols = -Sample, names_to = "Category", values_to = "Value")

# Add Treatment, Time and Replicate columns based on Category values
ASV_num_long <- ASV_num_long %>%
  mutate(
    # Assign Treatment based on the prefix of "Category"
    Treatment = case_when(
      str_starts(Category, "OSc") ~ "+CHX",
      str_starts(Category, "OSn") ~ "-CHX"
    ),
    
    # Extract Time based on the middle number in "Category"
    Time = case_when(
      str_detect(Category, "-0-")  ~ "Day0",
      str_detect(Category, "-14-") ~ "Day14",
      str_detect(Category, "-28-") ~ "Day28"
    ),
    
    # Extract Replicate based on the last number in "Category"
    Replicate = str_extract(Category, "[1-5]$")
  )%>%
  dplyr::select(c(Treatment,Time,Replicate,Value))

##Check normality (optional, but good practice)(If any group has p < 0.05, proceed with Kruskal–Wallis + Dunn’s test.)
# Shapiro test for each group (Treatment × Time)
ASV_num_long %>%
  group_by(Treatment, Time) %>%
  summarise(p_value = shapiro.test(Value)$p.value, .groups = "drop")

# Compute mean and standard error
summary_ASV <- ASV_num_long %>%
  group_by(Treatment, Time) %>%
  summarise(
    Mean = mean(Value),
    SE = sd(Value) / sqrt(n()),  # Standard error
    .groups = "drop"
  )


new_row <- summary_ASV %>%
  filter(Time == "Day0") %>% #Selects the row with "Day0"
  mutate(Treatment = "+CHX") #Changes "Treatment" from NA to "CHX"

# Bind the new row to the original dataframe
summary_ASV <- bind_rows(summary_ASV, new_row) #Appends the new row to the original dataframe.

# Convert 'Time' column to factor with specified levels
summary_ASV$Time <- factor(summary_ASV$Time, levels = c("Day0", "Day14", "Day28"))

# Modify treatment order: Move "CHX" before "NA"
summary_ASV$Treatment <- factor(summary_ASV$Treatment, levels = c("-CHX", "+CHX"))

# Perform Dunn’s test with Bonferroni correction
ASV_num_long$Group <- interaction(ASV_num_long$Treatment, ASV_num_long$Time)
dunn_result <- dunnTest(Value ~ Group, data = ASV_num_long, method = "bonferroni")

dunn_df <- as.data.frame(dunn_result$res)
colnames(dunn_df) <- c("Comparison", "Z", "p.unadj", "p.adj")

# Build pairwise matrix for multcompLetters
comparison_split <- strsplit(dunn_df$Comparison, " - ")
group1 <- sapply(comparison_split, `[`, 1)
group2 <- sapply(comparison_split, `[`, 2)

# Unique groups
all_groups <- unique(c(group1, group2))
p_mat <- matrix(1, nrow = length(all_groups), ncol = length(all_groups),
                dimnames = list(all_groups, all_groups))

for (i in seq_along(group1)) {
  p_val <- dunn_df$p.adj[i]
  p_mat[group1[i], group2[i]] <- p_val
  p_mat[group2[i], group1[i]] <- p_val
}

# Generate significance letters
letters_obj <- multcompLetters(p_mat, threshold = 0.05)
group_labels <- data.frame(
  Group = names(letters_obj$Letters),
  label = letters_obj$Letters
)

# Merge letters into summary_ASV
summary_ASV$Group <- interaction(summary_ASV$Treatment, summary_ASV$Time)
summary_ASV <- merge(summary_ASV, group_labels, by = "Group", all.x = TRUE)

# Set position for label display
summary_ASV <- summary_ASV %>%
  mutate(label_y = Mean + SE + 5)

# Plot with error bars
bar_os <- ggplot(summary_ASV, aes(x = Treatment, y = Mean, fill = Time)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +  # Bar plot
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), 
                width = 0.2, 
                position = position_dodge(width = 0.7)) +  # Error bars
  geom_text(aes(y = label_y, label = label),
            position = position_dodge(width = 0.7),
            vjust = 0, size = 5, color = "black", family = "Arial") +
  labs(title = "Organic Soil", x = "", y = "Observed ASVs") +  # Labels
  scale_fill_manual(values = c("Day0" = "#737373", "Day14" = "#1f77b4", "Day28" = "#ffcc00")) +  # Custom colors
  ylim(0, 200)+
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 15, color = "black", family = "Arial"),
    axis.text.y = element_text(size = 15, color = "black", family = "Arial"),
    axis.title.y = element_text(size = 15, color = "black", family = "Arial"),
    plot.title = element_text(size = 20, hjust = 0.5, family = "Arial"),
    legend.title = element_blank()
  )

print(bar_os)

# Save the merged plot as an image
ggsave("C:/Users/User/OneDrive/Desktop/qian/R/paper1_figs/Observed_ASVs_OS.tiff", plot = bar_os, device = "tiff",
       dpi = 600, width = 16, height = 12, units = "cm")


############Code to Merge Two Boxplots into One Figure#####################
library(ggpubr)
# Arrange plots side by side
merged_plot <- ggarrange(bar_os + 
                           theme(
                             #axis.title.x = element_blank(),
                             axis.title.y = element_text(size = 15, family = "Arial", color = "black",margin = margin(r = 30)),
                             plot.margin = margin(5, -5, 5, 5)  # Remove all margins
                           ),
                         bar_cs+
                           theme(
                             axis.title.y = element_blank(),
                             #axis.text.y = element_blank(),
                             #axis.ticks.y = element_blank(),
                             #axis.title.x = element_blank(),
                             plot.margin = margin(5, 5, 5, -5)  # Remove all margins
                           ), 
                         labels = c("A", "B"),  # Add labels for identification
                         font.label = list(size = 20), 
                         ncol = 2, nrow = 1, 
                         common.legend = TRUE,
                         legend = "bottom", # Arrange in one row, two columns
                         align = "hv")  # Align horizontally and vertically

print(merged_plot)

# Save the merged plot as an image
ggsave("C:/Users/User/OneDrive/Desktop/qian/R/paper1_figs/merged_Observed_ASVs.jpg", plot = merged_plot, device = "jpeg",
       dpi = 600, width = 25, height = 12, units = "cm")

ggsave("C:/Users/User/OneDrive/Desktop/qian/R/paper1_figs/merged_Observed_ASVs.tiff", plot = merged_plot, device = "tiff",
       dpi = 600, width = 25, height = 12, units = "cm")


