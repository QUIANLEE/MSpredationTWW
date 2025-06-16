# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)

# Clear the workspace
rm(list = ls())

###################### Organic soil without CHX ######################
# Set working directory
setwd("C:/Users/User/OneDrive/Desktop/qian/R/cyc+-relative_abundance(lefse)/18dna-os")

# Load relative abundance data
RA_data <- read.table("level-8.txt", sep = "\t", header = TRUE, check.names = FALSE)

# Reshape to long format
RA_data <- RA_data %>%
  pivot_longer(
    cols = -group,
    names_to = "Sample",
    values_to = "RelativeAbundance"
  ) %>%
  mutate(
    RelativeAbundance = RelativeAbundance * 100,  # Convert to percentage
    group = gsub("st__un_d[0-9]+__", "un_", group),
    group = gsub("d8__", "", group)
  ) %>%
  mutate(
    Treatment = case_when(
      grepl("^OSc-", Sample) ~ "+CHX",
      grepl("^OSn-", Sample) ~ "-CHX"
    ),
    Time = case_when(
      grepl("-14$", Sample) ~ "Day14",
      grepl("-28$", Sample) ~ "Day28",
      grepl("-0$", Sample)  ~ "Day0"
    )
  )

# Load trophic data
setwd("C:/Users/User/OneDrive/Desktop/qian/R/CYC-filter0.01/18dna-os")
trophicdata <- read.table("df_filtered-trophic.txt", sep = "\t", header = TRUE, check.names = FALSE) %>%
  rename(name = Rank8) %>%
  mutate(
    name = gsub("st__un_d[0-9]+__", "un_", name),
    name = gsub("d8__", "", name),
    name = trimws(name)
  ) %>%
  distinct(name, .keep_all = TRUE)

# Load DESeq2 results and filter
setwd("C:/Users/User/OneDrive/Desktop/qian/R/CYC-filter0.01/18dna-os/deseq_pic")
deseqdata <- read.table("DESeq2-level8(+cyc28.vs.-cyc28).txt", sep = "\t", header = TRUE, check.names = FALSE) %>%
  filter(sig == "Higher in CHX-free")

# Merge with trophic assignments
merged_trophic_deseq <- deseqdata %>%
  left_join(select(trophicdata, name, Trophic), by = "name") %>%
  distinct() %>%
  filter(Trophic == "Consumer")

# Set output directory
output_dir <- "C:/Users/User/OneDrive/Desktop/qian/R/paper1_figs/key protists RA change_scatterplot_18os/Deseq-Higher in CHX-free/consumer"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Extract and filter taxa
for (i in 1:nrow(merged_trophic_deseq)) {
  assign(paste0("taxa", i), merged_trophic_deseq$name[i])
  data <- RA_data %>% filter(group == get(paste0("taxa", i))) %>% select(-group, -Sample)
  if (max(data$RelativeAbundance, na.rm = TRUE) >= 1) {
    assign(paste0("taxa", i, "_data"), data)
  }
}

remaining_taxa_indices <- which(sapply(1:nrow(merged_trophic_deseq), function(i) exists(paste0("taxa", i, "_data"))))

# Augment data by adding Day0 +CHX
for (i in remaining_taxa_indices) {
  original_data <- get(paste0("taxa", i, "_data"))
  new_rows <- original_data %>% filter(Time == "Day0", Treatment == "-CHX") %>% mutate(Treatment = "+CHX")
  assign(paste0("updated_taxa", i, "_data"), bind_rows(original_data, new_rows))
}

# Convert Time to ordered factor
for (i in remaining_taxa_indices) {
  df <- get(paste0("updated_taxa", i, "_data"))
  df$Time <- factor(df$Time, levels = c("Day0", "Day14", "Day28"))
  assign(paste0("updated_taxa", i, "_data"), df)
}

# Plotting setup
num_cols <- 3
plot_list <- list()

for (i in remaining_taxa_indices) {
  data_to_plot <- get(paste0("updated_taxa", i, "_data"))
  plot_title <- get(paste0("taxa", i))
  ylim_max <- max(data_to_plot$RelativeAbundance, na.rm = TRUE) * 1.1
  is_leftmost <- ((which(remaining_taxa_indices == i) - 1) %% num_cols) == 0
  
  plot_for_saving <- ggplot(data_to_plot, aes(x = Treatment, y = RelativeAbundance, color = Time)) +
    geom_boxplot(position = position_dodge(width = 0.6), width = 0.5, size = 0.8, outlier.shape = NA) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6), size = 3, alpha = 0.7) +
    stat_summary(fun = median, geom = "crossbar", width = 0.4, size = 0.4,
                 position = position_dodge(width = 0.5), aes(color = Time), show.legend = FALSE) +
    labs(title = plot_title, x = "", y = "Relative abundance (%)") +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(0, ylim_max)) +
    scale_color_manual(values = c("Day0" = "#3A3737", "Day14" = "#1f77b4", "Day28" = "#ffcc00")) +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 20, color = "black", family = "Arial"),
      axis.text.y = element_text(size = 20, color = "black", family = "Arial"),
      axis.title.y = element_text(size = 20, color = "black", family = "Arial"),
      plot.title = element_text(size = 20, hjust = 0.5, family = "Arial"),
      legend.title = element_blank(),
      legend.position = "bottom"
    )
  
  ggsave(paste0(output_dir, "/taxa", i, "_", gsub("[^a-zA-Z0-9]", "_", plot_title), ".jpg"),
         plot = plot_for_saving, device = "jpg", width = 8, height = 10, units = "cm", dpi = 600)
  
  plot_for_merging <- plot_for_saving + labs(y = ifelse(is_leftmost, "Relative abundance (%)", ""))
  plot_list[[length(plot_list) + 1]] <- plot_for_merging
}

# Combine all plots
combined_plot <- wrap_plots(plotlist = plot_list, ncol = num_cols) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

combined_plot <- combined_plot +
  plot_annotation(title = "Organic soil", theme = theme(plot.title = element_text(size = 30, hjust = 0.5, family = "Arial")))

print(combined_plot)

# Save combined plot
ggsave(file.path(output_dir, "combined_taxa_plots-box.jpg"), plot = combined_plot,
       device = "jpg", width = 28, height = 22, units = "cm", dpi = 600)