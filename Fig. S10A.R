# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)
rm(list = ls())

####################### Step 1: Load and prepare relative abundance data #######################
setwd('C:/Users/User/OneDrive/Desktop/qian/R/cyc+-relative_abundance(lefse)/18dna-os')
RA_data <- read.table("level-8.txt", sep = "\t", header = TRUE, check.names = FALSE)

# Reshape and transform
RA_data <- RA_data %>%
  pivot_longer(cols = -group, names_to = "Sample", values_to = "RelativeAbundance") %>%
  mutate(
    RelativeAbundance = RelativeAbundance * 100,
    group = gsub("st__un_d[0-9]+__", "un_", group),
    group = gsub("d[0-9]+__", "", group),
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

####################### Step 2: Load and clean trophic data #######################
setwd('C:/Users/User/OneDrive/Desktop/qian/R/CYC-filter0.01/18dna-os')
trophicdata <- read.table("df_filtered-trophic.txt", sep = "\t", header = TRUE, check.names = FALSE)
colnames(trophicdata)[colnames(trophicdata) == "Rank8"] <- "name"
trophicdata$name <- gsub("st__un_d[0-9]+__", "un_", trophicdata$name)
trophicdata$name <- gsub("d[0-9]+__", "", trophicdata$name)
trophicdata_clean <- trophicdata %>% distinct(name, .keep_all = TRUE)
trophicdata_clean$name <- trimws(trophicdata_clean$name)

####################### Step 3: Load and clean network correlation data #######################
setwd('C:/Users/User/OneDrive/Desktop/qian/R/network/18os+cyc-ecoli')
ecoli_r <- read.csv("R.ecoli_level8_sig.csv", sep = ",", header = TRUE)
colnames(ecoli_r) <- trimws(colnames(ecoli_r))
ecoli_r <- ecoli_r %>% dplyr::select(name = X, rvalue = `E.coli`) %>% filter(name != "E.coli")
ecoli_r$name <- gsub("st__un_d[0-9]+__", "un_", ecoli_r$name)
ecoli_r$name <- gsub("d[0-9]+__", "", ecoli_r$name)
ecoli_r_less0 <- ecoli_r %>% filter(rvalue < 0)

####################### Step 4: Merge with trophic data and filter for Consumers #######################
merged_trophic_network_r_less0 <- ecoli_r_less0 %>%
  left_join(trophicdata_clean %>% select(name, Trophic), by = "name") %>%
  distinct() %>%
  filter(Trophic == "Consumer")

output_dir <- "C:/Users/User/OneDrive/Desktop/qian/R/paper1_figs/key protists RA change_scatterplot_18os/network+cyc_r_less0/Consumer"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

####################### Step 5: Filter RA tables by relative abundance >1% #######################
for (i in 1:nrow(merged_trophic_network_r_less0)) {
  assign(paste0("taxa_less0_", i), merged_trophic_network_r_less0$name[i])
  data <- RA_data %>% filter(group == get(paste0("taxa_less0_", i))) %>% select(-group, -Sample)
  if (max(data$RelativeAbundance, na.rm = TRUE) >= 1) {
    assign(paste0("taxa_less0_", i, "_data"), data)
  }
}
remaining_indices <- which(sapply(1:nrow(merged_trophic_network_r_less0), function(i) exists(paste0("taxa_less0_", i, "_data"))))

####################### Step 6: Add CHX rows for Day0 #######################
for (i in remaining_indices) {
  df <- get(paste0("taxa_less0_", i, "_data"))
  new_rows <- df %>% filter(Time == "Day0", Treatment == "-CHX") %>% mutate(Treatment = "+CHX")
  updated_df <- bind_rows(df, new_rows)
  updated_df$Time <- factor(updated_df$Time, levels = c("Day0", "Day14", "Day28"))
  assign(paste0("updated_taxa_less0_", i, "_data"), updated_df)
}

####################### Step 7: Plotting #######################
num_cols <- 3
plot_list <- list()

for (i in remaining_indices) {
  data_to_plot <- get(paste0("updated_taxa_less0_", i, "_data"))
  plot_title <- get(paste0("taxa_less0_", i))
  ylim_max <- max(data_to_plot$RelativeAbundance, na.rm = TRUE) * 1.1
  is_leftmost <- ((which(remaining_indices == i) - 1) %% num_cols) == 0
  
  p <- ggplot(data_to_plot, aes(x = Treatment, y = RelativeAbundance, color = Time)) +
    geom_boxplot(position = position_dodge(0.6), width = 0.5, size = 0.8, outlier.shape = NA) +
    geom_jitter(position = position_jitterdodge(0.15, 0.6), size = 3, alpha = 0.7) +
    stat_summary(fun = median, geom = "crossbar", width = 0.4, size = 0.4,
                 position = position_dodge(0.5), aes(color = Time), show.legend = FALSE) +
    labs(title = plot_title, x = "", y = ifelse(is_leftmost, "Relative abundance (%)", "")) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(0, ylim_max)) +
    scale_color_manual(values = c("Day0" = "#3A3737", "Day14" = "#1f77b4", "Day28" = "#ffcc00")) +
    theme_classic() +
    theme(
      axis.text = element_text(size = 20, color = "black", family = "Arial"),
      axis.title = element_text(size = 20, color = "black", family = "Arial"),
      plot.title = element_text(size = 20, hjust = 0.5, family = "Arial"),
      legend.title = element_blank(),
      legend.position = "bottom"
    )
  
  ggsave(filename = file.path(output_dir, paste0("taxa_less0_", i, "_", gsub("[^a-zA-Z0-9]", "_", plot_title), ".jpg")),
         plot = p, device = "jpg", width = 9, height = 8, units = "cm", dpi = 600)
  
  plot_list[[length(plot_list) + 1]] <- p
}

combined_plot <- wrap_plots(plotlist = plot_list, ncol = num_cols) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

combined_plot <- combined_plot + plot_annotation(
  title = "Organic soil",
  theme = theme(plot.title = element_text(size = 30, hjust = 0.5, family = "Arial"))
)

print(combined_plot)

ggsave(file.path(output_dir, "combined_taxa_plots_r_less0-box.jpg"),
       plot = combined_plot, device = "jpg", width = 28, height = 22, units = "cm", dpi = 600)
