# Load required libraries
library(ggplot2)
library(dplyr)
library(ggpubr)
library(gridExtra)
library(grid)

# Clear the workspace
rm(list = ls())

################################# Organic Soil #####################################

# Set working directory and load data
setwd("C:/Users/User/OneDrive/Desktop/qian/R")
data <- read.table("Soil Predation data for R-OS.txt", sep = "\t", header = TRUE)

# Standardize treatment labels
data$Treatment <- gsub("CHX-free", "-CHX", data$Treatment)
data$Treatment <- gsub("with CHX", "+CHX", data$Treatment)

# Prepare labels for line ends
label_data <- data %>%
  group_by(Treatment) %>%
  filter(Time == max(Time)) %>%
  mutate(
    label_y = Mean + 0.45,
    label_x = Time - 3.5
  )

# Create line plot for organic soil
os <- ggplot(data, aes(x = Time, y = Mean, color = Treatment, group = Treatment)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2, size = 0.8) +
  geom_text(
    data = label_data, 
    aes(x = Time, y = label_y, label = Treatment), 
    size = 3.5, family = "Arial", color = "black", vjust = -0.35
  ) +
  labs(
    title = "Organic soil",
    x = "Time (Days)",
    y = expression(atop(italic("E. coli")~"concentration", "(log10 CFU/g soil)"))
  ) +
  scale_color_manual(values = c("-CHX" = "#9AC9DB", "+CHX" = "#40588F")) +
  coord_cartesian(ylim = c(3, 7), xlim = c(0, 32)) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    legend.position = "none",
    axis.text = element_text(size = 15, color = "black", family = "Arial"),
    axis.title = element_text(size = 15, color = "black", family = "Arial"),
    plot.title = element_text(size = 20, hjust = 0.5, family = "Arial")
  )

print(os)

################################# Conventional Soil #####################################

# Load data
data <- read.table("Soil Predation data for R-CS.txt", sep = "\t", header = TRUE)

# Standardize treatment labels
data$Treatment <- gsub("CHX-free", "-CHX", data$Treatment)
data$Treatment <- gsub("with CHX", "+CHX", data$Treatment)

# Prepare labels
label_data <- data %>%
  group_by(Treatment) %>%
  filter(Time == max(Time)) %>%
  mutate(
    label_y = Mean + 0.45,
    label_x = Time - 3.5
  )

# Create line plot for conventional soil
cs <- ggplot(data, aes(x = Time, y = Mean, color = Treatment, group = Treatment)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2, size = 0.8) +
  geom_text(
    data = label_data, 
    aes(x = Time, y = label_y, label = Treatment), 
    size = 3.5, family = "Arial", color = "black", vjust = -0.35
  ) +
  labs(
    title = "Conventional soil",
    x = "Time (Days)",
    y = expression(atop(italic("E. coli")~"concentration", "(log10 CFU/g soil)"))
  ) +
  scale_color_manual(values = c("-CHX" = "#9AC9DB", "+CHX" = "#40588F")) +
  coord_cartesian(ylim = c(3, 7), xlim = c(0, 32)) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    legend.position = "none",
    axis.text = element_text(size = 15, color = "black", family = "Arial"),
    axis.title = element_text(size = 15, color = "black", family = "Arial"),
    plot.title = element_text(size = 20, hjust = 0.5, family = "Arial")
  )

print(cs)

#################### Combine and Export Plots ####################

merged_plot <- ggarrange(
  os +
    theme(
      axis.title.y = element_text(size = 15, family = "Arial", color = "black", margin = margin(r = 15)),
      plot.margin = margin(5, -5, 5, 5)
    ),
  cs +
    theme(
      axis.title.y = element_blank(),
      plot.margin = margin(5, 5, 5, -5)
    ),
  labels = c("A", "B"),
  font.label = list(size = 20),
  ncol = 2,
  nrow = 1,
  align = "hv"
)

print(merged_plot)

# Save as high-resolution TIFF and JPEG
ggsave("merged_ecolinumberplots3.tiff",
       plot = merged_plot, dpi = 600, width = 25, height = 12, units = "cm")

ggsave("merged_ecolinumberplots3.jpg",
       plot = merged_plot, dpi = 600, width = 25, height = 12, units = "cm")
