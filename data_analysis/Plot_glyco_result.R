# Load required packages
library(tidyverse)
library(readxl)
library(patchwork)

# Read NGlyco results
NGlyco_data <- read_excel('/Volumes/cos-lab-rwu60/Longping/OGlyco_DBA_12152025/Optimized_Extraction/NGlyco_HCD_Search_1/NGlyco_Comparison_Results.xlsx', 
                          sheet = "Individual_Results")

# Read OGlyco results
OGlyco_data <- read_excel('/Volumes/cos-lab-rwu60/Longping/OGlyco_DBA_12152025/Optimized_Extraction/OGlyco_HCD_Search_1/OGlyco_Comparison_Results_1.xlsx', 
                          sheet = "Individual_Results")

# Add condition labels to NGlyco data
NGlyco_data_long <- NGlyco_data |>
  mutate(Condition = case_when(
    Sample %in% c("NGlyco_1", "NGlyco_2", "NGlyco_3") ~ "96-well filter",
    Sample %in% c("NGlyco_4", "NGlyco_5", "NGlyco_6") ~ "Centrifuge"
  )) |>
  select(-Sample)

# Add condition labels to OGlyco data
OGlyco_data_long <- OGlyco_data |>
  mutate(Condition = case_when(
    Sample %in% c("OGlyco_1", "OGlyco_2", "OGlyco_3") ~ "96-well filter",
    Sample %in% c("OGlyco_4", "OGlyco_5", "OGlyco_6") ~ "Centrifuge"
  )) |>
  select(-Sample)

# Function to create bar plot with overlay points and error bars
create_barplot <- function(data, metrics, metric_labels, title, filename) {
  
  # Reshape data for plotting
  plot_data <- data |>
    select(Condition, all_of(metrics)) |>
    pivot_longer(cols = all_of(metrics), 
                 names_to = "Metric", 
                 values_to = "Value")
  
  # Calculate summary statistics
  summary_stats <- plot_data |>
    group_by(Condition, Metric) |>
    summarise(
      Mean = mean(Value),
      SD = sd(Value),
      .groups = "drop"
    )
  
  # Create factor with custom labels
  plot_data$Metric <- factor(plot_data$Metric, 
                             levels = metrics, 
                             labels = metric_labels)
  summary_stats$Metric <- factor(summary_stats$Metric, 
                                 levels = metrics, 
                                 labels = metric_labels)
  
  # Create plot
  p <- ggplot() +
    # Bar plot for means
    geom_bar(data = summary_stats, 
             aes(x = Metric, y = Mean, fill = Condition),
             stat = "identity", position = position_dodge(width = 0.9),
             width = 0.8) +
    # Error bars
    geom_errorbar(data = summary_stats,
                  aes(x = Metric, ymin = Mean - SD, ymax = Mean + SD, group = Condition),
                  position = position_dodge(width = 0.9),
                  width = 0.25, linewidth = 0.5) +
    # Individual points
    geom_point(data = plot_data,
               aes(x = Metric, y = Value, group = Condition),
               position = position_dodge(width = 0.9),
               size = 2, shape = 21, stroke = 1) +
    # Mean value labels
    geom_text(data = summary_stats,
              aes(x = Metric, y = Mean + SD, label = sprintf("%.0f", Mean), group = Condition),
              position = position_dodge(width = 0.9),
              vjust = -0.5, size = 3.5) +
    scale_fill_manual(values = c("96-well filter" = "#4DBBD5FF", 
                                 "Centrifuge" = "#E64B35FF")) +
    labs(title = title,
         x = NULL,
         y = "No. of glyco identifications",
         fill = "") +
    theme_bw() +
    theme(
      plot.title = element_text(size = 9, hjust = 0.5, color = 'black', family = 'Arial'),
      axis.text.x = element_text(angle = 30, hjust = 1, size = 8, color = 'black', family = 'Arial'),
      axis.text.y = element_text(size = 8, color = 'black', family = 'Arial'),
      axis.title.y = element_text(size = 8, color = 'black', family = 'Arial'),
      legend.position = "right",
      legend.title = element_text(size = 8, color = 'black', family = 'Arial'),
      legend.text = element_text(size = 8),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
  
  # Save plot
  ggsave(filename, plot = p, width = 5, height = 3, dpi = 300, device = "tiff")
  cat(sprintf("Saved: %s\n", filename))
  
  return(p)
}

# Function to create specificity bar plot (percentage format)
create_specificity_barplot <- function(data, metric, title, filename) {
  
  # Reshape data for plotting - convert to percentage
  plot_data <- data |>
    select(Condition, all_of(metric)) |>
    rename(Value = all_of(metric)) |>
    mutate(Value = Value * 100)  # Convert to percentage
  
  # Calculate summary statistics
  summary_stats <- plot_data |>
    group_by(Condition) |>
    summarise(
      Mean = mean(Value),
      SD = sd(Value),
      .groups = "drop"
    )
  
  # Create plot
  p <- ggplot() +
    # Bar plot for means
    geom_bar(data = summary_stats, 
             aes(x = Condition, y = Mean, fill = Condition),
             stat = "identity", position = position_dodge(width = 0.9),
             width = 0.7) +
    # Error bars
    geom_errorbar(data = summary_stats,
                  aes(x = Condition, ymin = Mean - SD, ymax = Mean + SD),
                  width = 0.25, linewidth = 0.5) +
    # Individual points
    geom_point(data = plot_data,
               aes(x = Condition, y = Value),
               size = 2, shape = 21, stroke = 1) +
    # Mean value labels
    geom_text(data = summary_stats,
              aes(x = Condition, y = Mean + SD, label = sprintf("%.2f", Mean)),
              vjust = -0.5, size = 3.5) +
    scale_fill_manual(values = c("96-well filter" = "#4DBBD5FF", 
                                 "Centrifuge" = "#E64B35FF")) +
    labs(title = title,
         x = NULL,
         y = "Specificity (%)",
         fill = "") +
    theme_bw() +
    theme(
      plot.title = element_text(size = 9, hjust = 0.5, color = 'black', family = 'Arial'),
      axis.text.x = element_text(angle = 30, hjust = 1, size = 8, color = 'black', family = 'Arial'),
      axis.text.y = element_text(size = 8, color = 'black', family = 'Arial'),
      axis.title.y = element_text(size = 8, color = 'black', family = 'Arial'),
      legend.position = "right",
      legend.title = element_text(size = 8, color = 'black', family = 'Arial'),
      legend.text = element_text(size = 8),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
  
  # Save plot
  ggsave(filename, plot = p, width = 3, height = 3, dpi = 300, device = "tiff")
  cat(sprintf("Saved: %s\n", filename))
  
  return(p)
}

# Generate NGlyco plots (without FDR - counts only)
cat("Generating NGlyco plots (without FDR - counts)...\n")
create_barplot(
  data = NGlyco_data_long,
  metrics = c("Num_glyco_psm", "Num_unique_glycopeptide", "Num_unique_glycoprotein"),
  metric_labels = c("Total GlycoPSMs", "Unique Glycopeptides", "Unique Glycoproteins"),
  title = "N-Glyco Results",
  filename = "/Volumes/cos-lab-rwu60/Longping/OGlyco_DBA_12152025/Optimized_Extraction/NGlyco_HCD_Search_1/NGlyco_Comparison_noFDR.tif"
)

# Generate NGlyco specificity plot
cat("Generating NGlyco specificity plot...\n")
create_specificity_barplot(
  data = NGlyco_data_long,
  metric = "NGlyco_specificity",
  title = "N-Glyco Results (Specificity)",
  filename = "/Volumes/cos-lab-rwu60/Longping/OGlyco_DBA_12152025/Optimized_Extraction/NGlyco_HCD_Search_1/NGlyco_Specificity.tif"
)

# Generate NGlyco plots (with FDR)
cat("Generating NGlyco plots (with FDR)...\n")
create_barplot(
  data = NGlyco_data_long,
  metrics = c("Num_glyco_psm_FDR005", "Num_unique_glycopeptide_FDR005", 
              "Num_unique_glycoprotein_FDR005"),
  metric_labels = c("Total GlycoPSMs", 
                    "Unique Glycopeptides", "Unique Glycoproteins"),
  title = "N-Glyco Results (Glycan FDR: 0.05)",
  filename = "/Volumes/cos-lab-rwu60/Longping/OGlyco_DBA_12152025/Optimized_Extraction/NGlyco_HCD_Search_1/NGlyco_Comparison_FDR005.tif"
)

# Generate OGlyco plots (without FDR - counts only)
cat("Generating OGlyco plots (without FDR - counts)...\n")
create_barplot(
  data = OGlyco_data_long,
  metrics = c("Num_glyco_psm", "Num_unique_glycopeptide", "Num_unique_glycoprotein"),
  metric_labels = c("Total GlycoPSMs", "Unique Glycopeptides", "Unique Glycoproteins"),
  title = "O-Glyco Results",
  filename = "/Volumes/cos-lab-rwu60/Longping/OGlyco_DBA_12152025/Optimized_Extraction/OGlyco_HCD_Search_1/OGlyco_Comparison_noFDR.tif"
)

# Generate OGlyco specificity plot
cat("Generating OGlyco specificity plot...\n")
create_specificity_barplot(
  data = OGlyco_data_long,
  metric = "OGlyco_specificity",
  title = "O-Glyco Results (Specificity)",
  filename = "/Volumes/cos-lab-rwu60/Longping/OGlyco_DBA_12152025/Optimized_Extraction/OGlyco_HCD_Search_1/OGlyco_Specificity.tif"
)

# Generate OGlyco plots (with FDR)
cat("Generating OGlyco plots (with FDR)...\n")
create_barplot(
  data = OGlyco_data_long,
  metrics = c("Num_glyco_psm_FDR005", "Num_unique_glycopeptide_FDR005", 
              "Num_unique_glycoprotein_FDR005"),
  metric_labels = c("Total GlycoPSMs", 
                    "Unique Glycopeptides", "Unique Glycoproteins"),
  title = "O-Glyco Results (Glycan FDR: 0.05)",
  filename = "/Volumes/cos-lab-rwu60/Longping/OGlyco_DBA_12152025/Optimized_Extraction/OGlyco_HCD_Search_1/OGlyco_Comparison_FDR005.tif"
)

cat("\nAll plots generated successfully!\n")
cat("Output files:\n")
cat("  - NGlyco_Comparison_noFDR.tif (count metrics)\n")
cat("  - NGlyco_Specificity.tif (specificity in %)\n")
cat("  - NGlyco_Comparison_FDR005.tif\n")
cat("  - OGlyco_Comparison_noFDR.tif (count metrics)\n")
cat("  - OGlyco_Specificity.tif (specificity in %)\n")
cat("  - OGlyco_Comparison_FDR005.tif\n")
