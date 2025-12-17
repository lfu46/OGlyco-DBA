# Load required packages
library(tidyverse)
library(writexl)

# Function to process NGlyco results
process_NGlyco <- function(file_path) {
  # Import NGlyco HCD searching result
  NGlyco_HCD_result <- read_tsv(
    file_path,
    col_names = TRUE,
    name_repair = 'universal'
  ) |> 
    select(
      Spectrum, Peptide, Charge, Observed.M.Z, Delta.Mass, Hyperscore, 
      Assigned.Modifications, Total.Glycan.Composition, Glycan.Score, Glycan.q.value,
      Is.Decoy, Is.Contaminant, Protein.ID, Entry.Name, Gene, Protein.Description
    ) |> 
    filter(
      str_detect(Entry.Name, 'HUMAN')
    ) |> 
    filter(
      Is.Decoy == FALSE & Is.Contaminant == FALSE
    )
  
  # Calculate specificity and the number of unique glycopeptide and glycoprotein
  Num_total_psm <- nrow(NGlyco_HCD_result)
  
  NGlyco_HCD_result_glycoPSM <- NGlyco_HCD_result |> 
    filter(
      !is.na(Total.Glycan.Composition)
    )
  
  Num_glyco_psm <- nrow(NGlyco_HCD_result_glycoPSM)
  
  NGlyco_specificity <- Num_glyco_psm/Num_total_psm
  
  NGlyco_unique_glycopeptide <- NGlyco_HCD_result_glycoPSM |> 
    select(Peptide, Total.Glycan.Composition) |> 
    distinct()
  
  Num_unique_glycopeptide <- nrow(NGlyco_unique_glycopeptide)
  
  NGlyco_unique_glycoprotein <- NGlyco_HCD_result_glycoPSM |> 
    distinct(Protein.ID)
  
  Num_unique_glycoprotein <- nrow(NGlyco_unique_glycoprotein)
  
  # Apply glycan FDR filter with 0.05
  NGlyco_HCD_result_glycoPSM_FDR005 <- NGlyco_HCD_result_glycoPSM |> 
    filter(Glycan.q.value <= 0.05)
  
  Num_glyco_psm_FDR005 <- nrow(NGlyco_HCD_result_glycoPSM_FDR005)
  
  NGlyco_unique_glycopeptide_FDR005 <- NGlyco_HCD_result_glycoPSM_FDR005 |> 
    select(Peptide, Total.Glycan.Composition) |> 
    distinct()
  
  Num_unique_glycopeptide_FDR005 <- nrow(NGlyco_unique_glycopeptide_FDR005)
  
  NGlyco_unique_glycoprotein_FDR005 <- NGlyco_HCD_result_glycoPSM_FDR005 |> 
    distinct(Protein.ID)
  
  Num_unique_glycoprotein_FDR005 <- nrow(NGlyco_unique_glycoprotein_FDR005)
  
  # Return results as a named list
  return(list(
    Num_total_psm = Num_total_psm,
    Num_glyco_psm = Num_glyco_psm,
    NGlyco_specificity = NGlyco_specificity,
    Num_unique_glycopeptide = Num_unique_glycopeptide,
    Num_unique_glycoprotein = Num_unique_glycoprotein,
    Num_glyco_psm_FDR005 = Num_glyco_psm_FDR005,
    Num_unique_glycopeptide_FDR005 = Num_unique_glycopeptide_FDR005,
    Num_unique_glycoprotein_FDR005 = Num_unique_glycoprotein_FDR005
  ))
}

# Function to process OGlyco results
process_OGlyco <- function(file_path) {
  # Import OGlyco HCD searching result
  OGlyco_HCD_result <- read_tsv(
    file_path,
    col_names = TRUE,
    name_repair = 'universal'
  ) |> 
    select(
      Spectrum, Peptide, Charge, Observed.M.Z, Delta.Mass, Hyperscore, 
      Assigned.Modifications, Total.Glycan.Composition, Glycan.Score, Glycan.q.value,
      Is.Decoy, Is.Contaminant, Protein.ID, Entry.Name, Gene, Protein.Description
    ) |> 
    filter(
      str_detect(Entry.Name, 'HUMAN')
    ) |> 
    filter(
      Is.Decoy == FALSE & Is.Contaminant == FALSE
    )
  
  # Calculate specificity and the number of unique glycopeptide and glycoprotein
  Num_total_psm <- nrow(OGlyco_HCD_result)
  
  OGlyco_HCD_result_glycoPSM <- OGlyco_HCD_result |> 
    filter(
      !is.na(Total.Glycan.Composition)
    )
  
  Num_glyco_psm <- nrow(OGlyco_HCD_result_glycoPSM)
  
  OGlyco_specificity <- Num_glyco_psm/Num_total_psm
  
  OGlyco_unique_glycopeptide <- OGlyco_HCD_result_glycoPSM |> 
    select(Peptide, Total.Glycan.Composition) |> 
    distinct()
  
  Num_unique_glycopeptide <- nrow(OGlyco_unique_glycopeptide)
  
  OGlyco_unique_glycoprotein <- OGlyco_HCD_result_glycoPSM |> 
    distinct(Protein.ID)
  
  Num_unique_glycoprotein <- nrow(OGlyco_unique_glycoprotein)
  
  # Apply glycan FDR filter with 0.05
  OGlyco_HCD_result_glycoPSM_FDR005 <- OGlyco_HCD_result_glycoPSM |> 
    filter(Glycan.q.value <= 0.05)
  
  Num_glyco_psm_FDR005 <- nrow(OGlyco_HCD_result_glycoPSM_FDR005)
  
  OGlyco_unique_glycopeptide_FDR005 <- OGlyco_HCD_result_glycoPSM_FDR005 |> 
    select(Peptide, Total.Glycan.Composition) |> 
    distinct()
  
  Num_unique_glycopeptide_FDR005 <- nrow(OGlyco_unique_glycopeptide_FDR005)
  
  OGlyco_unique_glycoprotein_FDR005 <- OGlyco_HCD_result_glycoPSM_FDR005 |> 
    distinct(Protein.ID)
  
  Num_unique_glycoprotein_FDR005 <- nrow(OGlyco_unique_glycoprotein_FDR005)
  
  # Return results as a named list
  return(list(
    Num_total_psm = Num_total_psm,
    Num_glyco_psm = Num_glyco_psm,
    OGlyco_specificity = OGlyco_specificity,
    Num_unique_glycopeptide = Num_unique_glycopeptide,
    Num_unique_glycoprotein = Num_unique_glycoprotein,
    Num_glyco_psm_FDR005 = Num_glyco_psm_FDR005,
    Num_unique_glycopeptide_FDR005 = Num_unique_glycopeptide_FDR005,
    Num_unique_glycoprotein_FDR005 = Num_unique_glycoprotein_FDR005
  ))
}

# Read file paths for NGlyco
NGlyco_file_paths <- read_lines('/Volumes/cos-lab-rwu60/Longping/OGlyco_DBA_12152025/Optimized_Extraction/NGlyco_HCD_Search_1/NGlyco_file_path.txt')

# Process all NGlyco files
NGlyco_results_list <- list()
for (i in seq_along(NGlyco_file_paths)) {
  cat(sprintf("Processing NGlyco file %d of %d: %s\n", i, length(NGlyco_file_paths), 
              basename(dirname(NGlyco_file_paths[i]))))
  NGlyco_results_list[[i]] <- process_NGlyco(NGlyco_file_paths[i])
}

# Create NGlyco comparison table
NGlyco_comparison <- tibble(
  Sample = paste0("NGlyco_", 1:6),
  Num_total_psm = sapply(NGlyco_results_list, function(x) x$Num_total_psm),
  Num_glyco_psm = sapply(NGlyco_results_list, function(x) x$Num_glyco_psm),
  NGlyco_specificity = sapply(NGlyco_results_list, function(x) x$NGlyco_specificity),
  Num_unique_glycopeptide = sapply(NGlyco_results_list, function(x) x$Num_unique_glycopeptide),
  Num_unique_glycoprotein = sapply(NGlyco_results_list, function(x) x$Num_unique_glycoprotein),
  Num_glyco_psm_FDR005 = sapply(NGlyco_results_list, function(x) x$Num_glyco_psm_FDR005),
  Num_unique_glycopeptide_FDR005 = sapply(NGlyco_results_list, function(x) x$Num_unique_glycopeptide_FDR005),
  Num_unique_glycoprotein_FDR005 = sapply(NGlyco_results_list, function(x) x$Num_unique_glycoprotein_FDR005)
)

# Add summary statistics for NGlyco
NGlyco_summary <- NGlyco_comparison |>
  summarise(
    Sample = "Mean ± SD",
    Num_total_psm = sprintf("%.0f ± %.0f", mean(Num_total_psm), sd(Num_total_psm)),
    Num_glyco_psm = sprintf("%.0f ± %.0f", mean(Num_glyco_psm), sd(Num_glyco_psm)),
    NGlyco_specificity = sprintf("%.3f ± %.3f", mean(NGlyco_specificity), sd(NGlyco_specificity)),
    Num_unique_glycopeptide = sprintf("%.0f ± %.0f", mean(Num_unique_glycopeptide), sd(Num_unique_glycopeptide)),
    Num_unique_glycoprotein = sprintf("%.0f ± %.0f", mean(Num_unique_glycoprotein), sd(Num_unique_glycoprotein)),
    Num_glyco_psm_FDR005 = sprintf("%.0f ± %.0f", mean(Num_glyco_psm_FDR005), sd(Num_glyco_psm_FDR005)),
    Num_unique_glycopeptide_FDR005 = sprintf("%.0f ± %.0f", mean(Num_unique_glycopeptide_FDR005), sd(Num_unique_glycopeptide_FDR005)),
    Num_unique_glycoprotein_FDR005 = sprintf("%.0f ± %.0f", mean(Num_unique_glycoprotein_FDR005), sd(Num_unique_glycoprotein_FDR005))
  )



# Export NGlyco results to Excel with separate sheets
write_xlsx(
  list(
    "Individual_Results" = NGlyco_comparison,
    "Summary_Statistics" = NGlyco_summary
  ), 
  '/Volumes/cos-lab-rwu60/Longping/OGlyco_DBA_12152025/Optimized_Extraction/NGlyco_HCD_Search_1/NGlyco_Comparison_Results.xlsx'
)
cat("NGlyco results exported to NGlyco_Comparison_Results.xlsx\n\n")

# Read file paths for OGlyco
OGlyco_file_paths <- read_lines('/Volumes/cos-lab-rwu60/Longping/OGlyco_DBA_12152025/Optimized_Extraction/OGlyco_HCD_Search_1/OGlyco_file_path.txt')

# Process all OGlyco files
OGlyco_results_list <- list()
for (i in seq_along(OGlyco_file_paths)) {
  cat(sprintf("Processing OGlyco file %d of %d: %s\n", i, length(OGlyco_file_paths), 
              basename(dirname(OGlyco_file_paths[i]))))
  OGlyco_results_list[[i]] <- process_OGlyco(OGlyco_file_paths[i])
}

# Create OGlyco comparison table
OGlyco_comparison <- tibble(
  Sample = paste0("OGlyco_", 1:6),
  Num_total_psm = sapply(OGlyco_results_list, function(x) x$Num_total_psm),
  Num_glyco_psm = sapply(OGlyco_results_list, function(x) x$Num_glyco_psm),
  OGlyco_specificity = sapply(OGlyco_results_list, function(x) x$OGlyco_specificity),
  Num_unique_glycopeptide = sapply(OGlyco_results_list, function(x) x$Num_unique_glycopeptide),
  Num_unique_glycoprotein = sapply(OGlyco_results_list, function(x) x$Num_unique_glycoprotein),
  Num_glyco_psm_FDR005 = sapply(OGlyco_results_list, function(x) x$Num_glyco_psm_FDR005),
  Num_unique_glycopeptide_FDR005 = sapply(OGlyco_results_list, function(x) x$Num_unique_glycopeptide_FDR005),
  Num_unique_glycoprotein_FDR005 = sapply(OGlyco_results_list, function(x) x$Num_unique_glycoprotein_FDR005)
)

# Add summary statistics for OGlyco
OGlyco_summary <- OGlyco_comparison |>
  summarise(
    Sample = "Mean ± SD",
    Num_total_psm = sprintf("%.0f ± %.0f", mean(Num_total_psm), sd(Num_total_psm)),
    Num_glyco_psm = sprintf("%.0f ± %.0f", mean(Num_glyco_psm), sd(Num_glyco_psm)),
    OGlyco_specificity = sprintf("%.3f ± %.3f", mean(OGlyco_specificity), sd(OGlyco_specificity)),
    Num_unique_glycopeptide = sprintf("%.0f ± %.0f", mean(Num_unique_glycopeptide), sd(Num_unique_glycopeptide)),
    Num_unique_glycoprotein = sprintf("%.0f ± %.0f", mean(Num_unique_glycoprotein), sd(Num_unique_glycoprotein)),
    Num_glyco_psm_FDR005 = sprintf("%.0f ± %.0f", mean(Num_glyco_psm_FDR005), sd(Num_glyco_psm_FDR005)),
    Num_unique_glycopeptide_FDR005 = sprintf("%.0f ± %.0f", mean(Num_unique_glycopeptide_FDR005), sd(Num_unique_glycopeptide_FDR005)),
    Num_unique_glycoprotein_FDR005 = sprintf("%.0f ± %.0f", mean(Num_unique_glycoprotein_FDR005), sd(Num_unique_glycoprotein_FDR005))
  )

# Export OGlyco results to Excel with separate sheets
write_xlsx(
  list(
    "Individual_Results" = OGlyco_comparison,
    "Summary_Statistics" = OGlyco_summary
  ), 
  '/Volumes/cos-lab-rwu60/Longping/OGlyco_DBA_12152025/Optimized_Extraction/OGlyco_HCD_Search_1/OGlyco_Comparison_Results.xlsx'
)
cat("OGlyco results exported to OGlyco_Comparison_Results.xlsx\n\n")

# 
