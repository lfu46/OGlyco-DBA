# import packages
library(tidyverse)

# import OGlyco HCD searching result
OGlyco_HCD_result <- read_tsv(
  '/Volumes/cos-lab-rwu60/Longping/OGlyco_DBA_12152025/Standard_Extraction/OGlyco_HCD_Search_1/E_LF_OGlyco_DBA_1_12152025/psm.tsv',
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

# calculate specificity and the number of unique glycopeptide and glycoprotein
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

# apply glycan FDR filter with 0.05
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


