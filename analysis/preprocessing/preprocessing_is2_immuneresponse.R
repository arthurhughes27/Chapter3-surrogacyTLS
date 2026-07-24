# File for preprocessing the raw immune response data files

# Load packages
library(fs)
library(dplyr)
library(stringr)
library(janitor)
library(readxl)
library(purrr)
library(tidyr)

# Specify folder within folder root where the raw data lives
raw_data_folder = fs::path("data-raw", "is2")
processed_data_folder = fs::path("data", "is2")

# Use fs::path() to specify the data paths robustly
p_load_elisa <- fs::path(raw_data_folder, "elisa_2025-01-10_01-14-21.xlsx")
p_load_elispot <- fs::path(raw_data_folder, "elispot_2025-01-10_01-14.xlsx")
p_load_nAb <- fs::path(raw_data_folder, "neut_ab_titer_2025-01-10_01-13-22.xlsx")
p_load_hai <- fs::path(raw_data_folder, "hai_2025-01-10_01-13-41.xlsx")
p_load_clinical <- fs::path(processed_data_folder, "is2_clinical.rds")

# HAI response data
response_hai = read_excel(p_load_hai) %>%
  clean_names()

# nAb response data
response_nAb = read_excel(p_load_nAb) %>%
  clean_names()

# ELISA response data
response_elisa = read_excel(p_load_elisa) %>%
  clean_names()

# ELISPOT response data
response_elispot = read_excel(p_load_elispot) %>%
  clean_names()

# Clinical data for filtration of participants
is2_clinical = readRDS(p_load_clinical)

# First, filter each dataframe to only contain information on participants for which we have gene expression data
# Find identifiers of participants with gene expression measurements
participants = is2_clinical %>%
  pull(participant_id) %>%
  unique()

# Filter immune response data by these participants
response_hai = response_hai %>%
  filter(participant_id %in% participants)

response_nAb = response_nAb %>%
  filter(participant_id %in% participants)

response_elisa = response_elisa %>%
  filter(participant_id %in% participants)

response_elispot = response_elispot %>%
  filter(participant_id %in% participants)

# Depending on the assay, there may be multiple viral strains or analytes measured.
# We rename some columns and add a column for the assay name, so that we can merge these data together.

response_hai = response_hai %>%
  rename(response_strain_analyte = virus) %>%
  mutate(assay = "hai") %>% 
  dplyr::select(participant_id, gender, race, cohort, study_time_collected, study_time_collected_unit, response_strain_analyte, value_preferred, assay)

response_nAb = response_nAb %>%
  rename(response_strain_analyte = virus) %>%
  mutate(assay = "nAb") %>% 
  dplyr::select(participant_id, gender, race, cohort, study_time_collected, study_time_collected_unit, response_strain_analyte, value_preferred, assay)

# For ELISPOT: filter to take only IgG analytes, convert measure into normalised cell count per 100,000 cells
response_elispot = response_elispot %>%
  filter(grepl("IgG", analyte)) %>% 
  rename(response_strain_analyte = analyte) %>%
  mutate(assay = "elispot",
         unit_preferred = "analyte count per 100000") %>% 
  mutate(value_preferred = 100000*spot_number_reported/cell_number_preferred) %>% 
  filter(value_preferred < 100000,
         value_preferred > 0) %>% 
  dplyr::select(participant_id, gender, race, cohort, study_time_collected, study_time_collected_unit, response_strain_analyte, value_preferred, assay)

# Filter ELISA data for IgG and Hepatitis B antibodies
response_elisa = response_elisa %>%
  filter(grepl("IgG", analyte) |
           analyte == "Hepatitis B Virus Surface Antibody") %>%
  rename(response_strain_analyte = analyte) %>%
  mutate(assay = "elisa") %>% 
  dplyr::select(participant_id, gender, race, cohort, study_time_collected, study_time_collected_unit, response_strain_analyte, value_preferred, assay)

# Now merge all the raw response data
response_raw_merged = bind_rows(response_nAb, response_elisa, response_hai, response_elispot) %>%
  arrange(participant_id) %>%
  distinct()

# First get the studies and other clinical data corresponding to each participant id
is2_studies = is2_clinical %>%
  dplyr::select(participant_id,
         study_accession,
         vaccine_name) %>%
  distinct()

# Merge the study names into the immune response data (it is not directly given)
response_raw_merged_studies = merge(x = response_raw_merged,
                                    y = is2_studies,
                                    by = "participant_id",
                                    all = F)


# First, summarise across strains/analytes within each study, vaccine, assay, and time
response_summary <- response_raw_merged_studies %>%
  group_by(
    study_time_collected,
    assay,
    study_accession,
    vaccine_name,
    response_strain_analyte
  ) %>%
  mutate(
    value_standardised = {
      x <- value_preferred
      s <- sd(x, na.rm = TRUE)
      
      if (n() == 1) {
        x
      } else if (is.na(s) || s == 0) {
        x
      } else {
        as.numeric(scale(x))
      }
    }
  ) %>%
  ungroup() %>%
  group_by(
    participant_id,
    study_time_collected,
    assay,
    study_accession,
    vaccine_name
  ) %>%
  summarise(
    response_standardised_mean = mean(value_standardised, na.rm = TRUE),
    n_analytes = n_distinct(response_strain_analyte),
    .groups = "drop"
  )


# Now choose one assay per study based on availability, and in the order nAb > hai > elisa > elispot
assay_priority <- c("nAb", "hai", "elisa", "elispot")

# Choose the highest-priority assay available within each study
study_assay_choice <- response_summary %>%
  filter(!is.na(response_standardised_mean)) %>%
  distinct(study_accession, assay) %>%
  mutate(priority = match(assay, assay_priority)) %>%
  filter(!is.na(priority)) %>%
  group_by(study_accession) %>%
  summarise(
    chosen_assay = assay[which.min(priority)],
    .groups = "drop"
  )

# Keep only the chosen assay per study and pivot to wide format
response_wide <- response_summary %>%
  left_join(study_assay_choice, by = "study_accession") %>%
  filter(assay == chosen_assay) %>%
  dplyr::select(participant_id, study_accession, vaccine_name, study_time_collected, response_standardised_mean) %>%
  distinct() %>%
  mutate(study_time_collected = as.character(study_time_collected)) %>%
  pivot_wider(
    names_from = study_time_collected,
    values_from = response_standardised_mean,
    names_prefix = "ab_p_"
  )

# Identify ab_p_ columns
ab_cols <- grep("^ab_p_", names(response_wide), value = TRUE)

# Extract numeric time values and sort
ab_cols_sorted <- ab_cols[order(as.numeric(sub("^ab_p_", "", ab_cols)))]

# Reorder dataframe columns
is2_immResp <- response_wide %>%
  dplyr::select(
    participant_id,
    study_accession,
    vaccine_name,
    all_of(ab_cols_sorted)
  )

# Use fs::path() to specify the data path robustly
p_save <- fs::path(processed_data_folder, "is2_immResp.rds")

# Save dataframe
saveRDS(is2_immResp, file = p_save)

rm(list = ls())
