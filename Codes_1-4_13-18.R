### Section 4 CODE 1 ###

library(tidyverse)

# List all TSV files
tsv_files <- list.files(pattern = "\\.tsv$")

### Section 4 CODE 2 ###

process_data_list <- lapply(tsv_files, function(file) {
  # Process each TSV file
  df <- read_tsv(file) %>%
    mutate(
      GFAP = `Multiplex - GFAP_D` / `Multiplex - Cy5`,
      CD31 = `Multiplex - CD31` / `Multiplex - TRITC`,
      CD4 = `Multiplex - CD4 EPR6854` / `Multiplex - Cy5`,
      CD8 = `Multiplex - CD8` / `Multiplex - TRITC`,
      CD68 = `Multiplex - CD68 PG-M1` / `Multiplex - TRITC`,
      CD11c = `Multiplex - CD11c` / `Multiplex - Cy5`,
      CD163 = `Multiplex - CD163` / `Multiplex - Cy5`,
      CD205 = `Multiplex - CD205_D` / `Multiplex - Cy5`,
      TMEM119 = `Multiplex - TMEM119_C` / `Multiplex - Cy5`,
      P2RY12 = `Multiplex - P2RY12` / `Multiplex - Cy5`,
      FOXP3 = `Multiplex - FOXP3 236A/E6` / `Multiplex - TRITC`,
      pSTAT3 = `Multiplex - p-STAT3_C` / `Multiplex - Cy5`,
      PD1 = `Multiplex - PD-1 EPR4877(2)` / `Multiplex - TRITC`,
      PDL1 = `Multiplex - PD-L1` / `Multiplex - Cy5`,
      TIM3 = `Multiplex - TIM-3` / `Multiplex - TRITC`,
      LAG3 = `Multiplex - LAG-3` / `Multiplex - Cy5`,
      GRZMB = `Multiplex - GRZMB` / `Multiplex - TRITC`,
      LCK = `Multiplex - LCK` / `Multiplex - Cy5`,
      HLADR = `Multiplex - HLA-DR_B` / `Multiplex - TRITC`,
    ) %>%
    select(Area, Distance, GFAP, CD31, CD4, CD8, CD68, CD11c, CD163, CD205, TMEM119, P2RY12, FOXP3, pSTAT3, PD1, PDL1, TIM3, LAG3, GRZMB, LCK, HLADR) %>%
    filter(!is.na(Distance))
  
  processed_file_name <- sub("\\.tsv$", "_Processed.tsv", file)
  
  # Write the processed data to a new TSV file
  write_tsv(df, processed_file_name)
})

### Section 5 CODE 3 ###
library(tidyverse)

# List all TSV files
tsv_files <- list.files(pattern = "\\.tsv$")

### Section 5 CODE 4 ###

process_data_list <- lapply(tsv_files, function(file) {
  # Process each TSV file
  df <- read_tsv(file) %>%
    rename(LabelName = `ObjectInfo - LabelName`) %>%
    mutate(
      GFAP = `Multiplex - GFAP_D` / `Multiplex - Cy5`,
      CD31 = `Multiplex - CD31` / `Multiplex - TRITC`,
      CD4 = `Multiplex - CD4 EPR6854` / `Multiplex - Cy5`,
      CD8 = `Multiplex - CD8` / `Multiplex - TRITC`,
      CD68 = `Multiplex - CD68 PG-M1` / `Multiplex - TRITC`,
      CD11c = `Multiplex - CD11c` / `Multiplex - Cy5`,
      CD163 = `Multiplex - CD163` / `Multiplex - Cy5`,
      CD205 = `Multiplex - CD205_D` / `Multiplex - Cy5`,
      TMEM119 = `Multiplex - TMEM119_C` / `Multiplex - Cy5`,
      P2RY12 = `Multiplex - P2RY12` / `Multiplex - Cy5`,
      FOXP3 = `Multiplex - FOXP3 236A/E6` / `Multiplex - TRITC`,
      pSTAT3 = `Multiplex - p-STAT3_C` / `Multiplex - Cy5`,
      PD1 = `Multiplex - PD-1 EPR4877(2)` / `Multiplex - TRITC`,
      PDL1 = `Multiplex - PD-L1` / `Multiplex - Cy5`,
      TIM3 = `Multiplex - TIM-3` / `Multiplex - TRITC`,
      LAG3 = `Multiplex - LAG-3` / `Multiplex - Cy5`,
      GRZMB = `Multiplex - GRZMB` / `Multiplex - TRITC`,
      LCK = `Multiplex - LCK` / `Multiplex - Cy5`,
      HLADR = `Multiplex - HLA-DR_B` / `Multiplex - TRITC`,
    ) %>%
    select(LabelName, Area, GFAP, CD31, CD4, CD8, CD68, CD11c, CD163, CD205, TMEM119, P2RY12, FOXP3, pSTAT3, PD1, PDL1, TIM3, LAG3, GRZMB, LCK, HLADR) %>%
    filter(LabelName != "Cell")
  
  processed_file_name <- sub("\\.tsv$", "_Processed.tsv", file)
  
  # Write the processed data to a new TSV file
  write_tsv(df, processed_file_name)
})

### Section 7 CODE 13 ###

library(tidyverse)

# List all TSV files
tsv_files <- list.files(pattern = "\\.tsv$")

### Section 7 CODE 14 ###

# Vectors phenotypes of interests based for different analyses
all_phenotypes1 <- c("CD4+", "CD8+", "CD68+", "CD11c+", "CD163+", "CD205+", 
                     "TMEM119+", "P2RY12+", "FOXP3+", "p-STAT3+", "PD-1+", 
                     "TIM-3+","LAG-3+", "GRZMB+", "LCK+", "HLADR+", "Other")
all_phenotypes2 <- c("CD11c+CD163+", "CD11c+CD205+", "CD11c+TMEM119+",
                     "CD11c+CD68+", "CD163+CD68+", "CD11c+CD163-CD205-TMEM119-CD68-",
                     "CD4+LCK+", "CD8+LCK+", "CD4+p-STAT3+", "CD8+p-STAT3+",
                     "CD4+p-STAT3-", "CD8+p-STAT3-", "CD8+FOXP3+", "Other")

# Initialize the count matrices with zeros
Phenotype1_Counts <- matrix(0, nrow = length(tsv_files), ncol = length(all_phenotypes1), dimnames = list(NULL, all_phenotypes1))
Phenotype2_Counts <- matrix(0, nrow = length(tsv_files), ncol = length(all_phenotypes2), dimnames = list(NULL, all_phenotypes2))

### Section 7 CODE 15 ###

# Process each TSV file
for (idx in seq_along(tsv_files)) {
  file <- tsv_files[idx]
  df <- read_tsv(file)
  
  # Manually define phenotypes based on positivity values
  df <- df %>%
    mutate(Phenotype1 = case_when(
      CD4_Positive == 1 ~ "CD4+",
      CD8_Positive == 1 ~ "CD8+",
      CD68_Positive == 1 ~ "CD68+",
      CD11c_Positive == 1 ~ "CD11c+",
      CD163_Positive == 1 ~ "CD163+",
      CD205_Positive == 1 ~ "CD205+",
      TMEM119_Positive == 1 ~ "TMEM119+",
      P2RY12_Positive == 1 ~ "P2RY12+",
      FOXP3_Positive == 1 ~ "FOXP3+",
      pSTAT3_Positive == 1 ~ "p-STAT3+",
      PD1_Positive == 1 ~ "PD-1+",
      TIM3_Positive == 1 ~ "TIM-3+",
      LAG3_Positive == 1 ~ "LAG-3+",
      GRZMB_Positive == 1 ~ "GRZMB+",
      LCK_Positive == 1 ~ "LCK+",
      HLADR_Positive == 1 ~ "HLADR+",
      # Add more conditions here for other phenotypes
      TRUE ~ "Other"),  # Default case if no other conditions are met
      Phenotype2 = case_when(
        CD11c_Positive == 1 & CD163_Positive == 1 ~ "CD11c+CD163+",
        CD11c_Positive == 1 & CD205_Positive == 1 ~ "CD11c+CD205+",
        CD11c_Positive == 1 & TMEM119_Positive == 1 ~ "CD11c+TMEM119+",
        CD11c_Positive == 1 & CD68_Positive == 1 ~ "CD11c+CD68+",
        CD163_Positive == 1 & CD68_Positive == 1 ~ "CD163+CD68+",
        CD11c_Positive == 1 & CD163_Positive == 0 & CD205_Positive == 0 & 
          TMEM119_Positive == 0 & CD68_Positive == 0 ~ "CD11c+CD163-CD205-TMEM119-CD68-",
        CD4_Positive == 1 & LCK_Positive == 1 ~ "CD4+LCK+",
        CD8_Positive == 1 & LCK_Positive == 1 ~ "CD8+LCK+",
        CD4_Positive == 1 & pSTAT3_Positive == 1 ~ "CD4+p-STAT3+",
        CD8_Positive == 1 & pSTAT3_Positive == 1 ~ "CD8+p-STAT3+",
        CD4_Positive == 1 & pSTAT3_Positive == 0 ~ "CD4+p-STAT3-",
        CD8_Positive == 1 & pSTAT3_Positive == 0 ~ "CD8+p-STAT3-",
        CD8_Positive == 1 & FOXP3_Positive == 1 ~ "CD8+FOXP3+",
        # Add more conditions here for other phenotypes
        TRUE ~ "Other")
    ) %>%
    select(Distance, Phenotype1, Phenotype2)
  
  
  # Count the unique values for Phenotype1 and Phenotype2
  phenotype1_counts <- table(factor(df$Phenotype1, levels = all_phenotypes1))
  phenotype2_counts <- table(factor(df$Phenotype2, levels = all_phenotypes2))
  
  # Update the count matrices
  Phenotype1_Counts[idx, names(phenotype1_counts)] <- as.numeric(phenotype1_counts)
  Phenotype2_Counts[idx, names(phenotype2_counts)] <- as.numeric(phenotype2_counts)
  
  # Write the processed data to a new TSV file
  processed_file_name <- sub("\\.tsv$", "_Processed.tsv", file)
  
  # Write the processed data to a new TSV file
  write_tsv(df, processed_file_name)
}

# The count matrices are now populated with counts from each TSV file
# You can access them as Phenotype1_Counts and Phenotype2_Counts

### Section 7 CODE 16 ###

# List all processed TSV files
processed_files <- list.files(pattern = "_Processed\\.tsv$")

# Read and combine all processed files into one big dataframe
Combined_Processed_Data <- processed_files %>% 
  map_dfr(~read_tsv(.x))

# Write the combined dataframe to a new TSV file
write_tsv(Combined_Processed_Data , "Combined_Processed_Data.tsv")

### Section 7 CODE 17 ###

### CHART 1

# Write the combined dataframe to a new TSV file
write_tsv(Combined_Processed_Data , "Combined_Processed_Data.tsv")


# Read the TSV file
df_visual <- read_tsv("Combined_Processed_Data.tsv")

# Filter for 'CD11c+' cells
cd11c_positive <- df_visual %>% 
  filter(Phenotype1 == "CD11c+")

# Define the maximum distance for creating bins
max_distance <- max(cd11c_positive$Distance, na.rm = TRUE)

# Create distance bins
bin_width <- 50
bins <- seq(0, max_distance + bin_width, by = bin_width)

# Plot the histogram
ggplot(cd11c_positive, aes(x = Distance)) +
  geom_histogram(breaks = bins, color = "black", fill = "blue") +
  labs(title = "Histogram of CD11c+ Cells Count by Distance",
       x = "Distance",
       y = "Count of CD11c+ Cells") +
  theme_minimal()

### Section 7 CODE 18 ###

### CHART 2


# Read the TSV file
df <- read_tsv("Combined_Processed_Data.tsv")

# Define the maximum distance for creating bins
max_distance <- max(df$Distance, na.rm = TRUE)
bin_width <- 50
bins <- seq(0, max_distance + bin_width, by = bin_width)

# Bin the data
df_binned <- df %>%
  mutate(Bin = cut(Distance, breaks = bins, include.lowest = TRUE, right = FALSE)) %>%
  group_by(Bin, Phenotype1) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  # Calculate the total count per bin
  mutate(Total = sum(Count)) %>%
  # Calculate the percentage
  mutate(Percent = (Count / Total) * 100)

# Plot the stacked bar chart
ggplot(df_binned, aes(x = Bin, y = Percent, fill = Phenotype1)) +
  geom_bar(stat = "identity") +
  labs(title = "Percentage Composition of CD11c+ vs CD11c- Cells by Distance",
       x = "Distance Bin",
       y = "Percentage (%)") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))