# ============================================================
# Title: Metadata preprocessing for 16S metagenome analysis
# # Date: 2025-10-16
# Description:
#   This script merges and filters metadata from multiple
#   BioProjects, performs downsampling of sequencing reads,
#   imputes invalid values, and formats metadata for QIIME2.
# ============================================================


# ============================================================
# 0. Install and load required packages
# ============================================================

if (!require("tidyverse")) install.packages("tidyverse")  # includes ggplot2, dplyr, readr
if (!require("readxl")) install.packages("readxl")

library(tidyverse)
library(dplyr)
library(ggplot2)
library(readxl)

cat("R version:", R.version.string, "\n")
cat("dplyr version:", packageVersion("dplyr"), "\n")


# ============================================================
# 1. Set working directory and import metadata
# ============================================================

setwd("~/Miguel/Investigación/Drivers_identification/data_bioprojects/bioproject_metadata")

# Load general metadata
df <- read.csv("SraRunTable_metagenome.csv")

# Load run information for each BioProject
df_spots1 <- read_csv("PRJNA1015591_runinfo.csv")
df_spots2 <- read_csv("PRJNA747886_runinfo.csv")
df_spots3 <- read_csv("PRJNA755178_runinfo.csv")
df_spots4 <- read_csv("PRJNA786222_runinfo.csv")
df_spots5 <- read_csv("PRJNA898075_runinfo.csv")

# Merge all run info into one table
df_spots <- bind_rows(df_spots1, df_spots2, df_spots3, df_spots4, df_spots5)


# ============================================================
# 2. Merge metadata and run info by "Run"
# ============================================================

df$Run <- as.character(df$Run)
df_spots$Run <- as.character(df_spots$Run)

df_merged <- df %>%
  left_join(df_spots, by = "Run")

head(df_merged)


# ============================================================
# 3. Filter metadata by sequencing type and BioProjects
# ============================================================

df_filtered <- df_merged %>%
  filter(
    GENE == "16S",
    BioProject.x %in% c("PRJNA1015591", "PRJNA747886",
                        "PRJNA755178", "PRJNA786222", "PRJNA898075"),
    toupper(LibraryLayout.x) == "PAIRED"
  )

# Export filtered run list
write.table(
  df_filtered$Run,
  file = "runs_list.txt",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

# Export filtered metadata
write.csv(df_filtered, "df_filtered.csv", row.names = FALSE)


# ============================================================
# 4. Summary and visualization of BioProjects
# ============================================================

# Number of samples per BioProject
df_filtered %>%
  group_by(BioProject.x) %>%
  summarise(num_samples = n_distinct(Run))

# Summarize sequencing statistics
df_summary <- df_filtered %>%
  group_by(BioProject.x) %>%
  summarise(
    total_bases   = sum(Bases, na.rm = TRUE),
    total_reads   = sum(spots, na.rm = TRUE),
    avg_read_len  = mean(AvgSpotLen, na.rm = TRUE),
    num_samples   = n_distinct(Run)
  )

# --- Visualization ---
# Total bases per BioProject
ggplot(df_summary, aes(x = BioProject.x, y = total_bases)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Total sequenced bases per BioProject",
       x = "BioProject", y = "Total bases") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Total reads per BioProject
ggplot(df_summary, aes(x = BioProject.x, y = total_reads)) +
  geom_bar(stat = "identity", fill = "seagreen") +
  labs(title = "Total reads per BioProject",
       x = "BioProject", y = "Number of reads") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Average read length
ggplot(df_summary, aes(x = BioProject.x, y = avg_read_len)) +
  geom_bar(stat = "identity", fill = "darkorange") +
  labs(title = "Average read length per BioProject",
       x = "BioProject", y = "Average length (bp)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Number of samples
ggplot(df_summary, aes(x = BioProject.x, y = num_samples)) +
  geom_bar(stat = "identity", fill = "purple") +
  labs(title = "Number of samples per BioProject",
       x = "BioProject", y = "Number of samples") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# ============================================================
# 5. Downsampling reads per sample (~50,000 reads)
# ============================================================

# Compute mean and SD of reads
stats <- df_filtered %>%
  group_by(BioProject.x) %>%
  summarise(
    mean_spots = mean(spots, na.rm = TRUE),
    sd_spots   = sd(spots, na.rm = TRUE),
    .groups = "drop"
  )

# Downsample preserving relative SD
df_downsampled <- df_filtered %>%
  left_join(stats, by = "BioProject.x") %>%
  mutate(
    spots_downsampled = round(((spots - mean_spots) / sd_spots) * sd_spots + 50000)
  ) %>%
  select(-mean_spots, -sd_spots)

# Check new distributions
df_downsampled %>%
  group_by(BioProject.x) %>%
  summarise(
    mean_new = mean(spots_downsampled, na.rm = TRUE),
    sd_new   = sd(spots_downsampled, na.rm = TRUE),
    .groups = "drop"
  )

# Visualization after downsampling
ggplot(df_downsampled, aes(x = BioProject.x, y = spots_downsampled, fill = BioProject.x)) +
  geom_boxplot(fill = "grey") +
  theme_minimal() +
  labs(
    title = "Read count distribution after downsampling",
    x = "BioProject",
    y = "Reads per sample"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# ============================================================
# 6. Impute negative values with global mean
# ============================================================

mean_valid <- mean(df_downsampled$spots_downsampled[df_downsampled$spots_downsampled > 0], na.rm = TRUE)

df_downsampled <- df_downsampled %>%
  mutate(spots_downsampled = ifelse(spots_downsampled < 0, round(mean_valid), spots_downsampled))

write.csv(df_downsampled, "df_downsampled.csv", row.names = FALSE)


# ============================================================
# 7. Prepare metadata for QIIME2
# ============================================================

df_qiime <- df_downsampled %>%
  select(Run, BioProject.x, AvgSpotLen, CenterName, Collection_Date,
         geo_loc_name_country_continent, geo_loc_name_country,
         geo_loc_name, ReleaseDate.x, HOST, env_broad_scale) %>%
  rename(
    `sample-id`         = Run,
    BioProject          = BioProject.x,
    `Center-Name`       = CenterName,
    `Collection-Date`   = Collection_Date,
    `Country-continent` = geo_loc_name_country_continent,
    `Loc-name`          = geo_loc_name_country,
    `Lat-lon`           = geo_loc_name,
    ReleaseDate         = ReleaseDate.x,
    `Env-broad`         = env_broad_scale
  )

# Replace missing or invalid HOST values
df_qiime <- df_qiime %>%
  mutate(HOST = ifelse(HOST == "" | HOST == "not applicable", "Theobroma cacao", HOST))

# Adjust environment categories
df_qiime <- df_qiime %>%
  mutate(
    `Env-broad` = case_when(
      `Env-broad` %in% c("ENVO_01000763", "Mineral", "cacao plantation") ~ "Cacao soil",
      `Loc-name` %in% c("Indonesia", "Bolivia") & `Env-broad` %in% c("", NA) ~ "Cacao soil",
      `Env-broad` == "tropical lowland evergreen broadleaf rain forest" ~ "Other soil",
      `Loc-name` == "USA" ~ "Greenhouse",
      TRUE ~ `Env-broad`
    )
  )


# ============================================================
# 8. Add random barcodes for QIIME2
# ============================================================

set.seed(123)

random_dna <- function(n, length = 12) {
  bases <- c("A", "T", "C", "G")
  replicate(n, paste0(sample(bases, length, replace = TRUE), collapse = ""))
}

df_qiime$Barcodes <- random_dna(nrow(df_qiime), length = 12)

# Reorder columns: place Barcodes after AvgSpotLen
df_qiime <- df_qiime %>%
  relocate(Barcodes, .after = AvgSpotLen)


# ============================================================
# 9. Export final metadata file
# ============================================================

write.table(
  df_qiime,
  file = "metadata_metagenomePE16S.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("Metadata processing completed successfully ✅\n")
