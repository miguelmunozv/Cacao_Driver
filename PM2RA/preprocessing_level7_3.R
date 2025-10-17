##############################################################
#  Title: Taxonomic Data Filtering and Normalization Pipeline
#  Author: Miguel A. Muñoz Valencia
#  Description: This R script processes microbial abundance data
#               (level-7 taxonomy), filters redundant taxa using
#               correlation rules, and normalizes relative abundances
#               for downstream analysis (e.g., PM2RA).
#  Date: 2025-09-17
##############################################################

# --- 0. Setup and Library Loading ------------------------------------------

# Set working directory
setwd("~/Miguel/Investigación/Drivers_identification/data_bioprojects/PM2RA")

# Install required packages if not already installed
if (!require("writexl")) install.packages("writexl")
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("Hmisc")) install.packages("Hmisc")
if (!require("openxlsx")) install.packages("openxlsx")

# Load libraries
library(tidyverse)
library(writexl)
library(Hmisc)
library(openxlsx)

# --- 1. Import Data --------------------------------------------------------

df <- read_csv("C:/Users/migue/OneDrive/Documentos/Miguel/Investigación/Drivers_identification/data_bioprojects/PM2RA/level-7-final.csv")

# --- 2. Normalize Taxonomic Names -----------------------------------------

normalize_taxonomy <- function(taxon_name) {
  if (!grepl("^d__", taxon_name)) return(taxon_name)
  parts <- unlist(strsplit(taxon_name, ";"))
  parts <- trimws(parts)
  levels <- c("d__", "p__", "c__", "o__", "f__", "g__", "s__")
  new_parts <- rep(NA, length(levels))
  
  for (p in parts) {
    for (lvl in levels) {
      if (startsWith(p, lvl)) {
        new_parts[which(levels == lvl)] <- p
      }
    }
  }
  
  # Fill missing taxonomic ranks as "unknown"
  new_parts[is.na(new_parts)] <- paste0(levels[is.na(new_parts)], "unknown")
  paste(new_parts, collapse = ";")
}

# Apply normalization to column names
colnames(df) <- sapply(colnames(df), normalize_taxonomy, USE.NAMES = FALSE)

# --- 3. Separate Metadata and Abundance Tables ----------------------------

taxa_cols <- grep("^d__", colnames(df), value = TRUE)
df_abund <- df[, taxa_cols]
df_meta  <- df[, !(colnames(df) %in% taxa_cols)]

# --- 4. Filter Samples with Total Abundance > 10 --------------------------

row_sums <- rowSums(df_abund, na.rm = TRUE)
df_abund <- df_abund[row_sums > 10, ]
df_meta  <- df_meta[row_sums > 10, ]

# --- 5. Remove Taxa with Zero Total Abundance ------------------------------

col_sums <- colSums(df_abund, na.rm = TRUE)
df_abund <- df_abund[, col_sums > 0]

# Combine metadata + filtered abundance data
df_filtered <- cbind(df_meta, df_abund)

# ========================================================================
# --- 6. Correlation-Based Filtering with Rules --------------------------
# ========================================================================

# Compute Pearson correlation among taxa
cor_matrix <- cor(df_abund, method = "pearson")

# Identify highly correlated taxa (r > 0.99 but < 1)
high_cor_pairs <- which(abs(cor_matrix) > 0.99 & abs(cor_matrix) < 1, arr.ind = TRUE)
high_cor_pairs <- as.data.frame(high_cor_pairs)
high_cor_pairs <- high_cor_pairs[high_cor_pairs[, 1] < high_cor_pairs[, 2], ]

# Decision rules for correlated taxa
decide_taxon <- function(t1, t2) {
  parts1 <- unlist(strsplit(t1, ";"))
  parts2 <- unlist(strsplit(t2, ";"))
  
  # Rule 1: Prefer taxa with fewer "unknown" levels
  unknown_same <- any(grepl("unknown", parts1) & grepl("unknown", parts2))
  if (unknown_same) return(list(keep = t2, remove = t1, rule = "Rule 1"))
  
  # Rule 2: Keep the most complete taxonomy
  complete1 <- !any(grepl("unknown", parts1))
  complete2 <- !any(grepl("unknown", parts2))
  if (complete1 & !complete2) return(list(keep = t1, remove = t2, rule = "Rule 2"))
  if (complete2 & !complete1) return(list(keep = t2, remove = t1, rule = "Rule 2"))
  
  # Rule 3: Prefer taxa with genus-level assignment
  has_genus1 <- grepl("g__", t1) & !grepl("g__unknown", t1)
  has_genus2 <- grepl("g__", t2) & !grepl("g__unknown", t2)
  if (has_genus1 & !has_genus2) return(list(keep = t1, remove = t2, rule = "Rule 3"))
  if (has_genus2 & !has_genus1) return(list(keep = t2, remove = t1, rule = "Rule 3"))
  
  # Default: Keep the second taxon arbitrarily
  return(list(keep = t2, remove = t1, rule = "Default (tie)"))
}

# Apply rules to all correlated pairs
decisions <- list()
for (i in seq_len(nrow(high_cor_pairs))) {
  t1 <- colnames(df_abund)[high_cor_pairs[i, 1]]
  t2 <- colnames(df_abund)[high_cor_pairs[i, 2]]
  decision <- decide_taxon(t1, t2)
  
  decisions[[i]] <- data.frame(
    Taxon_1 = t1,
    Taxon_2 = t2,
    Keep = decision$keep,
    Remove = decision$remove,
    Rule = decision$rule,
    stringsAsFactors = FALSE
  )
}

# Save rule decisions
decisions_df <- bind_rows(decisions)
write.xlsx(decisions_df, file = "decisions_rules.xlsx", sheetName = "Rules", overwrite = TRUE)

# Remove redundant taxa
to_remove <- unique(decisions_df$Remove)
df_abund <- df_abund[, !(colnames(df_abund) %in% to_remove)]

# --- 7. Rebuild Final Dataset --------------------------------------------

df_final <- cbind(df_meta, df_abund)

# Remove "Unassigned" column if present
if ("Unassigned;__;__;__;__;__;__" %in% colnames(df_final)) {
  df_final$`Unassigned;__;__;__;__;__;__` <- NULL
}

# Save filtered data
write.xlsx(df_final, file = "level-7_filtered.xlsx")

# --- 8. Normalize by Sample (Relative Abundance) -------------------------

options(scipen = 999)  # Avoid scientific notation

# Separate numeric part
taxa_cols <- grep("^d__", colnames(df_final), value = TRUE)
df_numeric <- df_final[, taxa_cols]

# Normalize rows (sample-wise)
df_numeric_norm <- sweep(df_numeric, 1, rowSums(df_numeric), FUN = "/")

# Rebuild full dataset
df_final_norm <- cbind(df_final[, !(colnames(df_final) %in% taxa_cols)], df_numeric_norm)

# --- 9. Filter Taxa with < 30 Non-Zero Values ----------------------------

nonzero_counts <- colSums(df_numeric_norm > 0)
df_numeric_norm <- df_numeric_norm[, nonzero_counts >= 30]

df_final_norm <- cbind(df_final[, !(colnames(df_final) %in% taxa_cols)], df_numeric_norm)

# --- 10. Shorten Taxonomic Names to "Genus_Species" -----------------------

shorten_taxonomy <- function(taxon_name) {
  if (!grepl("^d__", taxon_name)) return(taxon_name)
  parts <- unlist(strsplit(taxon_name, ";"))
  genus <- sub("^g__", "", parts[grep("^g__", parts)])
  if (length(genus) == 0) genus <- "unknown"
  species <- sub("^s__", "", parts[grep("^s__", parts)])
  if (length(species) == 0) species <- "unknown"
  paste0(genus, "_", species)
}

colnames(df_final_norm) <- sapply(colnames(df_final_norm), shorten_taxonomy, USE.NAMES = FALSE)

# --- 11. Keep Only "Loc-name" + Taxa -------------------------------------

df_final_norm <- cbind(
  df_final_norm[, "Loc-name", drop = FALSE],
  df_final_norm[, !(colnames(df_final_norm) %in% "Loc-name")]
)

# Remove unnecessary metadata columns
df_subset_final <- df_final_norm %>%
  select(-c(index, BioProject, AvgSpotLen, Barcodes, `Center-Name`,
            `Collection-Date`, `Country-continent`, `Lat-lon`, ReleaseDate,
            HOST, `Env-broad`))

# --- 12. Check Missing Values and Variance -------------------------------

sum(is.na(df_subset_final))            # Count NA values
sum(is.nan(as.matrix(df_subset_final))) # Count NaN values
apply(df_subset_final, 2, var)         # Check column variance

# --- 13. Export Final Table ----------------------------------------------

write.csv(df_subset_final, "df_pm2ra.csv", row.names = FALSE)

#####################################################################
# End of Script
#####################################################################
