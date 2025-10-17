install.packages("tidyverse")  # incluye ggplot2, dplyr, readr

R.version.string


# Cargar librería
library(dplyr)
library(tidyverse)
library(ggplot2)
library(readxl)

packageVersion("dplyr")

setwd("~/Miguel/Investigación/Drivers_identification/data_bioprojects/bioproject_metadata")

# Leer el archivo
df <- read.csv("SraRunTable_metagenome.csv")


# Cargar los excels con spots (ejemplo, si tienes uno por BioProject)
df_spots1 <- read_csv("PRJNA1015591_runinfo.csv")
df_spots2 <- read_csv("PRJNA747886_runinfo.csv")
df_spots3 <- read_csv("PRJNA755178_runinfo.csv")
df_spots4 <- read_csv("PRJNA786222_runinfo.csv")
df_spots5 <- read_csv("PRJNA898075_runinfo.csv")

# Unir todos los excels en una sola tabla
df_spots <- bind_rows(df_spots1, df_spots2, df_spots3, df_spots4, df_spots5)

# Asegurarnos de que la columna Run esté como carácter
df$Run <- as.character(df$Run)
df_spots$Run <- as.character(df_spots$Run)

# Hacer el merge por Run
df_merged <- df %>%
  left_join(df_spots, by = "Run")

# Ahora ya tienes spots en tu tabla principal
head(df_merged)


# Filtrar según condiciones
df_filtrado <- df_merged %>%
  filter(
    GENE == "16S",
    BioProject.x %in% c("PRJNA1015591", "PRJNA747886", "PRJNA755178", "PRJNA786222", "PRJNA898075"),
    toupper(LibraryLayout.x) %in% c("PAIRED")
  )

# Exportar listado de Run a un archivo .txt
write.table(
  df_filtrado$Run,
  file = "runs_list.txt",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)


write.csv(df_filtrado, "df_filtrado.csv", row.names = FALSE)


# Número de muestras por BioProject
df_filtrado %>%
  group_by(BioProject.x) %>%
  summarise(num_muestras = n_distinct(Run))


# Crear un resumen por BioProject
df_summary <- df_filtrado %>%
  group_by(BioProject.x) %>%
  summarise(
    total_bases = sum(Bases, na.rm = TRUE),
    num_lecturas = sum(spots, na.rm = TRUE),
    longitud_media = mean(AvgSpotLen, na.rm = TRUE),
    num_muestras = n_distinct(Run)
  )

# Gráfico 1: Total de bases por BioProject
ggplot(df_summary, aes(x = BioProject.x, y = total_bases)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Total de bases secuenciadas por BioProject",
       x = "BioProject", y = "Total de bases") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Gráfico 2: Número de lecturas
ggplot(df_summary, aes(x = BioProject.x, y = num_lecturas)) +
  geom_bar(stat = "identity", fill = "seagreen") +
  labs(title = "Número total de lecturas por BioProject",
       x = "BioProject", y = "Número de lecturas") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Gráfico 3: Longitud media de lectura
ggplot(df_summary, aes(x = BioProject.x, y = longitud_media)) +
  geom_bar(stat = "identity", fill = "darkorange") +
  labs(title = "Longitud media de lectura por BioProject",
       x = "BioProject", y = "Longitud media (bp)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Gráfico 4: Número de muestras
ggplot(df_summary, aes(x = BioProject.x, y = num_muestras)) +
  geom_bar(stat = "identity", fill = "purple") +
  labs(title = "Número de muestras por BioProject",
       x = "BioProject", y = "Número de muestras") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Gráfico 5: Distribución del número de muestras por BioProject
ggplot(df_filtrado, aes(x = BioProject.x, y = spots, fill = BioProject.x)) +
  geom_boxplot(fill = "grey") +
  theme_minimal() +
  labs(
    title = "Distribución del número de lecturas por muestra",
    x = "BioProject",
    y = "Número de lecturas"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



## # Downsampling a ~50,000 lecturas por muestra en cada BioProject calculando la desviación estándar

# Calcular media y sd original
estadisticas <- df_filtrado %>%
  group_by(BioProject.x) %>%
  summarise(
    mean_spots = mean(spots, na.rm = TRUE),
    sd_spots   = sd(spots, na.rm = TRUE),
    .groups = "drop"
  )
estadisticas

# Downsampling con preservación de la sd original
df_downsampled_1 <- df_filtrado %>%
  left_join(estadisticas, by = "BioProject.x") %>%
  mutate(
    spots_downsampled = round(((spots - mean_spots) / sd_spots) * sd_spots + 50000)
  ) %>%
  select(-mean_spots, -sd_spots)

#corroborar 
df_downsampled_1 %>%
  group_by(BioProject.x) %>%
  summarise(
    mean_new = mean(spots_downsampled, na.rm = TRUE),
    sd_new   = sd(spots_downsampled, na.rm = TRUE),
    .groups = "drop"
  )

# Gráfico
ggplot(df_downsampled_1, aes(x = BioProject.x, y = spots_downsampled, fill = BioProject.x)) +
  geom_boxplot(fill = "grey") +
  theme_minimal() +
  labs(
    title = "Distribución del número de lecturas por muestra post-downsampling",
    x = "BioProject",
    y = "Número de lecturas"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

df_downsampled_1$spots_downsampled

##Imputar datos negativos

# 1. Calcular la media global solo con valores positivos
media_valida <- mean(df_downsampled_1$spots_downsampled[df_downsampled_1$spots_downsampled > 0], na.rm = TRUE)

# 2. Reemplazar los valores negativos con esa media
df_downsampled_1 <- df_downsampled_1 %>%
  mutate(
    spots_downsampled = ifelse(spots_downsampled < 0, round(media_valida), spots_downsampled)
  )

df_downsampled_1$spots_downsampled

write.csv(df_downsampled_1, "df_downsampled_1.csv", row.names = FALSE)

# Gráfico
ggplot(df_downsampled_1, aes(x = BioProject.x, y = spots_downsampled, fill = BioProject.x)) +
  geom_boxplot(fill = "grey") +
  theme_minimal() +
  labs(
    title = "Distribución del número de lecturas por muestra post-downsampling",
    x = "BioProject",
    y = "Número de lecturas"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Seleccionar variables de interés para qiime

df_qiime <- df_downsampled_1 %>%
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

# Modificar el HOST
df_qiime <- df_qiime %>% 
  mutate(HOST = ifelse(HOST == "" | HOST == "not applicable", "Theobroma cacao", HOST)) %>% 
  mutate(`Env-broad` == "")

# Ajustar el env_broad
df_qiime <- df_qiime %>%
  mutate(
    `Env-broad` = case_when(
      # Reemplazar por cacao_soil
      `Env-broad` %in% c("ENVO_01000763", "Mineral", "cacao plantation") ~ "Cacao soil",
      
      # Si `Loc-name` es Indonesia o Bolivia Y `Env-broad` está vacío
      `Loc-name` %in% c("Indonesia", "Bolivia") & `Env-broad` %in% c("", NA) ~ "Cacao soil",
      
      # Reemplazar tropical lowland...
      `Env-broad` == "tropical lowland evergreen broadleaf rain forest" ~ "Other soil",
      
      # Reemplazar por greenhouse si es USA
      `Loc-name` == "USA" ~ "Greenhouse",
      
      # Mantener los demás valores sin cambios
      TRUE ~ `Env-broad`
    )
  )


# Eliminar columna llamada ``Env-broad` == ""`
df_qiime <- df_qiime %>% select(-`\`Env-broad\` == ""`)


set.seed(123) # reproducibilidad

# Función para generar secuencias aleatorias de ADN
random_dna <- function(n, length = 12) {
  bases <- c("A", "T", "C", "G")
  replicate(n, paste0(sample(bases, length, replace = TRUE), collapse = ""))
}

# Agregar columna Barcodes
df_qiime$Barcodes <- random_dna(nrow(df_qiime), length = 12)

# Reordenar columnas para que Barcodes quede después de AvgSpotLen
df_qiime <- df_qiime %>%
  relocate(Barcodes, .after = AvgSpotLen)


#Exportar

# Exportar como archivo .tsv (separado por tabuladores, sin comillas)
write.table(df_qiime,
            file = "metadata_metagenomePE16S.txt",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)
