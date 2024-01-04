library(readxl)
library(tidyverse)
library(ggplot2)
library(lme4)
library(lubridate)
library(zoo)
library(performance)
library(XNomial)
library(ggpubr)
library(rstatix)
library(janitor)
library(vegan)
library(factoextra)

# Delete all variables.
rm(list = ls())

# Set working directory
setwd("/home/julia/Documentos/PhD/Data_analysis/PNAE")
#setwd("/home/fbartu/Research/Julia_Rodriguez/PNAE")

# Read birds data and filter by study site
birds_data <- read.csv("Data/Birds_data/Dades_ICO_models_i_zonacions_WNV.csv", fileEncoding="latin1", sep=";")
colnames(birds_data)[colnames(birds_data) == "cell"] <- "UTM1X1"
birds_data <- birds_data %>% 
  filter(tipus.espècie == "Espècies positives per WNV i presents als Aiguamolls")

area_interes <- read_xlsx("Data/Birds_data/area_interes.xlsx")

birds_data_aiguamolls <- merge(birds_data,area_interes,by="UTM1X1",all = F)
birds_data_aiguamolls <- birds_data_aiguamolls %>%
  dplyr::select(UTM1X1, species_id, latin_name_ATLAS, predicted_value, 
         nom.científic, nom.castellà)

# Predicted value to 0-1 value
birds_data_aiguamolls$predicted_value <- birds_data_aiguamolls$predicted_value/1000

# Abundance by UTM1x1
birds_UTM1x1 <- birds_data_aiguamolls %>% 
  group_by(UTM1X1) %>% 
  summarise(predicted_value = sum(predicted_value))

# Normalized predicted value: (x-xmin)/(xmax-xmin)
birds_UTM1x1$predicted_value <- 
  (birds_UTM1x1$predicted_value - min(birds_UTM1x1$predicted_value)) / 
  (max(birds_UTM1x1$predicted_value) - min(birds_UTM1x1$predicted_value))

colnames(birds_UTM1x1)[colnames(birds_UTM1x1) == "predicted_value"] <- "predicted_value_birds"

write.csv(birds_UTM1x1, "Data/Birds_data/birds_UTM1x1.csv")

# Risk index calculation -------------------------------------------------------
pipiens_UTM1x1 <- read.csv("Prediction/Culex/culex_UTM1x1.csv")

pv <- merge(birds_UTM1x1, pipiens_UTM1x1, by="UTM1X1")
pv <- pv %>% dplyr::select("UTM1X1", "predicted_value_birds", "predicted_value_culex")

# Calculate risk index by summing predicted values
pv$risk_index <- pv$predicted_value_birds + pv$predicted_value_culex

pv$risk_index <- (pv$risk_index - min(pv$risk_index)) / 
  (max(pv$risk_index) - min(pv$risk_index)) #normalized value

write.csv(pv, "Data/Risk_data/predicted_values.csv")

################################################################################
# TAULA RDA 
################################################################################
# Pivot wider
rda_taula <- birds_data_aiguamolls %>%
  dplyr::select(UTM1X1, latin_name_ATLAS, predicted_value) %>%
  transform(predicted_value = as.numeric(predicted_value)) %>%
  pivot_wider(names_from = latin_name_ATLAS, values_from = predicted_value, values_fill = 0)

# Land uses UTM1x1
usosQUTM1k <- read_delim("Data/Land_uses/usos2017_qutm1k.csv", delim =";", escape_double = FALSE, trim_ws = TRUE)
colnames(usosQUTM1k)[colnames(usosQUTM1k) == "QUTM1K"] <- "UTM1X1"

usosQUTM1k$urba <- usosQUTM1k$`4` + usosQUTM1k$`5` + usosQUTM1k$`6` + usosQUTM1k$`7`
usosQUTM1k$herbacis <- usosQUTM1k$`8` + usosQUTM1k$`9` 
usosQUTM1k$wetland <- usosQUTM1k$`1` + usosQUTM1k$`20` 
usosQUTM1k$ricefield <- usosQUTM1k$`24`
usosQUTM1k$total <- rowSums(usosQUTM1k[,2:26])
usosQUTM1k$others <- usosQUTM1k$total - usosQUTM1k$urba - usosQUTM1k$herbacis - usosQUTM1k$wetland - usosQUTM1k$ricefield

usosQUTM1k <- usosQUTM1k %>% dplyr::select(UTM1X1,urba,herbacis,wetland,ricefield)

# Filter by area of interest
usosQUTM1k <- merge(usosQUTM1k,area_interes,by="UTM1X1",all = F) %>%
  dplyr::select(UTM1X1,urba,herbacis,wetland,ricefield)

# Merge both data frames
rda_taula <- merge(rda_taula, usosQUTM1k, by="UTM1X1")

# Add risk value
rda_taula <- merge(rda_taula, pv, by="UTM1X1")
rda_taula <- rda_taula %>% 
  dplyr::select(-"predicted_value_birds", -"predicted_value_culex") %>%
  clean_names()

#writexl::write_xlsx(rda_taula, "Data/Birds_data/rda_taula.xlsx")

# Delete data frames
rm(area_interes, birds_UTM1x1, birds_data, pv)

# Most abundant species in the high risk cells
high_risk_cells <- read_csv("Data/Risk_data/high_risk_cells.csv")
species_hrc <- rda_taula %>% dplyr::filter(utm1x1 %in% c(high_risk_cells$UTM1X1))
species_hrc <- species_hrc %>% dplyr::select(1:56)
species_hrc_ranking <- as.data.frame(colMeans(species_hrc[, -1]))

################################################################################
# RDA MODEL
################################################################################
species <- as.matrix(rda_taula[,2:56]) #matrix
rda_model <- rda(species ~ urba + herbacis + wetland + ricefield, 
                 data=rda_taula, 
                 scale = TRUE, 
                 na.rm = TRUE) #RDA model, land uses as explanatory variables

rda_plot <- plot(rda_model, type='n', scaling=1)
orditorp(rda_model, display='sp', cex=0.5, scaling=1, col='blue')
text(rda_model, display='cn', col='red') 

rda_plot <- plot(rda_model)

species_rda_axis <- as.data.frame(rda_plot$species) #dataframe of the species and its values on the RDA axes

rm(rda_model,rda_plot,species,species_hrc_ranking)

################################################################################
# KMEANS
################################################################################
# Using RDA axes
fviz_nbclust(species_rda_axis, kmeans, method = "wss") #elbow method
set.seed(123)
k <- 3  # Number of clusters
kmeans_rda_species <- kmeans(species_rda_axis, centers = k)
species_clustering_plot <- fviz_cluster(kmeans_rda_species, species_rda_axis, 
                     ellipse.type = "norm",
                     stand = TRUE,
                     axes = c(1, 2),
                     repel = TRUE,
                     pointsize = 1,
                     labelsize = 15,
                     main = "Bird species clustering based on RDA axes",
                     legend.title = "Clusters",
                     font.main = c(15,"bold","black"),
                     font.x = c(14,"plain","black"),
                     font.y = c(14,"plain","black"),
                     font.legend = c(14,"plain","black"))

#filename <- "Plots/species_clustering_plot.jpg"
#ggsave(filename, dpi = 300, width = 15, height = 8, units = "in", type = "jpg", quality = 100)

species_rda_axis$cluster <- kmeans_rda_species$cluster

for (i in 1:length(species_rda_axis$cluster)) {
  if (species_rda_axis$cluster[i] == 1) species_rda_axis$cluster[i]="urba"
  if (species_rda_axis$cluster[i] == 2) species_rda_axis$cluster[i]="wetland"
  if (species_rda_axis$cluster[i] == 3) species_rda_axis$cluster[i]="herbacis"}

rm(kmeans_rda_species)

################################################################################
# HABITAT SELECTION INDEX
################################################################################
habitat_selection <- read_excel("Data/Birds_data/habitat_selection.xlsx")

# Generalists vs specialists
habitat_selection <- habitat_selection %>% 
  dplyr::select(Especie, Zones_hum_i, urba_suburba_i, Heb_reg_i, Herb_seca_i)

# Bird species that select positively all habitats
generalistes <- habitat_selection[
  apply(habitat_selection[, 2:5] > 0, 1, all), ] 

# Bird species that select positively wet and urban areas
generalistes_wet_urba <- habitat_selection[
  (habitat_selection[, 2] > 0 & habitat_selection[, 3] > 0) | 
    (habitat_selection[, 3] > 0 & habitat_selection[, 4] > 0), ]

################################################################################
# WOAH
################################################################################

woah <- read_excel("WOAH.xlsx")

woah <- woah %>% 
  dplyr::select("Año", "País", "División administrativa", "Categoría animal", 
                "Especie", "Event_id", "Outbreak_id")

woah <- woah %>% filter(!(Especie %in% c("-", "Fauna silvestre (especie no especificada)", "Aves", 
                                         "Jabalí", "Perros", "Équidos", "Camélidos", "Corzo",
                                         "Accipitridae (no identificada)", "Columbidae (no identificada)",
                                         "Ardeidae (no identificada)", "Turdidae (no identificada)", 
                                         "Passeridae (no identificada)", "Charadriidae (no identificada)",
                                         "Psittacidae (no identificada)", "Phasianidae (no identificada)",
                                         "Podicipedidae (no identificada)", "Tytonidae (no identificada)",
                                         "Strigidae (no identificada)", "Falconidae (no identificada)",
                                         "Laridae (no identificada)", "Spheniscidae (no identificada)",
                                         "Phoenicopteridae (no identificada)", "Pelecanidae (no identificada)",
                                         "Bombycillidae (no identificada)", "Anserinae (no identificada)",
                                         "Corvidae (no identificada)")))

