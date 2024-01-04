library(readxl)
library(tidyverse)
library(ggplot2)
library(vegan)
library(factoextra)

# Delete all variables.
rm(list = ls())

# Set working directory
setwd("/home/julia/Documentos/PhD/Data_analysis/PNAE")
#setwd("/home/fbartu/Research/Julia_Rodriguez/PNAE")

# Read data
birds_UTM1x1 <- read.csv("Data/Birds_data/birds_UTM1x1.csv")
pipiens_UTM1x1 <- read.csv("Data/Mosquito_data/Culex_pipiens/pipiens_UTM1x1.csv")
pv <- read.csv("Data/Risk_data/predicted_values.csv")

# Scatter plot with regression line
ggplot(pv, 
       aes(x = predicted_value_birds, y = predicted_value_pipiens)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Gráfico de correlación", x = "Birds predicted value", y = "Culex pipiens predicted value")

rm(birds_UTM1x1, pipiens_UTM1x1, pv)

################################################################################
# LAND USES
################################################################################

# Land uses UTM1x1
usosQUTM1k <- read_delim("Data/Land_uses/usos2017_qutm1k.csv", delim =";", escape_double = FALSE, trim_ws = TRUE)
colnames(usosQUTM1k)[colnames(usosQUTM1k) == "QUTM1K"] <- "UTM1X1"

usosQUTM1k$urba <- usosQUTM1k$`4` + usosQUTM1k$`5` + usosQUTM1k$`6` + usosQUTM1k$`7`
usosQUTM1k$herbacis <- usosQUTM1k$`8` + usosQUTM1k$`9` 
usosQUTM1k$wetland <- usosQUTM1k$`1` + usosQUTM1k$`20` 
usosQUTM1k$ricefield <- usosQUTM1k$`24`
usosQUTM1k$total <- rowSums(usosQUTM1k[,2:26])
usosQUTM1k$others <- usosQUTM1k$total - usosQUTM1k$urba - usosQUTM1k$herbacis - usosQUTM1k$wetland - usosQUTM1k$ricefield

usosQUTM1k <- usosQUTM1k %>% dplyr::select(UTM1X1,urba,herbacis,wetland,ricefield,others,total)

# Comparison of land uses in high vs low risk cells
high_risk_cells <- read_csv("Data/Risk_data/high_risk_cells.csv")
high_risk_cells <- merge(high_risk_cells, usosQUTM1k, by="UTM1X1", all=F)
high_risk_cells <- high_risk_cells %>% 
  dplyr::select(urba, herbacis, wetland, ricefield)
high_risk_cells <- high_risk_cells %>% mutate(risk = 'H')
#high_risk_cells <- as.data.frame(colMeans(high_risk_cells))

low_risk_cells <- read_csv("Data/Risk_data/low_risk_cells.csv")
low_risk_cells <- merge(low_risk_cells, usosQUTM1k, by="UTM1X1", all=F)
low_risk_cells <- low_risk_cells %>% 
  dplyr::select(urba, herbacis, wetland, ricefield)
low_risk_cells <- low_risk_cells %>% mutate(risk = 'L')
#low_risk_cells <- as.data.frame(colMeans(low_risk_cells))

land_uses <- rbind(high_risk_cells, low_risk_cells) #all cells

# Pivot longer: multivariate analysis with two categorical variables
land_uses <- land_uses %>%
  pivot_longer(c(urba, herbacis, wetland, ricefield), 
               names_to = "land_uses", values_to = "area") 

# Normality
# Histograms per land use
land_uses %>% filter(land_uses == "urba") %>% ggplot() + 
  geom_histogram(aes(area))
land_uses %>% filter(land_uses == "herbacis") %>% ggplot() + 
  geom_histogram(aes(area))
land_uses %>% filter(land_uses == "wetland") %>% ggplot() + 
  geom_histogram(aes(area))
land_uses %>% filter(land_uses == "ricefield") %>% ggplot() + 
  geom_histogram(aes(area))

# Box plot to compare high and low risk cells: no normality
ggplot(land_uses) + geom_boxplot(aes(y=area, x=land_uses, color=risk))
#filename <- "Plots/boxplot_land_uses.jpg"
#ggsave(filename, dpi = 300, width = 10, height = 6, units = "in", type = "jpg", quality = 100)

# Pivot wider
land_uses <- rbind(high_risk_cells, low_risk_cells)

# Normality test
shapiro.test(land_uses$urba)
shapiro.test(land_uses$herbacis)
shapiro.test(land_uses$wetland)
shapiro.test(land_uses$ricefield)

# NP multivariate ANOVA 
Y <- land_uses[, c("urba", "herbacis", "wetland", "ricefield")]
result <- adonis2(Y ~ land_uses$risk, 
                  method = "euclidean", 
                  permutations = 999)

# Mann-Whitney (per 2 grups)
wilcox.test(urba ~ risk, data = land_uses)
wilcox.test(herbacis ~ risk, data = land_uses)
wilcox.test(wetland ~ risk, data = land_uses)
wilcox.test(ricefield ~ risk, data = land_uses)

rm(high_risk_cells, low_risk_cells, result, Y)

################################################################################
# LAND USES KMEANS
################################################################################
high_risk_cells <- read_csv("Data/Risk_data/high_risk_cells.csv")
usosQUTM1k_hrc <- usosQUTM1k %>% 
  filter(UTM1X1 %in% c(high_risk_cells$UTM1X1)) %>%
  select(UTM1X1, urba, herbacis, wetland, ricefield) 
usos <- as.matrix(usosQUTM1k_hrc[,2:5])

fviz_nbclust(usos, kmeans, method = "wss") #elbow method
set.seed(123)
k <- 3  # Number of clusters
kmeans_utm <- kmeans(usos, centers = k)
fviz_cluster(kmeans_utm, usos, ellipse.type = "norm") #plot clusters

usosQUTM1k_hrc$cluster <- kmeans_utm$cluster
#write_csv(usosQUTM1k_hrc, "Data/Risk_data/clusters_hrc.csv")

rm(kmeans_utm, usos, k)
