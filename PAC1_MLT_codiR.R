# Paquets necessaris
library(knitr)
library(ggplot2)
library(reshape2)
library(Biobase)
library(SummarizedExperiment)
library(dplyr)
library(tidyr)
library(stats)
library(tibble)
library(pheatmap)
library(FactoMineR)
library(factoextra)
library(ggExtra)
library(mixOmics)

# 1. Selecció d'un dataset de metabolòmica.
url <- "https://rest.xialab.ca/api/download/metaboanalyst/human_cachexia.csv"
cachexia_data <- read.csv(url, header = TRUE)
kable(head(cachexia_data, 10), format = "markdown")

# 2. Creació de SummarizedExperiment amb les dades i metadades.
expressionValues <- t(as.matrix(cachexia_data[, 3:65]))  # Dades metabolòmiques
colData <- data.frame(
  Patient.ID = cachexia_data$Patient.ID,
  Muscle.loss = cachexia_data$Muscle.loss,
  row.names = paste0("S", 1:77)
)
rowData <- data.frame(Metabolite = rownames(expressionValues))

# Assignar noms de files i columnes
colnames(expressionValues) <- paste0("S", 1:77)
rownames(expressionValues) <- paste0("M", 1:nrow(expressionValues))

# Crear objecte SummarizedExperiment
se_cachexia <- SummarizedExperiment(
  assays = list(raw = expressionValues),
  colData = colData,
  rowData = rowData
)

# Afegir metadades a l'objecte SummarizedExperiment
metadata(se_cachexia) <- list(study_description = "A metabolite concentration table from human urine samples with two groups")

# Mostrar objecte SummarizedExperiment creat
se_cachexia

# 3. Anàlisi exploratòria del dataset
# Dimensions de la matriu d'expressió
dim(assays(se_cachexia)$raw)

# Metadades de les files (metabolits)
head(rowData(se_cachexia))

# Metadades de les columnes (pacients)
head(colData(se_cachexia))

# Selecció de dades de les primeres 10 files i 5 columnes
assays(se_cachexia)[[1]][1:10, 1:5]

# 4. Exploració univariant: Boxplot de les concentracions dels metabolits
long_data <- melt(assays(se_cachexia)$raw)
colnames(long_data) <- c("Metabolite", "Sample", "Concentration")

ggplot(long_data, aes(x = Metabolite, y = Concentration)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 16) + 
  labs(title = "Distribució de les concentracions dels metabòlits", x = "Metabolits", y = "Concentració") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.title = element_text(hjust = 0.5))

# 5. Estadístiques univariants dels metabolits
data_long <- assays(se_cachexia)$raw %>%
  as.data.frame() %>%
  rownames_to_column("Metabolite") %>%
  pivot_longer(cols = -Metabolite, names_to = "Sample", values_to = "Concentration") %>%
  mutate(Group = colData(se_cachexia)$Muscle.loss[match(Sample, rownames(colData(se_cachexia)))])

pos_outcome <- "cachexic"

univar <- function(data, group_col, value_col, posclass) {
  data %>%
    group_by(Metabolite) %>%
    summarise(
      TotalMissing = sum(is.na(.data[[value_col]])),
      TTestPvalue = t.test(.data[[value_col]] ~ .data[[group_col]])$p.value,
      Cachexic_Mean = mean(.data[[value_col]][.data[[group_col]] == posclass], na.rm = TRUE),
      Control_Mean = mean(.data[[value_col]][.data[[group_col]] != posclass], na.rm = TRUE)
    )
}

stats_table <- univar(data_long, group_col = "Group", value_col = "Concentration", posclass = pos_outcome)

# Mostrar la taula d'estadístiques
kable(stats_table, caption = "Estadístiques dels Metabolits", digits = c(0, 0, 3, 2, 2), format = "markdown")

# Afegir la diferència entre les mitjanes de Cachexia i Control
stats_table <- stats_table %>%
  mutate(Difference = Cachexic_Mean - Control_Mean)

# Mostrar els 10 metabòlits amb la major diferència
stats_table_filtered <- stats_table %>%
  dplyr::select(Metabolite, Cachexic_Mean, Control_Mean, Difference)

top_10_diff <- stats_table_filtered %>%
  arrange(desc(abs(Difference))) %>%
  head(10)

# Mostrar per pantalla els 10 metabòlits amb la major diferència
print(top_10_diff)

# 6. Selecció de metabolits significatius (FDR < 0.05)
stats_table <- stats_table %>%
  mutate(FDR = p.adjust(TTestPvalue, method = "fdr"))

significant_metabolites <- stats_table %>%
  filter(FDR < 0.05) %>%
  arrange(FDR) %>%
  pull(Metabolite)

significant_metabolites
length(significant_metabolites)

# 7. Boxplot per als metabolits significatius
filtered_data <- assays(se_cachexia)$raw[rownames(assays(se_cachexia)$raw) %in% significant_metabolites, , drop = FALSE]
long_data_filtered <- melt(filtered_data)
colnames(long_data_filtered) <- c("Metabolite", "Sample", "Concentration")

long_data_filtered$Group <- rep(colData(se_cachexia)$Muscle.loss, times = nrow(filtered_data))

ggplot(long_data_filtered, aes(x = Metabolite, y = Concentration, fill = Group)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 16) +
  labs(title = "Comparació de concentració de metabòlits entre grups", x = "Metabolits", y = "Concentració") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  facet_wrap(~ Group, scales = "free_y", ncol = 1) + 
  ylim(0, 35000)

# 8. PCA
pca_result <- PCA(t(expressionValues), scale.unit = TRUE, graph = FALSE)

pca_plot <- fviz_pca_ind(
  pca_result,
  geom.ind = "point",
  col.ind = colData(se_cachexia)$Muscle.loss,
  addEllipses = TRUE,
  title = "Anàlisi PCA: Cachexia vs Control",
  legend.title = "Muscle Loss"
)

pca_ggplot <- pca_plot + theme_minimal()

# Afegir els eixos marginals
ggMarginal(
  pca_ggplot,
  type = "density",
  groupColour = TRUE,
  groupFill = TRUE
)

# 9. PCA amb 3 components (mixOmics)
df_cachexia <- t(assays(se_cachexia)$raw)
tune.pca.result <- tune.pca(df_cachexia, ncomp = 10, scale = TRUE)

plot(tune.pca.result, col = "lightblue")

# Realitzar el PCA amb 3 components
pca.final <- pca(df_cachexia, ncomp = 3, center = TRUE, scale = TRUE)

# Visualitzar les mostres projectades
plotIndiv(pca.final, comp = c(1, 2), legend = FALSE, ind.names = TRUE, 
          group = se_cachexia$Muscle.loss, title = "PCA mostres (PC1 i PC2)")

# Gràfic interactiu 3D
plotIndiv(pca.final, style = '3d', group = se_cachexia$Muscle.loss, title = 'PCA comp 1 - 3')

# 10. Gràfic de variables de PCA
plotVar(pca.final, comp = c(1, 2), var.names = TRUE, cex = 3, title = "PCA metabolits (PC1 i PC2)")

# 11. Mapa de calor de correlació entre metabolits
cor_matrix <- cor(t(assay(se_cachexia)))
pheatmap(cor_matrix, 
         main = "Mapa de calor de la correlació entre metabolits", 
         cluster_rows = TRUE, cluster_cols = TRUE, 
         show_rownames = FALSE, show_colnames = FALSE)
