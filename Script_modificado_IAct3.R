
######################################################
#      ACT 3: APREDIZAJE AUTOMÁTICO                  #
#    GRUPO 22: Laura Córdova, Xiomira Fiallos,       #
#  Lizbeth Navarrete, Soledad Ortega,                #
#  Mirkka Puente, Melanie Polo.                      #
######################################################

##--------------------------------------------------------------------
#             Limpieza, imputación, exploración de datos
##--------------------------------------------------------------------

rm(list=ls())
expression <- "./gene_expression.csv"
classes <- "./classes.csv"
names <- "./column_names.txt"

# Leer el documento de gene_expression 
gene_expression <- read.csv(expression, sep = ";", header = FALSE)

# Leer el documento del nombre de los genes
column_names <- readLines("column_names.txt")

# Leer las clasificacion y el sample number
classes <- read.csv(classes, sep = ";", header = FALSE, stringsAsFactors = TRUE)

# Asignar el nombre de las columnas 
colnames(gene_expression) <- column_names
colnames(classes) <- c("sample_ID", "class")

# Add los sample numbers and class
gene_expression$sample_ID <- classes$sample_ID  
gene_expression$class <- classes$class  

# Lista de paquetes necesarios
packs <- c(
  "glmnet",      # Elastic Net
  "tidyverse",   # Data wrangling, ggplot2, etc.
  "caret",       # Machine Learning
  "rpart",       # Decision Trees
  "rpart.plot",  # Plot de árboles
  "rattle",      # Visualización de árboles
  "pROC",        # ROC
  "PRROC",       # Precision-Recall Curve
  "MASS",        # LDA
  "klaR",        # RDA
  "gridExtra",   # Combinar gráficos
  "dplyr",       # Manipulación de datos
  "Rtsne",       # Reduccion de dimensionalidad
  "randomForest", # Machine Learnig 
  "e1071",       # SVM
  "factoextra"   # Clustering viz
)

# Comprobación automática y carga
lapply(packs, library, character.only = TRUE)

# Existe NA en nuestros datos?
cat("Total de NAs iniciales:", sum(is.na(gene_expression)), "\n")


##--------------------------------------------------------------------
#             Revision de normalidad de gene expression
##--------------------------------------------------------------------

prueba_shapiro<- lapply(column_names, function(g) {
  x <- as.numeric(gene_expression[[g]])
  if (length(na.omit(x)) < 3 || sd(x, na.rm = TRUE) == 0) {
    return(data.frame(gene = g, W = NA, p_value = NA))
  }
  test <- shapiro.test(x)
  data.frame(gene = g, W = test$statistic, p_value = test$p.value)
})

shapiro_table <- bind_rows(prueba_shapiro) %>%
  mutate(normalidad = ifelse(p_value >= 0.05, "Normalidad", "No normalidad")) %>%
  arrange(p_value)

# Data set sin los genes que no cumplen requisitos
columnas_a_excluir <- c("MIER3", "ZCCHC12", "RPL22L1", "sample_ID", "class")
df_genes <- gene_expression[, !names(gene_expression) %in% columnas_a_excluir]

##--------------------------------------------------------------------
#                      Escala de Datos
##--------------------------------------------------------------------

# Escalar Datos (Log2 + Scale)
datos_escalados <- log2(df_genes + 1) %>% scale()
datos_escalados <- as.data.frame(datos_escalados)

# Add los sample numbers and class a los datos escalados
datos_escalados$sample_ID <- classes$sample_ID  
datos_escalados$class <- classes$class

##-----------------------------------------------------------------------
#    1. Métodos no supervisados - Reduccion de dimensionalidad 
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
#                        Método 1: PCA
##-----------------------------------------------------------------------

columnas_no_numericas <- c("sample_ID", "class")
data_pca <- datos_escalados[, !names(datos_escalados) %in% columnas_no_numericas]

pca.results <- prcomp(data_pca, center = TRUE, scale. = FALSE) 

pca.df <- data.frame(pca.results$x)
pca.df$class <- datos_escalados$class 

varianzas <- pca.results$sdev^2
total.varianza <- sum(varianzas)
varianza.explicada <- varianzas / total.varianza

x_label <- paste0('PC1 (', round(varianza.explicada[1] * 100, 2), '%)')
y_label <- paste0('PC2 (', round(varianza.explicada[2] * 100, 2), '%)')

ggplot(pca.df, aes(x = PC1, y = PC2, color = class)) +
  geom_point(size = 3, alpha = 0.8) + 
  labs(title = 'PCA - Expresion Génica', x = x_label, y = y_label, color = 'Clase') +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "gray90"), plot.title = element_text(hjust = 0.5))

##--------------------------------------------------------------------
#    Método 2: t-Distributed Stochastic Neighbor Embedding (t-SNE)
##--------------------------------------------------------------------

datos_tsne_clean <- datos_escalados %>% distinct()
matrix_tsne <- as.matrix(datos_tsne_clean[, !names(datos_tsne_clean) %in% columnas_no_numericas])
labels_tsne <- datos_tsne_clean$class

set.seed(1995) 

tsne_out <- Rtsne(X = matrix_tsne, dims = 2, check_duplicates = FALSE)

tsne_result <- data.frame(tsne_out$Y)
colnames(tsne_result) <- c("Dim1", "Dim2")
tsne_result$class <- labels_tsne 

ggplot(tsne_result, aes(x = Dim1, y = Dim2, color = class)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = paste0("t-SNE - Datos expresion génica "), x = "Dimensión 1", y = "Dimensión 2", color = "Clase") +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "gray90"), plot.title = element_text(hjust = 0.5))

##--------------------------------------------------------------------
#         2.  Métodos no supervisados - Clusterización
##--------------------------------------------------------------------
##--------------------------------------------------------------------
#               Método 1: Clustering no jerárquico (K-means)
##--------------------------------------------------------------------

set.seed(1995) 
kmeans_result <- kmeans(data_pca, centers = 5, iter.max = 100, nstart = 25)

fviz_cluster(kmeans_result, data = data_pca, 
             geom = "point", ellipse.type = "convex",
             main = "Gráfico K-means (k=5)", ggtheme = theme_minimal())

table(Clúster = kmeans_result$cluster, Clase_Real = datos_escalados$class)

##--------------------------------------------------------------------
#               Método 2: Clustering Jerárquico 
##--------------------------------------------------------------------

dist_matrix <- dist(data_pca, method = "euclidean")
hclust_ward <- hclust(dist_matrix, method = "ward.D2")

fviz_dend(hclust_ward, k = 5, cex = 0.5, k_colors = "jco", rect = TRUE, 
          main = "Dendrograma (Método Ward)", xlab = "Índice", ylab = "Distancia") + theme_classic()

grupos_hclust <- cutree(hclust_ward, k = 5)
table(Clúster_Jerárquico = grupos_hclust, Clase_Real = datos_escalados$class)


##--------------------------------------------------------------------
#                      3. Métodos supervisados  
##--------------------------------------------------------------------

# ====================================================================
# DIVISIÓN DE DATOS (80/20)
# ====================================================================
# Quitamos sample_ID para modelado
datos_modelo <- datos_escalados %>%
  dplyr::select(-all_of("sample_ID")) %>%
  dplyr::mutate(class = factor(class))

set.seed(1995) # Semilla ÚNICA para todo el grupo
idx_global <- caret::createDataPartition(datos_modelo$class, p = 0.80, list = FALSE)

# Estos son los datos OFICIALES para todos los modelos
train_global <- datos_modelo[idx_global, , drop = FALSE]
test_global  <- datos_modelo[-idx_global, , drop = FALSE]

cat("--- Dimensiones UNIFICADAS (80/20) ---\n")
cat("Train:", dim(train_global), "\n")
cat("Test: ", dim(test_global), "\n")


##--------------------------------------------------------------------
#                        Método 1: Random Forest
##--------------------------------------------------------------------

set.seed(1995)

# 1) Asignación de datos globales 
#    split global para ser justos con los modelos
train_df <- train_global
test_df  <- test_global

# 2) Eliminar predictores de varianza ~0 
nzv_idx <- caret::nearZeroVar(train_df %>% dplyr::select(-class))
if (length(nzv_idx) > 0) {
  genes_keep <- setdiff(colnames(train_df), c("class"))[ -nzv_idx ]
  train_df <- dplyr::select(train_df, dplyr::all_of(c("class", genes_keep)))
  test_df  <- dplyr::select(test_df,  dplyr::all_of(c("class", genes_keep)))
}

# 3) Imputación mediana 
pp <- caret::preProcess(train_df %>% dplyr::select(-class), method = c("medianImpute"))
train_x <- predict(pp, train_df %>% dplyr::select(-class))
test_x  <- predict(pp, test_df  %>% dplyr::select(-class))
train_y <- train_df$class
test_y  <- test_df$class

# 4) Entrenamiento Random Forest
p <- ncol(train_x)
mtry_grid <- unique(pmax(1, round(c(sqrt(p)/2, sqrt(p), sqrt(p)*2))))

ctrl <- caret::trainControl(
  method = "repeatedcv", number = 5, repeats = 2,
  classProbs = FALSE, verboseIter = TRUE, allowParallel = TRUE
)

set.seed(1995)
rf_fit <- caret::train(
  x = train_x, y = train_y,
  method = "rf",
  trControl = ctrl,
  tuneGrid = data.frame(mtry = mtry_grid),
  ntree = 500, # Ajustado a 500 por eficiencia estándar
  importance = TRUE,
  metric = "Accuracy"
)

# 5) Evaluación RF
pred_rf <- predict(rf_fit, newdata = test_x)
cm_rf <- caret::confusionMatrix(pred_rf, test_y)

# Extracción de métricas para la tabla final
metricas_rf <- cm_rf$byClass
rf_metrics_vector <- c(
  Modelo        = "Random Forest",
  Precision     = round(mean(metricas_rf[, "Pos Pred Value"], na.rm=T), 4),
  Sensibilidad  = round(mean(metricas_rf[, "Sensitivity"], na.rm=T), 4),
  Especificidad = round(mean(metricas_rf[, "Specificity"], na.rm=T), 4),
  F1_Score      = round(mean(metricas_rf[, "F1"], na.rm=T), 4)
)

print(cm_rf)

# Gráfico de Importancia Top 30
vip <- as.data.frame(caret::varImp(rf_fit)$importance)
vip$Overall <- rowMeans(vip, na.rm=TRUE) # Promedio si hay múltiples clases
vip <- vip %>% tibble::rownames_to_column("gene") %>% dplyr::arrange(dplyr::desc(Overall))
plt_rf <- vip %>% dplyr::slice(1:30) %>%
  ggplot(aes(x = reorder(gene, Overall), y = Overall)) +
  geom_col(fill = "#2a9d8f") + coord_flip() + theme_classic() +
  labs(title="Top 30 Importancia (RF)", x="Gen", y="Importancia")
print(plt_rf)


##-----------------------------------------------------------------------
#                           Método 2: - Naive Bayes
##-----------------------------------------------------------------------

# 1) Usamos los datos globales
nb_train <- train_global
nb_test  <- test_global

# --- FIX 1: Nombres seguros
safe_names <- make.names(colnames(nb_train), unique = TRUE)
colnames(nb_train) <- safe_names
colnames(nb_test)  <- safe_names

# --- FIX 2: Varianza cero
genes_problematicos <- c()
clases <- unique(nb_train$class)
genes_list <- setdiff(names(nb_train), "class")

for (cl in clases) {
  subset_data <- nb_train[nb_train$class == cl, genes_list]
  vars <- apply(subset_data, 2, var)
  genes_problematicos <- c(genes_problematicos, names(vars[vars < 1e-9]))
}
genes_problematicos <- unique(genes_problematicos)

if(length(genes_problematicos) > 0){
  cat("Naive Bayes: Eliminando genes con varianza 0:", length(genes_problematicos), "\n")
  nb_train <- nb_train[, !names(nb_train) %in% genes_problematicos]
  nb_test  <- nb_test[,  !names(nb_test)  %in% genes_problematicos]
}

# 2) Entrenar Naive Bayes
library(klaR)
nb_model <- klaR::NaiveBayes(class ~ ., data = nb_train, usekernel = FALSE, fL = 0)

# 3) Evaluación NB
pred_nb <- suppressWarnings(predict(nb_model, nb_test)$class)
cm_nb <- caret::confusionMatrix(pred_nb, nb_test$class)

metricas_nb <- cm_nb$byClass
nb_metrics_vector <- c(
  Modelo        = "Naive Bayes",
  Precision     = round(mean(metricas_nb[, "Pos Pred Value"], na.rm=T), 4),
  Sensibilidad  = round(mean(metricas_nb[, "Sensitivity"], na.rm=T), 4),
  Especificidad = round(mean(metricas_nb[, "Specificity"], na.rm=T), 4),
  F1_Score      = round(mean(metricas_nb[, "F1"], na.rm=T), 4)
)

print(cm_nb)


##-----------------------------------------------------------------------
#               Método 3: SVM (Support Vector Machine)
##-----------------------------------------------------------------------

library(e1071)

# 1) Usamos los datos globales DIRECTAMENTE
# (SVM maneja bien la varianza cero, no requiere el filtro extra de NB)
svm_train <- train_global
svm_test  <- test_global

# 2) Entrenamiento (Kernel Lineal)
print("Entrenando SVM...")
svm_model <- svm(class ~ ., data = svm_train, kernel = "linear", cost = 1, scale = FALSE)

# 3) Evaluación SVM
svm_pred <- predict(svm_model, newdata = svm_test)
cm_svm <- confusionMatrix(svm_pred, svm_test$class)

metricas_svm <- cm_svm$byClass
svm_metrics_vector <- c(
  Modelo        = "SVM (Kernel Lineal)",
  Precision     = round(mean(metricas_svm[, "Pos Pred Value"], na.rm=T), 4),
  Sensibilidad  = round(mean(metricas_svm[, "Sensitivity"], na.rm=T), 4),
  Especificidad = round(mean(metricas_svm[, "Specificity"], na.rm=T), 4),
  F1_Score      = round(mean(metricas_svm[, "F1"], na.rm=T), 4)
)

print(cm_svm$table)


##-----------------------------------------------------------------------
#                          TABLA COMPARATIVA FINAL
##-----------------------------------------------------------------------

df_rf  <- as.data.frame(t(rf_metrics_vector))
df_nb  <- as.data.frame(t(nb_metrics_vector))
df_svm <- as.data.frame(t(svm_metrics_vector))

tabla_final <- rbind(df_rf, df_nb, df_svm)

print("=== TABLA COMPARATIVA DE MODELOS (80/20) ===")
library(knitr)
kable(tabla_final)

# Guardar
write.csv(tabla_final, "Tabla_Comparativa_Modelos.csv", row.names = FALSE)

# Aunque las tres técnicas demostraron un rendimiento sobresaliente (F1-Score > 0.97 en todos los casos), 
# el SVM con Kernel Lineal resultó ser el modelo superior para este conjunto de datos, logrando una 
# clasificación perfecta sin errores. Esto indica que la frontera de decisión entre los tipos de cáncer 
# es linealmente separable en el espacio de genes proporcionado. 
# Random Forest se posiciona como la segunda mejor opción, ofreciendo una alta interpretabilidad 
# biológica, mientras que Naive Bayes, pese a ser el de menor rendimiento relativo, 
# sigue siendo una opción válida por su eficiencia. 

