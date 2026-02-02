
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
  "randomForest" # Machine Learnig 
)


# Instalación en CRAN
#install.packages(packs, dependencies = TRUE)

# Comprobación automática
lapply(packs, library, character.only = TRUE)

library(glmnet) # ElasticNet
library(tidyverse)
library(caret) # ML
library(rpart) # DT
library(rpart.plot) # DT plot
library(rattle) # DT plot
library(pROC) # ROC
library(PRROC) # PR-Curve
library(MASS) # LDA
library(klaR) # RDA
library(gridExtra) # juntar los gráficos
library(dplyr)
library(Rtsne)#Métodos no supervidados t-SNE
library(randomForest) # Metodo supervisado

# Existe NA en nuestros datos?
cat("Total de NAs:", sum(is.na(gene_expression)))


##--------------------------------------------------------------------
#             Revision de normalidad de gene expression
##--------------------------------------------------------------------

prueba_shapiro<- lapply(column_names, function(g) {
  
  x <- as.numeric(gene_expression[[g]])
  
  if (length(na.omit(x)) < 3 || sd(x, na.rm = TRUE) == 0) {
    return(data.frame(
      gene = g,
      W = NA,
      p_value = NA
    ))
  }
  
  test <- shapiro.test(x)
  
  data.frame(
    gene = g,
    W = test$statistic,
    p_value = test$p.value
  )
})

shapiro_table <- bind_rows(prueba_shapiro) %>%
  mutate(
    normalidad = ifelse(p_value >= 0.05,
                        "Normalidad",
                        "No normalidad")
  ) %>%
  arrange(p_value)

tail(shapiro_table)

#Genes con resultado NA
shapiro_table$gene[(is.na(shapiro_table$p_value))] #No cumple con requisitos de distribucion normal

# Data set sin los genes que no cumplen requisitos y sin las columnas de sampleID y class para analisis
columnas_a_excluir <- c("MIER3", "ZCCHC12", "RPL22L1", "sample_ID", "class")

df_genes <- gene_expression[, !names(gene_expression) %in% columnas_a_excluir]

##--------------------------------------------------------------------
#                      Escala de Datos
##--------------------------------------------------------------------

# Porcentage de genes con distribucion normal
(sum(shapiro_table$p_value >= 0.05, na.rm = TRUE)/nrow(shapiro_table))*100

# >>> Se recomienda trabajar con métodos no paramétricos como t-SNE o PCA con datos escalados

# Escalar Datos -- Usar esta data base para futuros analisis
# Preprocesamiento en base de log
datos_escalados <- log2(df_genes + 1) %>% scale()
datos_escalados <- as.data.frame(datos_escalados)

# Add los sample numbers and class a los datos escalados
datos_escalados$sample_ID <- classes$sample_ID  
datos_escalados$class <- classes$class

##-----------------------------------------------------------------------
#    1. Métodos no supervisados - Reduccion de dimensionalidad de datos 
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
#                        Método 1: PCA
##-----------------------------------------------------------------------

# Excluimos sample_ID y class para el cálculo matemático
columnas_no_numericas <- c("sample_ID", "class")
data_pca <- datos_escalados[, !names(datos_escalados) %in% columnas_no_numericas]

# Calculo PCA
pca.results <- prcomp(data_pca, center = TRUE, scale. = FALSE) 

# Resultados para el grafico.
pca.df <- data.frame(pca.results$x)
pca.df$class <- datos_escalados$class 

# Varianza
varianzas <- pca.results$sdev^2
total.varianza <- sum(varianzas)
varianza.explicada <- varianzas / total.varianza

# Ejes
x_label <- paste0('PC1 (', round(varianza.explicada[1] * 100, 2), '%)')
y_label <- paste0('PC2 (', round(varianza.explicada[2] * 100, 2), '%)')

# 5. Graficar PCA
ggplot(pca.df, aes(x = PC1, y = PC2, color = class)) +
  geom_point(size = 3, alpha = 0.8) + 
  labs(title = 'PCA - Expresion Génica', x = x_label, y = y_label, color = 'Clase') +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "gray90"), 
        plot.title = element_text(hjust = 0.5))

##--------------------------------------------------------------------
#    Método 2: t-Distributed Stochastic Neighbor Embedding (t-SNE)
##--------------------------------------------------------------------

# 1. Eliminacion posibles datos duplicados para t-SNE
datos_tsne_clean <- datos_escalados %>% distinct()

# Diferenciacion matriz numerica y etiquetado.
matrix_tsne <- as.matrix(datos_tsne_clean[, !names(datos_tsne_clean) %in% columnas_no_numericas])
labels_tsne <- datos_tsne_clean$class

set.seed(1234) 

# 3. Ejecución algoritmo t-SNE
tsne_out <- Rtsne(X = matrix_tsne, 
                  dims = 2, 
                  check_duplicates = FALSE) # Ya los limpiamos arriba

# 4. Dataframe para uso del ggplot.
tsne_result <- data.frame(tsne_out$Y)
colnames(tsne_result) <- c("Dim1", "Dim2")
tsne_result$class <- labels_tsne # Añadimos las etiquetas correctas

# 5. Graficar t-SNE
ggplot(tsne_result, aes(x = Dim1, y = Dim2, color = class)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = paste0("t-SNE - Datos expresion génica "), 
                      x = "Dimensión 1", y = "Dimensión 2", color = "Clase") +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "gray90"), 
        plot.title = element_text(hjust = 0.5))
##--------------------------------------------------------------------
#         2.  Métodos no supervisados - Clusterización
##--------------------------------------------------------------------

##--------------------------------------------------------------------
#               Método 1: Clustering no jerárquico (K-means)
##--------------------------------------------------------------------

# Dado que el dataset tiene 5 clases biológicas, usamos k=5
set.seed(1234) 
kmeans_result <- kmeans(data_pca, centers = 5, iter.max = 100, nstart = 25)

# Visualización del clustering 

fviz_cluster(kmeans_result, data = data_pca, 
             geom = "point",
             ellipse.type = "convex",
             main = "Gráfico K-means (k=5)",
             ggtheme = theme_minimal())

# Comparación con las etiquetas reales para validar la agrupación
table(Clúster = kmeans_result$cluster, Clase_Real = datos_escalados$class)


##--------------------------------------------------------------------
#               Método 2: Clustering Jerárquico 
##--------------------------------------------------------------------

# 1. Cálculo de la matriz de distancias 
dist_matrix <- dist(data_pca, method = "euclidean")

# 2. Aplicación del método de Ward.D2 (minimiza la varianza interna)
hclust_ward <- hclust(dist_matrix, method = "ward.D2")

# 3. Representación del Dendrograma 
fviz_dend(hclust_ward, k = 5, 
          cex = 0.5, 
          k_colors = "jco",
          rect = TRUE, 
          main = "Dendrograma de Expresión Génica (Método Ward)",
          xlab = "Índice de Observaciones",
          ylab = "Distancia") + 
  theme_classic()

# 4. Corte del árbol para asignar muestras a grupos
grupos_hclust <- cutree(hclust_ward, k = 5)

# Comparación de resultados
table(Clúster_Jerárquico = grupos_hclust, Clase_Real = datos_escalados$class)



##--------------------------------------------------------------------
#                     3. Métodos supervisados  
##--------------------------------------------------------------------
##--------------------------------------------------------------------
#                       Método 1: Random Forest
##--------------------------------------------------------------------

set.seed(12345)  # Fijar semilla para que resultados sean reproducibles.

# 1) Conjunto de modelado 
#    - Filtramos columnas no numéricas y aseguramos la clase como factor.
columnas_no_numericas <- c("sample_ID", "class")
datos_modelo <- datos_escalados %>%
  dplyr::select(-all_of("sample_ID")) %>%      # dejamos 'class' y los genes
  dplyr::mutate(class = factor(class))

# 2) Partición estratificada 80/20
idx <- caret::createDataPartition(datos_modelo$class, p = 0.80, list = FALSE)
train_df <- datos_modelo[idx, , drop = FALSE]
test_df  <- datos_modelo[-idx, , drop = FALSE]

# 3) Eliminar predictores de varianza ~0 usando SOLO el train, identifica genes con muy poca variación o muy pocos valores únicos.
nzv_idx <- caret::nearZeroVar(train_df %>% dplyr::select(-class))
if (length(nzv_idx) > 0) {
  genes_keep <- setdiff(colnames(train_df), c("class"))[ -nzv_idx ]
  train_df <- dplyr::select(train_df, dplyr::all_of(c("class", genes_keep)))
  test_df  <- dplyr::select(test_df,  dplyr::all_of(c("class", genes_keep)))
}

# 4) (Opcional) Imputación mediana si hubiera NA. Para blindar el flujo por si aparece algún NA.
pp <- caret::preProcess(train_df %>% dplyr::select(-class), method = c("medianImpute"))
train_x <- predict(pp, train_df %>% dplyr::select(-class))
test_x  <- predict(pp, test_df  %>% dplyr::select(-class))
train_y <- train_df$class
test_y  <- test_df$class

# 5) Entrenamiento Random Forest con CV 5x2 y búsqueda simple de mtry
p <- ncol(train_x)
mtry_grid <- unique(pmax(1, round(c(sqrt(p)/2, sqrt(p), sqrt(p)*2))))

ctrl <- caret::trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 2,
  classProbs = FALSE,
  verboseIter = TRUE,
  allowParallel = TRUE
)

set.seed(12345)
rf_fit <- caret::train(
  x = train_x, y = train_y,
  method = "rf",
  trControl = ctrl,
  tuneGrid = data.frame(mtry = mtry_grid),
  ntree = 1000,
  importance = TRUE,
  metric = "Accuracy"
)

print(rf_fit)
cat("Mejor mtry seleccionado:", rf_fit$bestTune$mtry, "\n")

# 6) Evaluación en test: Matriz de confusión y métricas.
pred_test <- predict(rf_fit, newdata = test_x)
cm <- caret::confusionMatrix(pred_test, test_y)
cm

# Extraer métricas pedidas
accuracy <- unname(cm$overall["Accuracy"])
byClass <- as.data.frame(cm$byClass)

# F1 por clase = 2 * (PPV * Sens) / (PPV + Sens)
if ("Pos Pred Value" %in% colnames(byClass)) {
  byClass$F1 <- 2 * (byClass$Sensitivity * byClass$`Pos Pred Value`) /
    pmax(byClass$Sensitivity + byClass$`Pos Pred Value`, .Machine$double.eps)
} else {
  byClass$F1 <- NA_real_
}

macro_metrics <- c(
  Macro_Sensitivity = mean(byClass$Sensitivity, na.rm = TRUE),
  Macro_Specificity = mean(byClass$Specificity, na.rm = TRUE),
  Macro_Precision   = mean(byClass$`Pos Pred Value`, na.rm = TRUE),
  Macro_F1          = mean(byClass$F1, na.rm = TRUE)
)

cat("\n===== RESULTADOS TEST (Random Forest) =====\n")
cat(sprintf("Accuracy global: %.4f\n", accuracy))
print(round(macro_metrics, 4))
cat("\nMétricas por clase:\n")
print(round(byClass[, c("Sensitivity","Specificity","Pos Pred Value","F1")], 4))

# 7) Guardar salidas a disco
readr::write_csv(as.data.frame(cm$table), "rf_confusion_matrix.csv")
readr::write_csv(
  tibble::tibble(Metric = names(macro_metrics), Value = as.numeric(macro_metrics)),
  "rf_macro_metrics.csv"
)

# 8) Importancia de genes y gráfico Top 30 (versión robusta)
vip_raw <- caret::varImp(rf_fit)$importance

# Asegurar que es data.frame
vip_tbl <- as.data.frame(vip_raw)

if (!"Overall" %in% colnames(vip_tbl)) {
  if (ncol(vip_tbl) > 1) {
    # varias columnas (p. ej., una por clase) -> promedio como Overall
    vip_tbl$Overall <- rowMeans(vip_tbl, na.rm = TRUE)
  } else {
    # una sola columna 
    only_col <- colnames(vip_tbl)[1]
    vip_tbl$Overall <- vip_tbl[[only_col]]
  }
}

vip <- vip_tbl %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::arrange(dplyr::desc(Overall))

# Guardar tabla completa de importancias
readr::write_csv(vip, "rf_feature_importance.csv")

# Gráfico Top-30
topN <- 30
plt_rf <- vip %>%
  dplyr::slice(1:topN) %>%
  ggplot2::ggplot(ggplot2::aes(x = reorder(gene, Overall), y = Overall)) +
  ggplot2::geom_col(fill = "#2a9d8f") +
  ggplot2::coord_flip() +
  ggplot2::labs(title = sprintf("Top %d genes más importantes (Random Forest)", topN),
                x = "Gen", y = "Importancia (Overall)") +
  ggplot2::theme_classic(base_size = 12)
ggplot2::ggsave("rf_top_genes_importance.png", plt_rf, width = 7, height = 6, dpi = 300)

plt_rf


# 9) Guardar el modelo
saveRDS(rf_fit, file = "rf_model.rds")

##-----------------------------------------------------------------------
#            Método 2: - Naive Bayes
##-----------------------------------------------------------------------

# Usamos los datos ya escalados (datos procesados/limpios)
data_sup <- datos_escalados

# La clase debe ser factor (clasificación)
data_sup$class <- as.factor(data_sup$class)

# Quitamos sample_ID porque no es predictor
data_sup <- data_sup[, !names(data_sup) %in% c("sample_ID")]

# Dividir en entrenamiento (80%) y prueba (20%)
set.seed(1995)
trainIndex <- createDataPartition(data_sup$class, p = 0.8, list = FALSE)
trainData <- data_sup[trainIndex, ]
testData  <- data_sup[-trainIndex, ]

# Quitar genes que tienen varianza 0 dentro de alguna clase (pueden romper NB)
y <- trainData$class
X <- trainData[, setdiff(names(trainData), "class"), drop = FALSE]
X <- X[, sapply(X, is.numeric), drop = FALSE]

keep <- sapply(X, function(v) {
  vars_por_clase <- tapply(v, y, function(z) var(z, na.rm = TRUE))
  all(vars_por_clase > 0, na.rm = TRUE)
})

cat("Genes eliminados por varianza 0 en alguna clase:", sum(!keep), "\n")

X_train <- X[, keep, drop = FALSE]

# Hacer nombres seguros (algunos genes tienen caracteres raros)
orig_names <- colnames(X_train)
safe_names <- make.names(orig_names, unique = TRUE)

colnames(X_train) <- safe_names
X_test <- testData[, orig_names, drop = FALSE]
colnames(X_test) <- safe_names

# Entrenar Naive Bayes
library(klaR)
library(caret)

nb_model <- klaR::NaiveBayes(
  x = X_train,
  grouping = y,
  usekernel = FALSE,
  fL = 0
)

# Predicción (se suprimen warnings porque NB en alta dimensión suele generar avisos)
pred_nb <- suppressWarnings(predict(nb_model, X_test)$class)

# Matriz de confusión
cm_nb <- confusionMatrix(pred_nb, testData$class)
cm_nb

# Métricas pedidas (por clase, porque es multiclase)
byc <- cm_nb$byClass

metrics_nb <- data.frame(
  Clase = rownames(byc),
  Precision = byc[, "Pos Pred Value"],
  Sensibilidad = byc[, "Sensitivity"],
  Especificidad = byc[, "Specificity"],
  F1 = byc[, "F1"],
  row.names = NULL
)

metrics_nb

##-----------------------------------------------------------------------
#   Método 3: SVM (Support Vector Machine) + MÉTRICAS
##-----------------------------------------------------------------------

#SVM es excelente encontrando el "hiperplano" (la línea divisoria) óptima en espacios con muchas dimensiones
#Busca maximizar el margen de separación entre clases, lo que lo hace robusto y generaliza bien.
#Nota: Es obligatorio escalar los datos (StandardScaler)

# Carga de librería específica para SVM
if(!require(e1071)) install.packages("e1071")
library(e1071)

# Preparación del Dataset para el Modelo
# Se elimina la columna sample_ID ya que no aporta información biológica predictiva
# y podría causar sobreajuste ya que solo memorizaría los IDs en lugar de aprender de los genes.
df_model <- datos_escalados[, !names(datos_escalados) %in% c("sample_ID")]

# División de Datos (Train / Test)
# Se usa una partición 70% Entrenamiento - 30% Prueba
set.seed(1995) # Semilla para reproducibilidad
trainIndex <- createDataPartition(df_model$class, p = 0.7, list = FALSE)

train_data <- df_model[trainIndex, ]
test_data  <- df_model[-trainIndex, ]

cat("\nDimensiones del set de Entrenamiento:", dim(train_data))
cat("\nDimensiones del set de Prueba:", dim(test_data), "\n")

# Entrenamiento del Modelo SVM
# Se utiliza un Kernel Lineal, ideal para datos de alta dimensionalidad (muchos genes)
print("Entrenando modelo SVM")
svm_model <- svm(class ~ ., 
                 data = train_data, 
                 kernel = "linear", 
                 cost = 1,      # Penalización estándar
                 scale = FALSE) # Los datos ya vienen escalados de la fase anterior


# Predicción y Evaluación
# Predicción sobre datos nuevos (Test set)
svm_pred <- predict(svm_model, newdata = test_data)

# Generación de la Matriz de Confusión
# Se asume que la clase de interés (positiva) es la primera o la patológica.
# caret detecta automáticamente los niveles, pero se puede forzar 'positive'
cm_svm <- confusionMatrix(data = svm_pred, 
                          reference = test_data$class)

# Visualización de Resultados
cat("\n-------------------------------------------")
cat("\n RESULTADOS DEL MODELO SVM (SUPERVISADO)")
cat("\n-------------------------------------------\n")

cm_svm


cat("\n-------------------------------------------")
cat("\n MATRIZ DE CONFUSIÓN DEL MODELO SVM (SUPERVISADO)")
cat("\n-------------------------------------------\n")

print(cm_svm$table)

# Extracción de métricas para Múltiples Clases
metricas <- cm_svm$byClass

cat("\n--- Métricas de Desempeño (Promedio Global) ---\n")

# Como metricas es una matriz, accedemos a la columna [ , "NombreColumna"]
# y calculamos la media (mean) para tener un solo número representativo.

# Precisión
precision_global <- mean(metricas[, "Precision"], na.rm = TRUE)
cat("Precisión:      ", round(precision_global, 4), "\n")

# Sensibilidad
sensibilidad_global <- mean(metricas[, "Sensitivity"], na.rm = TRUE)
cat("Sensibilidad:   ", round(sensibilidad_global, 4), "\n")

# Especificidad
especificidad_global <- mean(metricas[, "Specificity"], na.rm = TRUE)
cat("Especificidad:  ", round(especificidad_global, 4), "\n")

# F1-Score
f1_global <- mean(metricas[, "F1"], na.rm = TRUE)
cat("F1-Score:       ", round(f1_global, 4), "\n")

# Interpretación:
# - F1-Score cercano a 1 indica un balance excelente entre precisión y sensibilidad.
# - SVM Lineal suele separar muy bien clases en datos de expresión génica. 


# --- CREACIÓN DE LA TABLA COMPARATIVA FINAL (SUPERVISADO) ---

#  Definimos los resultados del modelo SVM
# Métricas globales salieron 1.0000, las ponemos directas.
modelo_svm <- c(
  Modelo        = "SVM (Kernel Lineal)",
  Precision     = round(precision_global, 4),    
  Sensibilidad  = round(sensibilidad_global, 4),
  Especificidad = round(especificidad_global, 4),
  F1_Score      = round(f1_global, 4)
)

###AGREGAR SUS MODELOS
#modelo_rf <- c(
#Modelo        = "Random Forest",
#Precision     = NA, 
#Sensibilidad  = NA, 
#Especificidad = NA, 
#F1_Score      = NA
#)

#modelo_otro <- c(
#  Modelo        = "KNN / Otro",  # Cambiar por el nombre real
#  Precision     = NA, 
#  Sensibilidad  = NA, 
#  Especificidad = NA, 
#  F1_Score      = NA
#)

#Unir todo en un solo Dataframe
#tabla_comparativa <- rbind(modelo_svm, modelo_rf, modelo_otro)
#tabla_comparativa <- as.data.frame(tabla_comparativa)
tabla_comparativa <- as.data.frame(modelo_svm)

# Fin del script  
