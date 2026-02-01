##-----------------------------------------------------------------------
#            MÉTODO SUPERVISADO - Naive Bayes
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