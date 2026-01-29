
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
  "dplyr"
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
#             Reduccion de dimensionalidad de datos
##--------------------------------------------------------------------

# Porcentage de genes con distribucion normal
(sum(shapiro_table$p_value >= 0.05, na.rm = TRUE)/nrow(shapiro_table))*100

# >>> Se recomienda trabajar con métodos no paramétricos como t-SNE o PCA con datos escalados

# Escalar Datos -- Usar esta data base para futuros analisis
# Preprocesamiento en base de log
datos_escalados <- log2(df_genes + 1) %>% scale()







