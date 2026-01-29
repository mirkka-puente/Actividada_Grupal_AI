
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
colnames(classes) <- c("sample", "class")

# Add los sample numbers and class
gene_expression$sample <- classes$sample  
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
  "gridExtra"    # Combinar gráficos
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

