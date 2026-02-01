

##--------------------------------------------------------------------
#           Métodos de Clusterización (Aprendizaje No Supervisado)
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



