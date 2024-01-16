
### Librerías

library(readr)
library(data.table)
library(dplyr)
library(tidyr)
library(caret)
library(tibble)
library(reshape)


# Ruta al archivo TSV con identificadores y datos clinicos
ruta_tsv <-  "C:/Yo/MÁSTER/UOC/TFM/rabajo/RNAseq KIRC KIRP KICH/carpetas/gdc_sample_sheet.2023-11-10.tsv"
ruta_carpeta <- "C:/Yo/MÁSTER/UOC/TFM/rabajo/RNAseq KIRC KIRP KICH/carpetas/"
ruta_clinical <-
  'C:/Yo/MÁSTER/UOC/TFM/rabajo/RNAseq KIRC KIRP KICH/carpetas/clinical.tsv'
ruta_sample <-
  'C:/Yo/MÁSTER/UOC/TFM/rabajo/RNAseq KIRC KIRP KICH/carpetas/sample.tsv'


# Leer solo las columnas especificas del archivo TSV

columnas_labels <- c("case_submitter_id", "primary_diagnosis",
                     "ajcc_pathologic_stage")

columnas_tsv <- c("File ID", "File Name", "Case ID", "Sample ID")
columnas_sample <- c("sample_submitter_id", "tissue_type")

## Leer TSV sobre muestra
sample_tsv <- fread(ruta_sample, select = columnas_sample, header = TRUE, sep = "\t")
names(sample_tsv)[names(sample_tsv) == "sample_submitter_id"] <- "Sample ID"

### Datos RNA
datos_tsv <- fread(ruta_tsv, select = columnas_tsv, header = TRUE, sep = "\t")

### Datos LABELS. Eliminando Stage desconocido y Tumor desconocido, uniendo a Tipo Tumor
clinical_tsv <- fread(ruta_clinical, select = columnas_labels, header = TRUE, sep = "\t")

clinical_tsv <- clinical_tsv[!clinical_tsv$ajcc_pathologic_stage == "'--", ]

clinical_tsv <- clinical_tsv[!clinical_tsv$primary_diagnosis == "Renal cell carcinoma, NOS", ]


# Cambiar los Stage I II III IV a Early y Late
clinical_tsv$ajcc_pathologic_stage  <- as.factor(clinical_tsv$ajcc_pathologic_stage)

levels(clinical_tsv$ajcc_pathologic_stage) <- c('Early', 'Early', 'Late', 'Late')


# Función para leer los datos de los archivos de texto
leer_datos_texto <- function(File_ID, File_Name, ruta_carpeta) {
  ruta_archivo <- file.path(ruta_carpeta, paste0(File_ID, "/", File_Name))

  # Verificar si el archivo de texto existe
  if (file.exists(ruta_archivo)) {
    # Leer solo las columnas especificas del archivo de texto
    columnas_txt <- c("gene_id", "unstranded")
    # Omitimos Na values, nos saltamos una fila porque es no informativa, los
    # nombres de preguntas estan en la segunda fila
    datos_texto <- na.omit(fread(ruta_archivo, skip = 1, select = columnas_txt,
                                 header = TRUE, sep = "\t"))
    datos_texto <- t(datos_texto)
    colnames(datos_texto) <- unlist(datos_texto[row.names(datos_texto)=="gene_id",])

    datos_texto <- datos_texto[!row.names(datos_texto)=="gene_id",]
  } else {
    # Crear un data.table con un mensaje si el archivo no existe
    datos_texto <- data.table("Sin datos de archivo de texto")
  }

  return(datos_texto)
}


# Aplicar la función a cada fila del DataFrame del archivo TSV
# que nos dice nombre de carpeta y nombre de archivo
datos_completos <- mapply(leer_datos_texto, datos_tsv$`File ID`, datos_tsv$`File Name`,
                          ruta_carpeta)

# Unimos todo en un dataframe, haciendo la traspuesta de datos_completos
datos_tsv <- cbind(datos_tsv, t(datos_completos))

datos_SAMPLE <- merge(datos_tsv, sample_tsv, by="Sample ID",all.x=TRUE)
table(is.na(datos_SAMPLE))
# No hay datos perdidos
datos_tsv <- datos_SAMPLE
table(datos_tsv$tissue_type)
# Normal 129 Tumor 899

### unir con clinical. Ahora tendremos columnas Diagnostico, Estadio, Tumor/Control
names(clinical_tsv)[names(clinical_tsv) == "case_submitter_id"] <- "Case ID"
datos_tsv <- inner_join(datos_tsv, clinical_tsv, by="Case ID",multiple='any')
table(datos_tsv$tissue_type)
# Normal 129 Tumor 852

table(is.na(datos_tsv))
# NO hay datos faltantes.

# Eliminamos las columnas File ID, File Name, Case ID, N_unmapped, N_multimapping,
# N_noFeature, N_ambiguous
datos_tsv <- datos_tsv[,-c(2,3,4,5,6,7,8)]
tissue_type <- datos_tsv$tissue_type


## Unir en la misma columna subtipo tumor + estadio
### Esto se puede hacer antes o despues del DESeq2, depende de si queremos
## que compare entre tipo de tejido (tumor/control) o entre estado early/late
## en este caso como quería encontrar genes expresados diferencialmente entre
## early y late, lo hice después del DESeq2. Sin embargo, después de hacer el
## DESeq2, es más cómodo hacer todos estos pasos de nuevo y reorganizar los datos
## desde aquí.
## El código está puesto según la segunda vuelta (cuando ya habia hecho DESeq2)
# Para hacer DESeq2 simplemente se hace el unite() más tarde, una vez hecho
datos_tsv2 <- datos_tsv
datos_tsv <- unite(datos_tsv, Diagnosis, c("primary_diagnosis", "ajcc_pathologic_stage"))

# Hacemos que la columna Sample ID tenga valores unicos para hacerla rownames
#datos_tsv[,1] <- lapply(datos_tsv[,1], function(x) make.unique(as.character(x)))
#datos_tsv[,1] <- lapply(datos_tsv, function(x) make.unique(as.character(x)))
datos_tsv$`Sample ID` <- make.unique(as.character(datos_tsv[,1]))

rownames(datos_tsv) <- datos_tsv$`Sample ID`
datos_tsv <- datos_tsv[,-1]
# Si no se transformaron las columnas a Diagnosis con unite(), sustituir 2 por 3.
datos_tsv[,1:(ncol(datos_tsv)-2)] <-
  apply(datos_tsv[,1:(ncol(datos_tsv)-2)],2,as.numeric)

### Solo si se hizo unite()
# Cambiar los nombres de los niveles de Diagnosis para mas facil visualizacion
diagnosis <- datos_tsv$Diagnosis
datos_tsv$Diagnosis <- factor(datos_tsv$Diagnosis,levels=unique(datos_tsv$Diagnosis))
levels(datos_tsv$Diagnosis)[levels(datos_tsv$Diagnosis)==
                              "Papillary adenocarcinoma, NOS_Early"] <- "pRCC, Early"
levels(datos_tsv$Diagnosis)[levels(datos_tsv$Diagnosis)==
                              "Papillary adenocarcinoma, NOS_Late"] <- "pRCC, Late"

levels(datos_tsv$Diagnosis)[levels(datos_tsv$Diagnosis)==
                              "Clear cell adenocarcinoma, NOS_Early"] <- "ccRCC, Early"
levels(datos_tsv$Diagnosis)[levels(datos_tsv$Diagnosis)==
                              "Clear cell adenocarcinoma, NOS_Late"] <- "ccRCC, Late"

levels(datos_tsv$Diagnosis)[levels(datos_tsv$Diagnosis)==
                              "Renal cell carcinoma, chromophobe type_Early"] <- "chRCC, Early"
levels(datos_tsv$Diagnosis)[levels(datos_tsv$Diagnosis)==
                              "Renal cell carcinoma, chromophobe type_Late"] <- "chRCC, Late"


# Hacemos los boxplots para ver cómo se distribuyen los datos por clase
# Reorganizar los datos para graficar los boxplots
datos_melted <- gather(datos_tsv, key = "Gen", value = "Expresion", -c(Diagnosis, tissue_type))

# Crear el boxplot con ggplot2
ggplot(datos_melted, aes(x = Diagnosis, y = Expresion, fill = Diagnosis)) +
  geom_boxplot() +
  labs(title = "Boxplot de Expresion Genica por Diagnostico",
       x = "Diagnostico", y = "Expresion Genica") +
  theme_minimal()

#Los datos tal como están son ilegibles. Hace fala una transformación logarítmica ya que las expresiones tienen rangos muy grandes.

#Transformacion logaritmica
datos_melted[,4] <- log2(datos_melted[,4]+1)


# Crear el boxplot con ggplot2
ggplot(datos_melted, aes(x = Diagnosis, y = Expresion, fill = Diagnosis)) +
  geom_boxplot() +
  labs(title = "Boxplot de Expresion Genica por Diagnostico",
       x = "Diagnostico", y = "Expresion Genica") +
  theme_minimal()


# Preparación de los datos para DESeq2
#colData son las dos ultimas columnas de datos_tsv. Depende de si se hizo unite
# o no, serán 2 o 3, porque estarán tipo de tumor + estado unidas o no.
colData <- datos_tsv[,(ncol(datos_tsv)-1):ncol(datos_tsv)]
countData <- as.data.frame(t(datos_tsv[,1:(ncol(datos_tsv)-2)]))
countData <- countData[ , order(names(countData))]
colData <- colData[ order(rownames(colData)) , ]
head(countData[,1:3],3)
head(colData,5)
# Los nombres de filas de colData deben ser los mismos que los nombres de
# columnas de countData
table(rownames(colData)==colnames(countData))


##### DESeq2
library(DESeq2)
library(tibble)

colData2 <- tibble::rownames_to_column(colData, "Sample_ID")
#countData3 <- tibble::rownames_to_column(countData2, "Sample_ID")
#countData3 <- countData3 %>% mutate_at(2:ncol(countData3), as.numeric)
#countData2 <- as.matrix(countData2)
#colData <- as.matrix(colData)
colData2 <- data.frame(colData2, row.names=NULL)
countData3 <- countData
dds <- DESeqDataSetFromMatrix(countData=countData3,
                              colData=colData2,
                              # diseño tipo tejido + Diagnosis
                              # design=~tissue_type+Diagnosis, tidy = TRUE)
                              # diseño estado tumor + tipo tumor
                 design=~ajcc_pathologic_stage+primary_diagnosis, tidy = TRUE)
                  # diseño estado tumor + tipo tumor + tipo_tejido
                  #design=~ajcc_pathologic_stage+primary_diagnosis+tissue_type, tidy = TRUE)



####  Se obtienen los resultados del analisis de expresion
## por defecto te hace la comparacion Estado
res <- results(dds)
head(results(dds, tidy=TRUE))

res <- res[order(res$padj),]
head(res)

res_filtered <- subset(res, padj < 0.05)
head(res_filtered)

# Análisis gráfico de los resultados obtenidos

### Volcano Plot

with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot",
               xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20,
                                    col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange,
                                                           -log10(pvalue), pch=20, col="red"))

### PCA

vsdata <- vst(dds, blind=FALSE)

plotPCA(vsdata, intgroup=c("ajcc_pathologic_stage", "primary_diagnosis"))

# Seleccionamos los genes obtenidos tras el analisis

selected_genes <- rownames(res_filtered[1])
colData3 <- colData2
rownames(colData3) <- colData3[,1]
colData3 <- colData3[,-1]
countData4 <- t(countData3)

colnames(countData4) <- countData4[1,]
countData4<- countData4[-1,]

countData4 <- select(countData4, all_of(selected_genes))
datos_final <- cbind(data.frame(countData4),colData2)

# si no se hizo el unite antes...
#datos_final <- unite(datos_final, Diagnosis, c("primary_diagnosis", "ajcc_pathologic_stage"))

table(datos_final$Diagnosis)


# Preprocesamiento


### Hay que transformar las labels- Las que sean en control, que digan control
# en Diagnosis.

datos_final <- within(datos_final, Diagnosis[tissue_type=='Normal']<-'Control')
datos_final <- datos_final[,!names(datos_final)==c('tissue_type','Sample_ID')]
#### Tenemos 15498 Variables. 1 es Diagnosis. 15498 genes.
labels <- datos_final[,ncol(datos_final)]
datos_final <- datos_final %>% mutate_at(1:(ncol(datos_final)-1), as.numeric)

table(datos_final$Diagnosis)


# Los datos no estan balanceados así que se balancean
data_new <- data.frame()
clases <- c('pRCC, Late','pRCC, Early','chRCC, Late', 'chRCC, Early',
            'ccRCC, Early','ccRCC, Late', 'Control')

library(scutr)

for (i in 1:length(clases)) {
  # Realizar balanceo con SMOTE
  datos_balanceados <- oversample_smote(data=datos_final, cls=clases[i],
                                        # cls_col='Diagnosis',m=400,k=NA)
                                        cls_col='Diagnosis',m=320,k=NA)
  print(table(datos_balanceados$Diagnosis))
  data_new <- rbind(data_new,datos_balanceados)
}
# VCuantas clases hay despues del balanceo
table(data_new$Diagnosis)
saveRDS(data_new,'datos_balanceados.Rda')
datos_final <- readRDS(file='datos_balanceados.Rda')
# Primero se hizo SMOTE con 320 muestras, para que todas tengan como la mayoritaria
# Luego se hizo que todas tuvieran 400.

### Aplicar transformacion logaritmica porque son datos de expresion
labels <- datos_final$Diagnosis
datos_final <- log2(datos_final[,1:(ncol(datos_final)-1)]+1)

#data <- cbind(datos_final, labels)

## Aqui fue mi error. Deberia haber hecho esta partición antes de balancear.

seed <- 144
perc <- 0.8
set.seed(seed)

# Con sample sacamos los indices con los que obtenemos los datos de train
train_indices <- sample(nrow(datos_final), floor(perc * nrow(datos_final)))
data_train <- datos_final[train_indices, 1:(ncol(datos_final))]  # Datos para train
data_test <- datos_final[-train_indices, 1:(ncol(datos_final))]  # Datos para test

#Etiquetas
data_train_labels <- labels[train_indices]
data_test_labels <- labels[-train_indices]

#### Normalizar

maximo <- max(data_train)
minimo <- min(data_train)


saveRDS(maximo, file="maximo_train.Rda")
saveRDS(minimo, file="minimo_train.Rda")

# Función para normalizar
normalize <- function(x) {
  return((x - minimo) / (maximo - minimo))
}

# Aplicamos la función
data_train <- as.data.frame(lapply(data_train, normalize))
data_test <- as.data.frame(lapply(data_test, normalize))
table(is.na(data_train))
# No hay datos faltantes






######################################
#####################MODELOS



library(ROCR)  # para el KNN
library(neuralnet) # Para Redes Neuronales
library(kernlab)  # Para el SVM
library(C50) # Decision Tree
library(randomForest) # Random Forest


####################### KNN ################################

# valores de k que usaremos en el algoritmo KNN
k_values <- c(1, 3, 5, 7, 11)


resultados_knn <- data.frame('Accuracy'=numeric(), 'Kappa' =numeric(),
                             'Sensitivity'=numeric(),'Specificity' = numeric())


for (k in k_values) {
  # ajustar el modelo para cada k en el loop for
  modelo_knn <- class::knn(train = data_train, test = data_test,
                           cl = as.factor(data_train_labels), k=k, prob=FALSE)

  # Metricas
  matriz <- confusionMatrix(modelo_knn, as.factor(data_test_labels))
  Accuracy <- matriz$overall["Accuracy"]
  Kappa <- matriz$overall["Kappa"]
  Sensitivity <- mean(matriz$byClass[,1])
  Specificity <- mean(matriz$byClass[,2])

  resultados_knn <- rbind(resultados_knn, data.frame(row.names=paste0('k',k),
                                                     Accuracy=Accuracy,
                                                     Kappa=Kappa,
                                                     Sensitivity=Sensitivity,
                                                     Specificity=Specificity))

}

print(resultados_knn)

saveRDS(resultados_knn, file="resultados_knn.Rda")


######################### NEURAL NETWORK ################################



data_train_2 <- cbind(Diagnosis = data_train_labels, data_train)

formula <- as.formula("Diagnosis ~ .")

modelo_nn_1 <- neuralnet(formula, data = data_train_2, hidden = c(8))
modelo_nn_2 <- neuralnet(formula, data = data_train_2, hidden = c(10,5))

saveRDS(modelo_nn_1, file="modelo_nn_1.Rda")
saveRDS(modelo_nn_2, file="modelo_nn_2.Rda")

modelo_nn_1 <- readRDS(file="modelo_nn_1.Rda")
modelo_nn_2 <- readRDS(file="modelo_nn_2.Rda")


resultados_nn <- data.frame('Accuracy'=numeric(), 'Kappa' =numeric(),
                            'Sensitivity'=numeric(),'Specificity' = numeric())

for (i in 1:2) {
  modelo_nn_results <- compute(get(paste0('modelo_nn_',i)), data_test)$net.result

  predict_nn <- as.factor(apply(modelo_nn_results, 1,
                                function(x) levels(as.factor(labels))[which.max(x)]))

  # Calcular la precisión del modelo
  matriz <- confusionMatrix(predict_nn, as.factor(data_test_labels))
  Accuracy <- matriz$overall["Accuracy"]
  Kappa <- matriz$overall["Kappa"]
  Sensitivity <- mean(matriz$byClass[,1])
  Specificity <- mean(matriz$byClass[,2])


  resultados_nn <- rbind(resultados_nn, data.frame(row.names=paste0('Modelo ', i),
                                                   Accuracy=Accuracy,
                                                   Kappa=Kappa,
                                                   Sensitivity=Sensitivity,
                                                   Specificity=Specificity))

}

print(resultados_nn)



saveRDS(resultados_nn, file="resultados_nn.Rda")






########################### SVM ###########################################
seed <- 144

kernel <- c("vanilladot", "rbfdot")


resultados_SVM <- data.frame('Accuracy'=numeric(), 'Kappa' =numeric(),
                             'Sensitivity'=numeric(),'Specificity' = numeric(),
                             'Media'=numeric())


for (i in kernel) {

  set.seed(seed)
  # Modelo de SVM con función RBF
  modelo_SVM <- ksvm(as.factor(data_train_labels) ~ .,
                     data = data_train, kernel = i)

  predict_SVM <- predict(modelo_SVM, data_test)
  # Valores de interés
  matriz <- confusionMatrix(predict_SVM, as.factor(data_test_labels))
  Accuracy <- matriz$overall["Accuracy"]
  Kappa <- matriz$overall["Kappa"]
  Sensitivity <- mean(matriz$byClass[,1])
  Specificity <- mean(matriz$byClass[,2])



  resultados_SVM <- rbind(resultados_SVM, data.frame(row.names=i,
                                                     Accuracy=Accuracy,
                                                     Kappa=Kappa,
                                                     Sensitivity=Sensitivity,
                                                     Specificity=Specificity,
                                                     Media=mean(Accuracy,Kappa,
                                                                Sensitivity,
                                                                Specificity)))


}

# Imprimir la tabla de resultados
print(resultados_SVM)


saveRDS(resultados_SVM, file="resultados_SVM.Rda")





####################################################################
####################################################################
## Naive Bayes.
#data_train_labels <- data_train[,1]
#data_train <- data_train[,2:ncol(data_train)]


resultados_NB <- data.frame('Accuracy'=numeric(), 'Kappa' =numeric(),
                            'Sensitivity'=numeric(),'Specificity' = numeric(),
                            'Media'=numeric())
lapl <- c(0)

for (i in lapl) {
  # Crear el clasificador Naive Bayes
  modelo_NB <- naiveBayes(data_train, as.factor(data_train_labels), laplace=i)

  # Realizar predicciones en el conjunto de prueba
  predict_Bayes <- predict(modelo_NB, data_test)

  # Calcular la precisión del modelo
  matriz <- confusionMatrix(predict_Bayes, as.factor(data_test_labels))
  Accuracy <- matriz$overall["Accuracy"]
  Kappa <- matriz$overall["Kappa"]
  Sensitivity <- mean(matriz$byClass[,1])
  Specificity <- mean(matriz$byClass[,2])

  resultados_NB <- rbind(resultados_NB, data.frame(row.names=paste0('Laplace: ',i),
                                                   Accuracy=Accuracy,
                                                   Kappa=Kappa,
                                                   Sensitivity=Sensitivity,
                                                   Specificity=Specificity,
                                                   Media=mean(Accuracy,Kappa,
                                                              Sensitivity,
                                                              Specificity)))
}

print(resultados_NB)


saveRDS(resultados_NB, file="resultados_NB.Rda")





## Árbol de Clasificación.


boost <- c(1,10)


resultados_C50 <- data.frame('Accuracy'=numeric(), 'Kappa' =numeric(),
                             'Sensitivity'=numeric(),'Specificity' = numeric(),
                             'Media'=numeric())
for (i in boost) {
  modelo_C50 <- C5.0(data_train, as.factor(data_train_labels), trials=i,costs=NULL)

  predict_C50 <- predict(modelo_C50, data_test, type='class')

  # Calcular la precisión del modelo
  matriz <- confusionMatrix(predict_C50, as.factor(data_test_labels))
  Accuracy <- matriz$overall["Accuracy"]
  Kappa <- matriz$overall["Kappa"]
  Sensitivity <- mean(matriz$byClass[,1])
  Specificity <- mean(matriz$byClass[,2])

  resultados_C50 <- rbind(resultados_C50, data.frame(row.names=paste0('Boost: ',i),
                                                     Accuracy=Accuracy,
                                                     Kappa=Kappa,
                                                     Sensitivity=Sensitivity,
                                                     Specificity=Specificity,
                                                     Media=mean(Accuracy,Kappa,
                                                                Sensitivity,
                                                                Specificity)))


}

print(resultados_C50)



saveRDS(resultados_C50, file="resultados_C50.Rda")










##### Random Forest


ntree <- c(100,200,500)


resultados_randomF <- data.frame('Accuracy'=numeric(), 'Kappa' =numeric(),
                                 'Sensitivity'=numeric(),'Specificity' = numeric())


for (i in ntree) {

  modelo_randomF <- randomForest(data_train, as.factor(data_train_labels), ntree=i)

  predict_randomF <- predict(modelo_randomF, data_test, type='class')





  # Calcular la precisión del modelo
  matriz <- confusionMatrix(predict_randomF, as.factor(data_test_labels))
  Accuracy <- matriz$overall["Accuracy"]
  Kappa <- matriz$overall["Kappa"]
  Sensitivity <- mean(matriz$byClass[,1])
  Specificity <- mean(matriz$byClass[,2])

  resultados_randomF <- rbind(resultados_randomF,
                              data.frame(row.names=paste0('N Arboles: ', i),
                                         Accuracy=Accuracy,
                                         Kappa=Kappa,
                                         Sensitivity=Sensitivity,
                                         Specificity=Specificity))

}

print(resultados_randomF)



saveRDS(resultados_randomF, file="resultados_randomF.Rda")






modelos <- c('k-Nearest Neighbour', 'Naive Bayes', 'Artificial Neural Network',
             'Support Vector Machine', 'Árbol de Clasificación', 'Random Forest')
media_knn <- max(apply(resultados_knn, 1, mean))
media_NB <- max(apply(resultados_NB, 1, mean))
media_nn <- max(apply(resultados_nn, 1, mean))
media_SVM <- max(apply(resultados_SVM, 1, mean))
media_C50 <- max(apply(resultados_C50, 1, mean))
media_randomF <- max(apply(resultados_randomF, 1, mean))

medias <- list()
medias <- rbind(media_knn, media_NB, media_nn, media_SVM, media_C50, media_randomF)


mejor_modelo <- which.max(medias)
mejor_modelo


##### MODELOS CON CARET


### MODELOS CARET
data_train_2 <- cbind(Diagnosis = data_train_labels, data_train)
#modelos <- c('multinom','neuralnet')
modelos <- c('knn','svmLinear','svmRadial','rpart','rf','naive_bayes')
resultados_caret <- data.frame('Accuracy'=numeric(), 'Kappa' =numeric(),
                               'Sensitivity'=numeric(),'Specificity' = numeric(),
                               'Media'=numeric())
control <- trainControl(method = "cv", number = 5, allowParallel = TRUE)  # 5-fold cross-validation

for (i in modelos) {

  modelo <- train(Diagnosis ~ ., data = data_train_2, method = modelo,
                  trControl = control)

  # Ver resultados
  print(modelo)

  predict <- predict(modelo, data_test)
  matriz <- confusionMatrix(predict, as.factor(data_test_labels))

  Accuracy <- matriz$overall["Accuracy"]
  Kappa <- matriz$overall["Kappa"]
  Sensitivity <- mean(matriz$byClass[,1])
  Specificity <- mean(matriz$byClass[,2])

  resultados_caret <- rbind(resultados_nb,
                            data.frame(row.names=paste0('Modelo: ', 'Naive Bayes'),
                                       Accuracy=Accuracy,
                                       Kappa=Kappa,
                                       Sensitivity=Sensitivity,
                                       Specificity=Specificity,
                                       Media=mean(Accuracy, Kappa,
                                                  Sensitivity, Specificity)))
}
print(resultados_caret)





##########################################
###### optimizacion de los modelos


### analisis de importancia


# Separar los datos en variables predictoras (X) y variable objetivo (y)
X <- datos_final[, 1:(ncol(datos_final)-1)] # Variables predictoras (genes)
y <- as.factor(datos_final$Diagnosis) # Variable objetivo

# Algoritmo randomForest con importance=tRUE
rf_model <- randomForest(X, y, ntree = 500, importance = TRUE)

importance <- as.data.frame(importance(rf_model))
head(importance,5)

# Ordenamos los genes por importancia
importance<-importance[order(importance$MeanDecreaseAccuracy, decreasing=TRUE),]
head(importance,10)

# Grafico de importancia
varImpPlot(rf_model, sort=TRUE, n.var=20,main='Gráfico de importancia de las variables')

nrow(importance[importance$MeanDecreaseAccuracy>0,])
# Un total de 15233 genes no restan accuracy al modelo. el resto si. nos quedamos con los nombres:
best_genes <- row.names(importance[importance$MeanDecreaseAccuracy>0,])
head(best_genes,5)


data_train <- data_train[,best_genes]
data_test <- data_test[,best_genes]



#####################################################################
#### se vuelve a ejecutar el código de antes con los modelos de antes.


#### reduccion de dimensionalidad con PCA


# Primero hacemos el PCA con el conjunto de entrenamiento
pca_model <- prcomp(data_train, scale. = TRUE)

# Se obtiene la varianza explicada
varianza_explicada <- pca_model$sdev^2 / sum(pca_model$sdev^2)

# Crear un gráfico de barras para la varianza explicada por las primeras 10 componentes
df_varianza <- data.frame(Componente = 1:length(pca_model$sdev), Varianza_Explicada = varianza_explicada)

ggplot(head(df_varianza, 10), aes(x = factor(Componente), y = Varianza_Explicada)) +
  geom_bar(stat = "identity", fill = "darkturquoise") +
  labs(title = "Varianza explicada por las primeras 10 PC",
       x = "Componente Principal",
       y = "Varianza Explicada") +
  theme_minimal()


# Nos quedamos con las componentes que explican el 95% de la varianza
varianza_acumulada <- cumsum(pca_model$sdev^2) / sum(pca_model$sdev^2)
num_pcs_a_retener <- which(varianza_acumulada >= 0.95)[1]

# Número de PCs a retener: 212 COMPONENTES

# Transformamos los datos de test y train con este modelo
data_train_pca <- predict(pca_model, data_train)[, 1:num_pcs_a_retener]
data_test_pca <- predict(pca_model, data_test)[, 1:num_pcs_a_retener]

# Y ahora se usan los nuevos data train y data test

# Modelo SVM con PCA

modelo_SVM <- ksvm(as.factor(data_train_labels) ~ .,
                   data = data_train_pca, kernel = "vanilladot")

predict_SVM <- predict(modelo_SVM, data_test_pca)

matriz <- confusionMatrix(predict_SVM, as.factor(data_test_labels))
Accuracy <- matriz$overall["Accuracy"]
Kappa <- matriz$overall["Kappa"]
Sensitivity <- mean(matriz$byClass[,1])
Specificity <- mean(matriz$byClass[,2])



resultados_SVM <- rbind(resultados_SVM, data.frame(row.names=i, Accuracy=Accuracy,
                                                   Kappa=Kappa,
                                                   Sensitivity=Sensitivity,
                                                   Specificity=Specificity,
                                                   Media=mean(Accuracy,Kappa,
                                                              Sensitivity,
                                                              Specificity)))

print(resultados_SVM)




### Probar PCA con otros parametros
# Grid para buscar hiperparametros
parametros_svm <- expand.grid(C = c(0.001, 0.01, 0.1, 1, 10))


modelo_SVM <- train(
  Diagnosis ~ .,
  data = data_train_2,
  method = "svmLinear",
  trControl = trainControl(method = "cv", number = 5),
  tuneGrid = parametros_svm
)

predict_SVM <- predict(modelo_SVM, data_test)

matriz <- confusionMatrix(predict_SVM, as.factor(data_test_labels))
Accuracy <- matriz$overall["Accuracy"]
Kappa <- matriz$overall["Kappa"]
Sensitivity <- mean(matriz$byClass[,1])
Specificity <- mean(matriz$byClass[,2])



resultados_SVM <- rbind(resultados_SVM, data.frame(row.names=i, Accuracy=Accuracy,
                                                   Kappa=Kappa,
                                                   Sensitivity=Sensitivity,
                                                   Specificity=Specificity,
                                                   Media=mean(Accuracy,Kappa,
                                                              Sensitivity,
                                                              Specificity)))

print(resultados_SVM)


### No mejoró el rendimiento, es el mismo.
print(matriz)



#### grafico de metricas de evaluacion y guardar los modelos finales

library(ggplot2)
medidas_rendimiento <- data.frame(
  metrica = c("Exactitud", "Valor Kappa","Sensibilidad", "Especificidad",
              "Media (%)"),
  valor = cbind(as.numeric(resultados_SVM[1,]))
)


# Crear el gráfico de barras
ggplot(medidas_rendimiento, aes(x = metrica, y = valor)) +
  geom_bar(stat = "identity", fill = "turquoise") +
  ylim(0, 1) +
  labs(title = "Medidas de Rendimiento del Modelo",
       x = "Medida de Rendimiento",
       y = "Valor") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Se guardaron las variables usadas al final (los genes escogidos), los valores máximos y mínimos de los conjuntos de train y test utilizados con esas variables, y el modelo.

saveRDS(modelo_SVM,file='modelo_SVM_843.Rda')
genes_modelo <- colnames(data_test)
saveRDS(genes_modelo,file='genes_modelo_SVM.Rda')







########################################################
################## INTENTO CORRECCION #################
################## datos balanceados #################



### DATA TEST REAL

### Hay que excluir los datos generados por SMOTE del conjunto TEST.
## Tenia guardados en otra carpeta los datos balanceados usados en el modelo
# Los divido en train y test, como antes
# Y ahora lo que hago es excluir de test las filas cuyo nombre no empiece por
# TCGA, ya que esas son las generadas por SMOTE.


datos_balanceados <- readRDS(file='C:/Users/AlbiusBlack/Documents/NaiveBayesHS/datos_balanceados.Rda')


labels <- datos_balanceados$Diagnosis

# Si ejecutamos table(labels) veremos que tenemos 320 datos de cada uno,
# igual que usamos en el modelo SVM1
datos_balanceados <- log2(datos_balanceados[,-ncol(datos_balanceados)]+1)



seed <- 144
perc <- 0.8
set.seed(seed)

# Con sample sacamos los indices con los que obtenemos los datos de train
train_indices <- sample(nrow(datos_balanceados), floor(perc * nrow(datos_balanceados)))
data_train <- datos_balanceados[train_indices, 1:(ncol(datos_balanceados))]  # Datos para train
data_test <- datos_balanceados[-train_indices, 1:(ncol(datos_balanceados))]  # Datos para test

#Etiquetas
data_test_labels <- labels[-train_indices]
data_train_labels <- labels[train_indices]

# Borro de Test y Test labels las filas que no quiero
data_test_2 <- cbind(data_test, Diagnosis=data_test_labels)
data_test_2 <- data_test_2[rownames(data_test) %like% "TCGA", ]
data_test <- data_test_2[,-ncol(data_test_2)]
data_test_labels <- data_test_2$Diagnosis
table(data_test_labels)
maximo <- max(data_train)
minimo <- min(data_train)


# Función para normalizar
normalize <- function(x) {
  return((x - minimo) / (maximo - minimo))
}



# Aplicamos la función
data_train <- as.data.frame(lapply(data_train, normalize))
data_test <- as.data.frame(lapply(data_test, normalize))


# Y obtengo las metricas del modelo

library(kernlab)
library(caret)


########################### SVM ###########################################



resultados_SVM <- data.frame('Accuracy'=numeric(), 'Kappa' =numeric(),
                             'Sensitivity'=numeric(),'Specificity' = numeric(),
                             'Media'=numeric())




set.seed(seed)
# Usando el modelo que ya tenía guardado. El que está en la aplicación
modelo_SVM <- readRDS('inst/extdata/modelo_SVM_843.Rda')

predict_SVM <- predict(modelo_SVM, data_test)
# Valores de interés
matriz <- confusionMatrix(predict_SVM, as.factor(data_test_labels))
Accuracy <- matriz$overall["Accuracy"]
Kappa <- matriz$overall["Kappa"]
Sensitivity <- mean(matriz$byClass[,1])
Specificity <- mean(matriz$byClass[,2])



resultados_SVM <- rbind(resultados_SVM, data.frame(row.names='test_real', Accuracy=Accuracy,
                                                   Kappa=Kappa,
                                                   Sensitivity=Sensitivity,
                                                   Specificity=Specificity,
                                                   Media=mean(Accuracy,Kappa,
                                                              Sensitivity,
                                                              Specificity)))




# Imprimir la tabla de resultados
print(resultados_SVM)
matriz


evaluacion <- as.matrix(matriz$byClass)
write.csv(evaluacion, file='C:/Yo/MÁSTER/UOC/TFM/rabajo/REALES.csv',
          col.names = TRUE, row.names = TRUE, append = FALSE)









############################################################
####### CODIGO SCRIPT CON MI MODELO

### Script que leerá la app para normalizar y obtener predicciones
modelo_SVM <- function(data, data_labels, modelo,genes_modelo) {



  minimo <- 0
  maximo <- 19.15297


  data <- data[,genes_modelo]


  ### Normalizacion
  data <- log2(data[,1:(ncol(data))]+1)


  # Función para normalizar
  normalize <- function(x) {
    return((x - minimo) / (maximo - minimo))
  }



  # Aplicamos la función
  data <- as.data.frame(lapply(data, normalize))


  if (prop.table(table(is.na(data)))['FALSE']!=1){
    return('Inserta datos sin valores faltantes')
  }else {


    predict <- predict(modelo, data)
    results <- data.frame('Etiqueta'=character(), 'Tumor'=character())

    for (i in 1:length(predict)) {
      results <- rbind(results, data.frame(Etiqueta=data_labels[i,1],
                                           Tumor=as.character(predict[i])))
      #print('Los resultados de las muestras son: ')
      #print(results)

    }
    return(results)
  }
}



#############################################################
 ################# Aplicacion de Shiny##################

library(shiny)





# R script:
source("modelos/SVM_RCC.R")



options(shiny.maxRequestSize = 60*1024^2)




ui <- fluidPage(



  titlePanel("Clasificación de Tipos de Tumor"),
  sidebarLayout(
    sidebarPanel(


      fileInput("dataFile", "Cargar archivo CSV de datos de expresión de genes"),
      fileInput("labelFile", "Cargar archivo de etiquetas"),
      actionButton("predictButton", "Ejecutar Predicción"),

    ),





    mainPanel(
      tabsetPanel(
        tabPanel("Resultado",
                 fluidRow(
                   column(6,tableOutput("prediction"))
                 )),
        tabPanel("How to Use",
                 fluidRow(
                  column(7,
                          h4("Ejemplo de CSV de datos:"),imageOutput("dataImage")),
                   column(5, h4("Ejemplo de CSV de etiquetas:"),imageOutput("labelsImage")),

                   column(7,
                          downloadLink("downloadData", "Descargar ejemplo de datos")),

                   column(5,

                          downloadLink("downloadLabels", "Descargar ejemplo de etiquetas"))
                 )
        ),
        tabPanel("Rendimiento",
               imageOutput("appImage")
        )


      )
    )
  )
)



server <- function(input, output, session) {

  # Leer datos
  data <- reactive({
    inFile <- input$dataFile
    if (is.null(inFile))
      return(NULL)
    read.csv(inFile$datapath)
  })

  # Leer etiquetas
  labels <- reactive({
    inFile <- input$labelFile
    if (is.null(inFile))
      return(NULL)
    read.csv(inFile$datapath)
  })






  output$appImage <- renderImage({
    list(
      src = "www/medidas_modelo_SVM.png",
      width = "70%",
      height = "auto",
      class = "centered-image",
      alt = "Imagen que muestra el rendimiento del modelo. Accuracy
      0.66 Kappa value 0.55 Sensitivity 0.57 Specificity 0.92 Precision 0.54
      Recall 0.57"

    )
  },deleteFile = FALSE)

  output$dataImage <- renderImage({
    list(
      src = "www/data_ejemplo.png",
      width = "100%",
      height = "auto",
      class = "centered-image",
      alt = "Imagen que muestra un ejemplo de cómo debe ser el CSV de datos"

    )
  },deleteFile = FALSE)

  output$labelsImage <- renderImage({
    list(
      src = "www/data_labels_ejemplo.png",
      width = "60%",
      height = "auto",
      class = "centered-image",
      alt = "Imagen que muestra un ejemplo de cómo debe ser el CSV de etiquetas"

    )
  },deleteFile = FALSE)

  ###### Descargar ejemplos de datos y etiquetas
  output$downloadData <- downloadHandler(
    filename = function() {
      "data_test2.csv"
    },
    content = function(file) {
      file.copy("www/data_test2.csv", file)
    }
  )

  # Crea el ejemplo de CSV de etiquetas
  output$downloadLabels <- downloadHandler(
    filename = function() {
      "labels2.csv"
    },
    content = function(file) {
      file.copy("www/labels2.csv", file)
    }
  )


  observeEvent(input$predictButton, {
    # Mensaje de depuración para verificar si se está ejecutando el evento
    print("Ejecutando...")

    # Mensaje de error
    tryCatch({
      # Cargar el modelo
      genes_modelo <- readRDS("modelos/genes_modelo_SVM.Rda")
      modelo <- readRDS("modelos/modelo_SVM_843.Rda")

      # Mensaje de depuración para verificar si se cargó el modelo
      print("Modelo cargado")

      output$prediction <- renderTable({
        # Realizar predicciones con el algoritmo
        prediction <- modelo_SVM(data(), labels(), modelo, genes_modelo)

        # Mensaje de depuración para verificar si se realizó la predicción
        print("Prediccion hecha")

        # Devolver el resultado de la predicción
        return(prediction)
      }, striped = TRUE)

    }, error = function(e) {
      # Manejar el error y mostrar un mensaje
      output$prediction <- renderTable({
        return(data.frame(Error = c("Error al cargar el modelo. Verifica los archivos de modelo.")))
      }, striped = TRUE)
      # Mensaje de depuración para imprimir el mensaje de error
      print(paste("Error: ", e))
    })
  })
}

shinyApp(ui, server)




