modelo_SVM <- function(data, data_labels, ruta_modelo, ruta_genes) {
  modelo <- readRDS(file = ruta_modelo)
  genes_modelo <- readRDS(file = ruta_genes)


  minimo <- 0
  maximo <- 19.15297
#### ARCHIVOS
#
#
# # URL del archivo en GitHub
# url_modelo <- "https://github.com/hugosuarezglez/svmRCC/raw/main/inst/extdata/modelo_SVM_843.Rda"
# url_maximo <- "https://github.com/hugosuarezglez/svmRCC/raw/main/inst/extdata/maximo_train.Rda"
# url_minimo <- "https://github.com/hugosuarezglez/svmRCC/raw/main/inst/extdata/minimo_train.Rda"
# url_genes <- "https://github.com/hugosuarezglez/svmRCC/raw/main/inst/extdata/genes_modelo_SVM.Rda"
#
# # Descargar el archivo RDA
# maximo <- httr::GET(url_maximo, httr::write_disk("maximo_train.Rda", overwrite = TRUE))
# minimo <- httr::GET(url_minimo, httr::write_disk("minimo_train.Rda", overwrite = TRUE))
# modelo <- httr::GET(url_modelo, httr::write_disk("modelo_SVM_843.Rda", overwrite = TRUE))
# genes_modelo <- httr::GET(url_genes, httr::write_disk("genes_modelo_SVM.Rda", overwrite = TRUE))
# # Cargar el contenido del archivo RDA en R
# maximo <- readRDS("maximo_train.Rda")
# minimo <- readRDS("minimo_train.Rda")
# modelo <- readRDS("modelo_SVM_843.Rda")
# genes_modelo <- readRDS("genes_modelo_SVM.Rda")



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
