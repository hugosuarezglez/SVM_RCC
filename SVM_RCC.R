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
