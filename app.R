library(shiny)


# Instalar y cargar el paquete desde GitHub

if (!requireNamespace("svmRCC", quietly = TRUE)) {
  # Si no está instalado, instalarlo
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  remotes::install_github("hugosuarezglez/svmRCC")
}


# Cargar el paquete
library(svmRCC)



options(shiny.maxRequestSize = 60*1024^2)




ui <- fluidPage(
  #shinyjs::useShinyjs(), #Para mensaje de estado


  titlePanel("Clasificación de Tipos de Tumor"),
  sidebarLayout(
    sidebarPanel(

      #tabsetPanel(
                 fileInput("dataFile", "Cargar archivo CSV de datos de expresión de genes"),
                 fileInput("labelFile", "Cargar archivo de etiquetas"),
                 actionButton("predictButton", "Ejecutar Predicción"),
                 # Imagen en la pagina principal?
                 #br(),
                 #imageOutput("appImage")
    ),
# cortar desde tabPanel("how to use") hasta antes de mainpanel




        mainPanel(
          tabsetPanel(
            tabPanel("Resultado",
                     fluidRow(
                       column(6,tableOutput("prediction"))
                     )),
             tabPanel("How to Use",
                      fluidRow(
                        #column(12,
                        #      tags$img(src = "C:/Yo/MÁSTER/UOC/TFM/imgs/data_ejemplo.png", width = "100%")),
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
                      # fluidRow(
                      #   column(12,
                      #          imageOutput("appImage", width = "70%"))
                      #
                      imageOutput("appImage")
             )
             #)

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


    output$prediction <- renderTable({
      # Realizar predicciones con el algoritmo Naive Bayes
      prediction <- svmRCC::modelo_SVM(data(), labels(),
                                       ruta_modelo = "modelos/modelo_SVM_843.Rda",
                                       ruta_genes = "modelos/genes_modelo_SVM.Rda")


      # Devolver el resultado de la predicción
      return(prediction)
    }, striped = TRUE)
  })
}

shinyApp(ui, server)

