library(shiny)
library(ggplot2)
library(dplyr)

# upload and merge 2 datasets with cultures and subcultures

#cultures <- read.table("/home/rostyslav/Desktop/R_analysis/PhD_project/Quiescence_mutation_detection/data/output/4cultures_preprocessed.txt", header = T)
#subcultures <- read.table("/home/rostyslav/Desktop/R_analysis/PhD_project/Quiescence_mutation_detection/data/output/6cultures_preprocessed.txt", header = T)
#final <- rbind(cultures, subcultures)


#final <- read.table("/home/rostyslav/Desktop/R_analysis/PhD_project/NGS_sequencing/final_variants_and_table_total.txt", header = T)
#final < df

ui <- fluidPage(
  titlePanel("Frequency distribution of the variants"),
  sidebarPanel(selectInput("sampleInput", "Sample name",
                           sort(unique(final$Sample)),
                           selected = "culture_0"),
               radioButtons("timeInput", "Time in quiescence",
                            choices = c("2_months", "3_months"),
                            selected = "2_months"),
               radioButtons("mutation_typeInput", "Mutation type",
                            choices = c("SNP", "Indel"),
                            selected = "SNP"),
               sliderInput("CoverageInput", "Coverage", 0, 150000, c(0, 150000)), # Coverage: reads/position, need to update
#               numericInput("VarFreqInput", "Variant Frequency", 0, 1, 0.001),
               sliderInput("VarFreqInput", "Variant Frequency", min=0, max=0.1,  value=c(0, 0.1), step=0.001), # 
               sliderInput("Qual2Input", "Base Quality", 0, 39, c(0, 39)), # Quality of an alternative variants, need to update
               sliderInput("ReadProportionInput", "Distribution within the PairEnd Reads", min=0, max=100, value=c(0, 100), step=0.01), # need to update
               sliderInput("PvalueInput", "p-value (Fisher's test)", 0, 1, c(0, 1))),
  mainPanel(plotOutput("coolplot"),
            br(), br(),
            textOutput("selected_var"))
  )

server <- function(input, output) {
  output$coolplot <- renderPlot({
    filtered <- reactive({
      final %>%
      filter(Sample == input$sampleInput,
             time_in_quiescence == input$timeInput,
             mutation_type == input$mutation_typeInput,
             Coverage >= input$CoverageInput[1],
             Coverage <= input$CoverageInput[2],
             VarFreq >= input$VarFreqInput[1],
             VarFreq <= input$VarFreqInput[2],
             Qual2 >= input$Qual2Input[1],
             Qual2 <= input$Qual2Input[2],
             ReadProportion >= input$ReadProportionInput[1],
             ReadProportion <= input$ReadProportionInput[2],
             Pvalue >= input$PvalueInput[1],
             Pvalue <= input$PvalueInput[2]
      )
    })
      
    ggplot(data = filtered(), aes(VarFreq))+
      geom_histogram(binwidth = 0.001, show.legend = T)+
      coord_cartesian(xlim = c(0, 0.1), ylim=c(0,130))+
#      scale_x_continuous(labels = VarFreq, breaks = seq(0, 0.1, 0.01), 
#                         minor_breaks = seq(0, 0.1, 0.001), 
#                         name = "Variant frequency")+
      scale_y_continuous(breaks = seq(0, 200, 10), 
                         minor_breaks = seq(0, 200, 1), 
                         name = "Number of variants")+
      theme(axis.title=element_text(family = "Times New Roman", 
                                    face="bold", 
                                    size=14,
                                    color="black"), 
            axis.text=element_text(family = "Times New Roman", 
                                   face="bold", 
                                   size=14,
                                   color="black"))
  })
  
#  output$selected_var <- renderText({ 
#    paste("The final number of variants (mutations) is", input$nrow(final)
#)
#  })
}

shinyApp(ui = ui, server = server)
