library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Exploring effects of LD"),
  downloadLink("downloadDocs", "Download documentation"),
  fluidRow(
    column(4,
           wellPanel(
             
             h4("Time"),
             sliderInput("time",
                  "t:",
                  min = 1,
                  max = 120,
                  value = c(30,105))),
           wellPanel(
             h4("Proportion hazard"),
             checkboxInput("constrained", "Freeze P1,P2", value = FALSE),
             conditionalPanel("input.constrained == 0",
              sliderInput("slider3","m10:", min = 0,max = 1,value = 0.1),
              sliderInput("slider1", "m00:", min = 0, max = 1, value = 0.48)),
             sliderInput("slider4", "m11:", min = 0, max = 1, value = 0.34),
             conditionalPanel("input.constrained == 0",
              sliderInput("slider2","m01:", min = 0, max = 1, value = 0.08))
             ),
    wellPanel(
      h4("Mortality"),
      sliderInput("mu00","mu00:",min = 0,max = 0.1,value = 0.013),
      checkboxInput("dcase", "D1,D2", value = FALSE),
      conditionalPanel("input.dcase == 0",
        sliderInput("H1","H1:",min = 0,max = 5,value = 0.5, step=0.05),
        sliderInput("H2","H2:",min = 0,max = 5,value = 1, step=0.05)),
      conditionalPanel("input.dcase == 1",
        sliderInput("D1","D1:",min = -1,max = 1,value = -0.1, step=0.05),
        sliderInput("D2","D2:",min = -1,max = 1,value = 0.2, step=0.05))
    )),
    
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(
        tabPanel("Main Plots", 
                 plotOutput("distPlot",height = 1024), 
                 downloadButton("downloadPlot", label="Save"), 
                 checkboxInput("notitle_main", "No title", value = FALSE)),
        tabPanel("Mortality", 
                 plotOutput("mortalityPlot",height = 480), 
                 downloadButton("downloadPlotMu", label="Save"),
                 checkboxInput("notitle_mortality", "No title", value = FALSE)),
        tabPanel("MAF", 
                 plotOutput("mafPlot",height = 480), 
                 downloadButton("downloadPlotMAF", label="Save"),
                 checkboxInput("notitle_maf", "No title", value = FALSE)),
        tabPanel("LD", 
                 plotOutput("ldPlot",height = 480), 
                 downloadButton("downloadPlotLD", label="Save"),
                 checkboxInput("notitle_ld", "No title", value = FALSE))
      )
      
    )
  )
))