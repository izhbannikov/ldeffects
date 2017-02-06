# Utility functions #

outputDir <- "." #"saved"

getState<-function(input) c(input$slider1, input$slider2, input$slider3, input$slider4)


saveData <- function(data) {
  data <- t(data)
  # Create a unique file name
  fileName <- sprintf("%s_%s.csv", as.integer(Sys.time()), digest::digest(data))
  # Write the file to the local system
  write.csv(
    x = data,
    file = file.path(outputDir, fileName), 
    row.names = FALSE, quote = TRUE
  )
}

loadData <- function() {
  # Read all the files into a list
  files <- list.files(outputDir, full.names = TRUE,pattern = "*.csv")
  #print(files)
  data <- lapply(files, read.csv, stringsAsFactors = FALSE,row.names=NULL) 
  #print(data)
  # Concatenate all data together into one data.frame
  data <- do.call(rbind, data)
  data
}
