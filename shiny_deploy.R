# A simple script to deploy on shinyapps:
install.packages("devtools")
library(devtools)
devtools::install_github("rstudio/shinyapps")
library(rsconnect)
rsconnect::setAccountInfo(name='ilmadester',token='1D5950B47F9FFD38E2ED2609985BAF8A',secret='fbL28KfBVVOGAG2NNR6kFpUW6gIcKlZZ04Fxpmnx')
rsconnect::deployApp("~/Projects/ld_effects/", account="ilmadester")
