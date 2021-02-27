appDir <- getwd()
monte.carlo.default <- 1000
MonteCarlo <- getShinyOption("MonteCarlo", monte.carlo.default)

shinyApp(
  ui = fluidPage(paste("Chosen parameter:", MonteCarlo)),
  server = function(input, output, session){
    oldwd <- setwd(appDir)
    on.exit(setwd(oldwd))
    ## put more server logic here
  }
)
