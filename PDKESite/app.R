library(shiny)

# Source the ui and server scripts
source("ui.R")
source("server.R")

# Run the application
shinyApp(ui = ui, server = server)