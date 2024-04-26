#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(leaflet)

# # Define the user interface
# ui <- fluidPage(
#   titlePanel("Leaflet Map in Shiny with Drawing Tools"),
#   
#   # Leaflet map with drawing tools
#   leafletOutput("map", height = 600),
#   
#   # Add a button to save the selected area as a shapefile
#   actionButton("save_button", "Save Selected Area as Shapefile")
# )

ui <- fluidPage(
  # Application title
  titlePanel("Leaflet Map with Shapefile Selector"),
  
  # Sidebar layout with input and output definitions
  sidebarLayout(
    # Sidebar panel for inputs
    sidebarPanel(
      # Dropdown menu to select a shapefile
      uiOutput("shapefile_select")
    ),
    
    # Main panel for displaying outputs
    mainPanel(
      # Leaflet map
      leafletOutput("map", height = 600),
      
      # Add a button to save the selected area as a shapefile
      actionButton("save_button", "Save Selected Area as Shapefile")
    )
  )
)