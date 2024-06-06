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

ui <- fluidPage(
  # Application title
  titlePanel("Leaflet Map with Shapefile Selector"),
  
  # Sidebar layout with input and output definitions
  sidebarLayout(
    # Sidebar panel for inputs
    sidebarPanel(
      # Dropdown menu to select a shapefile
      uiOutput("shapefile_select"),
      # Text input for email address
      textInput("email", "Enter your email address:", value = "", placeholder = "your.email@example.com"),
      
      # Warning message for email validation
      tags$div(id = "email_warning", style = "color: red; margin-top: 12px;",
        tags$p("An email address is required to submit.")
      ),
      
      # JavaScript code snippet to handle email validation
      tags$script(HTML('
    $(document).on("shiny:inputchanged", function(event) {
      if (event.name === "email_input") {
        var email = event.value;
        if (email.length === 0) {
          $("#email_warning").show();
          $("#save_button").prop("disabled", true);
        } else {
          $("#email_warning").hide();
          $("#save_button").prop("disabled", false);
        }
      }
    });
  '))
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
