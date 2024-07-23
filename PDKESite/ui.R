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
  tags$head(
    tags$style(HTML("
      .info-circle {
        background-color: #007bff;
        color: white;
        border-radius: 50%;
        width: 20px;
        height: 20px;
        display: inline-block;
        text-align: center;
        line-height: 20px;
        margin-left: 5px;
        cursor: pointer;
      }
    "))
  ),
  
  # Application title
  titlePanel("Leaflet Map with Shapefile Selector"),
  
  # Sidebar layout with input and output definitions
  sidebarLayout(
    # Sidebar panel for inputs
    sidebarPanel(
      # Text input for email address
      textInput("email", "Enter your email address:", value = "", placeholder = "your.email@example.com"),
      # Warning message for email validation
      tags$div(id = "email_warning", style = "color: red; margin-top: 12px;",
        tags$p("An email address is required since results will be emailed to you once generated.  The submit button will not be enabled until a shape is selected or drawn and all values are provided.")
      ),
      
      # Dropdown menu to select a shapefile
      uiOutput("shapefile_select"),
      #textInput("polygon_name", "Polygon name:", ""),
      #textInput("polygon_short_name", "Polygon short name:", ""),
      
      div(
        style = "display: flex; align-items: center;",
        textInput("polygon_name", "Polygon name:", ""),
        div(id = "info_polygon_name", class = "info-circle", "?")
      ),
      # Add text input for polygon short name with an info button
      div(
        style = "display: flex; align-items: center;",
        textInput("polygon_short_name", "Polygon short name:", ""),
        div(id = "info_polygon_short_name", class = "info-circle", "?")
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

    Shiny.addCustomMessageHandler("bindInfoCircleClick", function(message) {
      document.getElementById(message.id).addEventListener("click", function() {
        Shiny.setInputValue(message.event, Math.random());
      });
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
