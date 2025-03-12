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

logo_files_funders <- c(
  "logo-nidis.png",
  "logo-NRCS.png",
  "logo-ewcsig-sqr_10inprnt.jpg",
  "logo-PI-CASC-marble-w-full-name.jpg",
  "logo-cwrm.png"
)

logo_files_partners <- c(
  "logo-USFS.png",
  "logo-UH_Manoa.png",
  "logo-SeaGrant.jpg",
  "logo-clark-university.png",
  "logo-mesonet-Ringshaped.jpg",
  "logo-PDKE_Logo_Color_Rounded_Type-03.jpg",
  "logo-HCDP_Text_Bottom_Color.jpg"
)

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
      .inline-container {
        display: flex;
        justify-content: center;
        align-items: center;
      }
      .inline-container img {
        margin-right: 15px;
      }
      .inline-container h1 {
        margin: 0;
      }
      .logo-container {
        display: flex;
        flex-wrap: wrap;
        justify-content: center; /* Center horizontally */
      }
      .logo-item {
        margin: 10px; /* Adjust spacing between logos */
        max-width: 150px; /* Adjust maximum logo width */
        max-height: 100px; /* Adjust maximum logo height */
        flex: 0 0 auto; /* Prevent logos from stretching */
      }
      .logo-item img {
        max-width: 100%;
        max-height: 100%;
        height: auto;
        width: auto;
        display: block; /* Remove extra space below images */
      }
    "))
  ),

  # Display an image above the main content
  #img(src = "logo-PDKE_Logo_Color_Rounded_Type-03.jpg", height = "200px", alt = "PDKE Logo", align = "center"),
  
  # Application title
  #titlePanel("Hawaiʻi Climate Portfolios"),

  # Container for the inline image and title
  div(
    class = "inline-container",
    img(src = "logo-PDKE_Logo_Color_Rounded_Type-03.jpg", height = "60px", alt = "PDKE Logo"),
    h1("Hawaiʻi Climate Portfolios")
  ),
  
  # Sidebar layout with input and output definitions
  sidebarLayout(
    # Sidebar panel for inputs
    sidebarPanel(

      # Dropdown menu to select a shapefile
      uiOutput("shapefile_select"),
      
      div(
        style = "display: flex; align-items: center;",
        textInput("polygon_name", "Location name:", ""),
        div(id = "info_polygon_name", class = "info-circle", "?")
      ),
      # Add text input for polygon short name with an info button
      div(
        style = "display: flex; align-items: center;",
        textInput("polygon_short_name", "Location short name:", ""),
        div(id = "info_polygon_short_name", class = "info-circle", "?")
      ),
      # Text input for email address
      textInput("email", "Your email address:", value = "", placeholder = "your.email@example.com"),
      # Warning message for email validation
      tags$div(id = "email_warning", style = "color: red; margin-top: 12px;",
        tags$p("An email address is required since results will be emailed to you once generated.  The submit button will not be enabled until a shape is selected or drawn and all values are provided.")
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
      actionButton("save_button", "Use Selected Area"),
      
      # Placeholder for submission confirmation message
      textOutput("status_messages")
    )
  ),
  
  # Header for the first row of images
  h3(""),
  h3("Funders:"),
  div(class = "logo-container",
    lapply(logo_files_funders, function(file) {
      div(class = "logo-item", img(src = file))
    })
  ),
  h3("Partners:"),
  div(class = "logo-container",
    lapply(logo_files_partners, function(file) {
      div(class = "logo-item", img(src = file))
    })
  )
)
