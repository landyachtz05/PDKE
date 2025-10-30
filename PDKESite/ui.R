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
  "logo-PI-CASC-marble-w-full name.jpg",
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
      },
      .btn-primary {
        background-color: #007bff; /* Blue background color */
        border-color: #007bff; /* Blue border color */
        color: white; /* White text color */
      }
      .btn-primary:hover {
        background-color: #0056b3; /* Darker blue on hover */
        border-color: #0056b3;
      }
      .main-header {
        display: flex;
        flex-direction: column; /* Stack centered row and left-aligned contact */
        align-items: center; /* Center the first row horizontally */
      }
      .centered-row {
        display: flex;
        align-items: center; /* Vertically center the logo and title */
      }
      .contact-line {
        align-self: flex-end; /* Align the contact line to the left */
        margin-top: 5px; /* Add some spacing */
      }
      .explanation-text {
        margin-top: 15px;
        font-size: 16px;
        color: #555;
        line-height: 1.5;
      }
    "))
  ),

  # Display an image above the main content
  #img(src = "logo-PDKE_Logo_Color_Rounded_Type-03.jpg", height = "200px", alt = "PDKE Logo", align = "center"),
  
  # Application title
  #titlePanel("Hawaiʻi Climate Portfolios"),

  # Container for the inline image and title
  div(
    class = "main-header",
    div(
      class = "centered-row",
      img(src = "logo-PDKE_Logo_Color_Rounded_Type-03.jpg", height = "60px", alt = "PDKE Logo"),
      h1("Hawaiʻi Climate Portfolios")
    ),
    div(
      class = "contact-line",
      tags$p(
        "Contact Us: ",
        tags$a(href = "mailto:djford@hawaii.edu", "djford@hawaii.edu")
      )
    ),
    # Add this new div for your explanatory text
    div(
      class = "explanation-text",
      style = "text-align: left; margin: 15px 0 25px 0; max-width: 800px;",
      tags$p(
       HTML("This <strong>CCVD portfolio</strong> is a comprehensive synthesis of climate and drought information developed for any area of interest within the State of Hawaii. It is designed to provide relevant climate and drought information for a wide variety of users at the scale of one-hundred acres to whole islands. When generated, each CCVD portfolio is delivered via email as a Microsoft PowerPoint presentation of approximately 40 slides and 16MB in size. To learn more or to view sample portfolios, please visit the <a href='https://www.soest.hawaii.edu/pdke/'>PDKE website</a>.<br><br>
    <strong>Instructions for use:</strong><br>
    <strong>1)</strong> Select your area of interest (two options):<br>
    <span style='margin-left: 20px;'>a) Select a predefined location by choosing from the dropdown menu then selecting on the map.<br>
    <span style='margin-left: 20px;'>b) Create a custom area of interest by using the area drawing tools found at the top-left of the map.<br>
    <strong>2)</strong> Enter the Location name and Location short name. If you selected a pre-defined location, this will be auto-populated but can be edited if needed. Do not use any diacriticals or special characters as this will result in a portfolio download error.<br>
    <strong>3)</strong> Enter the email address where the portfolio will be sent once generated.<br>
    <strong>4)</strong> Click the “Generate  Report” button.<br><br>
    <strong>Disclaimer</strong>: Some data used within this portfolio are produced in near-real-time, and may be subject to futher quality control measures. For any questions or concerns regarding the data or content of a CCVD portfolio, please email the contact listed in the top right corner of this webpage.")
      )
    )
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
      tags$div(id = "email_warning", style = "color: grey; margin-top: 12px;",
        tags$p("An email address is required since results will be emailed to you once generated.  The submit button will not be enabled until a shape is selected or drawn and all values are provided.")
      ),
      div(
        # Add a button to save the selected area as a shapefile
        actionButton("save_button", "Generate Report", class = "btn-primary"),
        style = "text-align: center;" # Inline CSS for centering
      ),
      # Placeholder for submission confirmation message
      # textOutput("status_messages"),
      uiOutput("status_messages"),
      
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

        $(document).on("shiny:inputchanged", function(event) {
          if (event.name === "polygon_name" && event.value.length > 150) {
            $("#polygon_name").val(event.value.substring(0, 150));
          }
          if (event.name === "polygon_short_name" && event.value.length > 100) {
            $("#polygon_short_name").val(event.value.substring(0, 100));
          }
          if (event.name === "email" && event.value.length > 150) {
            $("#email").val(event.value.substring(0, 150));
          }
        });
    '))
    ),
    
    # Main panel for displaying outputs
    mainPanel(
      # Leaflet map
      leafletOutput("map", height = 600)
    )
  ),
  
  # Header for the first row of images
  h3(""),
  h3("Funders:"),
  div(
    class = "logo-container",
    lapply(seq_along(logo_files_funders), function(i) {
      urls_funders <- c(
        "https://www.drought.gov/",
        "https://www.nrcs.usda.gov/",
        "https://www.eastwestcenter.org/",
        "https://pi-casc.soest.hawaii.edu/",
        "https://dlnr.hawaii.gov/cwrm/"
      )
      div(
        class = "logo-item",
        tags$a(
          href = urls_funders[i],
          target = "_blank",
          img(src = logo_files_funders[i], alt = basename(logo_files_funders[i]))
        )
      )
    })
  ),
  
  h3("Partners:"),
  div(
    class = "logo-container",
    lapply(seq_along(logo_files_partners), function(i) {
      urls_partners <- c(
        "https://www.fs.usda.gov/",
        "https://manoa.hawaii.edu/",
        "https://seagrant.soest.hawaii.edu/",
        "https://www.clarku.edu/",
        "https://www.hawaii.edu/climate-data-portal/hawaii-mesonet/",
        "https://www.soest.hawaii.edu/pdke/",
        "https://www.hawaii.edu/climate-data-portal/"
      )
      div(
        class = "logo-item",
        tags$a(
          href = urls_partners[i],
          target = "_blank",
          img(src = logo_files_partners[i], alt = basename(logo_files_partners[i]))
        )
      )
    })
  )
)