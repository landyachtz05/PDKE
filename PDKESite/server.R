#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#
# using https://maps.equatorstudios.com to test if saved shapefiles are legit
# watershed shapefile from https://geoportal.hawaii.gov/datasets/HiStateGIS::watersheds/explore 
#
# --------------------------------------------------------------------------
# 
# Before running server.R, make sure you have a file called credentials.json
# in the base level dir of the project dir and edit the values to match your 
# environment.  The file should look like this: 
# {
#  
#  "PROJ_LIB_VAL": "/opt/anaconda3/share/proj/",
#  "RSCRIPT_PATH": "/Library/Frameworks/R.framework/Resources/bin/Rscript"
# }
# Set ENV_TYPE to either "linux" or "windows", depending on the box you will be running this on.
# You also need to make an empty file called .here in your project's 
# base level directory.
#

#install.packages("jsonlite")
#install.packages("processx") # Install if you haven't already

library(shiny)
library(leaflet)
library(leaflet.extras)
library(sf)
library(jsonlite)
library(here)
#library(processx)

read_credentials <- function(filepath) {
  tryCatch({
    credentials <- fromJSON(filepath)
    return(credentials)
  }, error = function(e) {
    print(paste("Error reading credentials file:", e$message))
    return(NULL) # Or handle the error as you see fit
  })
}

# default values for prod
ENV_TYPE = "windows"
RSCRIPT_PATH = "/bin/Rscript" 
#RSCRIPT_PATH = paste0(Sys.getenv("R_HOME"), 'Rscript')

#BASE_DIR = "/srv/shiny-server/"
BASE_DIR <- paste0(here(), "/") # Gets the project root which is whereever the .here file is located

status_message = "Your request has been submitted.  In about 45 minutes you will receive an email with a link to the download.  The link expires in 7 days."


# Load credentials
credentials_file <- paste0(BASE_DIR, "/credentials.json")
creds <- read_credentials(credentials_file)
if (!is.null(creds)) {
  ENV_TYPE <- creds$ENV_TYPE
  RSCRIPT_PATH = creds$RSCRIPT_PATH
  #BASE_DIR <- creds$BASE_DIR
} else {
  print("Credentials could not be loaded")
}
print(paste("BASE_DIR:", BASE_DIR))
print(paste("RSCRIPT_PATH:", RSCRIPT_PATH))

location_file = paste0(BASE_DIR, "locations.csv")
myscript_path = paste0(BASE_DIR, "CCVD_portfolio_content.R")
run_string = paste(RSCRIPT_PATH, myscript_path)

validate_text_inputs <- function(email, polygon_name, polygon_short_name) {
  print(paste0("validate_text_inputs: ", email, ", ", polygon_name, ", ", polygon_short_name))
  if (grepl("^[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\\.[A-Za-z]{2,}$", email)) {
    print("validate_text_inputs: 2")
    # make sure none are empty
    if ((nchar(email) > 0) && (nchar(polygon_name) > 0) && (nchar(polygon_short_name))) {
      return(TRUE)
    }
  }
  print("validate_text_inputs: 5, returning false")
  return(FALSE)
}

# write "content" data out to a file, this is intended to keep a record of areas 
# getting used to so we can see if it's valuable to have a few pre-created results.
record_location_to_file <- function(content) {
  # Open the file in append mode
  con <- file(location_file, open = "a")
  
  # Write the content to the file
  writeLines(content, con)
  
  # Close the connection
  close(con)
}

# Function to update the state of the submit button
update_submit_button <- function(session, valid_text_inputs, polygon_selected) {
  print(paste0("update_submit_button: ", valid_text_inputs, ", ", polygon_selected))
  updateActionButton(session, "save_button", label = "Generate Report", disabled = !(valid_text_inputs && polygon_selected))
}

display_intersected_islands <- function(intersected_island_names) {
  island_names_string <- ""
  # Display intersected island names
  if (length(intersected_island_names) > 0) {
    island_names_string <- paste(intersected_island_names, collapse = ", ")
    cat("\nisland_names_string: ", island_names_string, "\n")
    #island_names_text <- paste("The selected area intersects with the following island(s):", paste(intersected_island_names, collapse = ", "))
    #showNotification(island_names_text, type = "message")
  } else {
    showNotification("The selected area does not intersect with any islands.", type = "warning")
  }
}

get_intersected_islands <- function(sf_object, island_boundaries) {
  # Check for intersection with island boundaries
  intersected_islands <- sf::st_intersects(sf_object, island_boundaries)
  
  # Extract island names or identifiers
  intersected_island_names <- unique(island_boundaries$isle[intersected_islands[[1]]])
  display_intersected_islands(intersected_island_names)
  return(intersected_island_names)
}

run_ccvd <- function(session, sf_object, island_boundaries, shapefile_full_path, polygon_name, polygon_short_name, email) {
  ISLAND_FULL_NAMES <- c("Hawaii", "Maui", "Kahoolawe", "Lanai", "Molokai", "Oahu", "Kauai")
  ISLAND_SHORT_NAMES <- c("BI", "MN", "KO", "LA","MO","OA","KA")
  
  # get the short name for the given island name
  island_full_name <- get_intersected_islands(sf_object, island_boundaries)
  index <- match(island_full_name, ISLAND_FULL_NAMES)
  island_short_name <- ""
  # Check if the island was found
  if (!is.na(index)) {
    # Find the associated short name using the index
    island_short_name <- ISLAND_SHORT_NAMES[index]
  }
  cat("shapefile_full_path: ", shapefile_full_path, "\n")
  
  if (ENV_TYPE == "linux") {
    # start of code for mac/linux box
    full_run_string <- paste0(run_string, " ",
      shQuote(email), " ",
      shQuote(paste0(BASE_DIR, "PDKESite/", shapefile_full_path)), " ",
      shQuote(polygon_name), " ",
      shQuote(polygon_short_name), " ",
      shQuote(island_full_name), " ",
      shQuote(island_short_name))
    print(paste0("full_run_string: ", full_run_string))
    # /Library/Frameworks/R.framework/Resources/Rscript /Users/jgeis/Work/PDKE/CCVD_portfolio_content.R 'jgeis@hawaii.edu' '/Users/jgeis/Work/PDKE/PDKESite/Shapefiles/SelectedPolygon/Kaa_2024_07_25_08_21_36.shp' 'Kaa' 'Kaa' 'Lanai' 'LA'"
  
    # this works, gets temporarily commented out so I can test w/o invoking the other stuff
    system(full_run_string, wait = FALSE)
    
    # print(paste0("RSCRIPT_PATH: ", RSCRIPT_PATH))
    # print(paste0("myscript_path: ", myscript_path))
    # print(paste0("BASE_DIR: ", BASE_DIR))
    # print(paste0("shapefile_full_path: ", shapefile_full_path))
    # print(paste0("polygon_name: ", polygon_name))
    # print(paste0("polygon_short_name: ", polygon_short_name))
    # print(paste0("island_full_name: ", island_full_name))
    # print(paste0("island_short_name: ", island_short_name))
  }
  else if (ENV_TYPE == "windows") {
    # not actually used, but useful to print out in case we need to re-run via command line
    full_run_string <-
      paste0(
        shQuote(run_string), " ",
        shQuote(email), " ",
        shQuote(paste0(BASE_DIR, "PDKESite/", shapefile_full_path)), " ",
        shQuote(polygon_name), " ",
        shQuote(polygon_short_name), " ",
        shQuote(island_full_name), " ",
        shQuote(island_short_name))
    print(paste0("full_run_string: ", full_run_string))

    args <- c(
      myscript_path,
      email,
      file.path(BASE_DIR, "PDKESite", shapefile_full_path),
      polygon_name,
      polygon_short_name,
      island_full_name,
      island_short_name
    )

    print("run ccvd")
    result <- processx::run(
      command = RSCRIPT_PATH,  # The Rscript or other command
      args = args,            # Arguments as a character vector
      echo = TRUE              # Print the output to the console
    )

    cat(result$stdout)
    if (result$stderr != "") {
      cat("Error:", result$stderr)
    }
  }

  # Disable the save button to prevent a re-submission.
  updateActionButton(session, "save_button", disabled = TRUE)
  
  datetime_str <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  csv_output_string <- paste0(datetime_str, ",", shapefile_full_path, ", ", island_full_name, ", ", polygon_name)
  print(paste0("csv_output_string: ", csv_output_string))
  record_location_to_file(csv_output_string)
}


# Define server logic
server <- function(input, output, session) {
  # Initialize a reactive value to store the drawn feature
  drawnFeature <- reactiveVal(NULL)
  
  # Load Hawaiian island boundaries
  # Replace 'Coastline.shp' with the path to your island boundaries shapefile
  # data from: https://geoportal.hawaii.gov/datasets/045b1d5147634e2380566668e04094c6/explore  
  island_boundaries <- sf::st_read("Coastline/Coastline.shp")
  
  ########## stuff for drop-down menu - start   ##########  
  # Get a list of shapefiles in the "shapefiles" directory
  shapefile_paths <- list.files("Shapefiles", pattern = "\\.shp$", full.names = TRUE)
  
  # Create a reactive value to store the selected shapefile
  selected_shapefile <- reactiveVal(NULL)
  selected_shapefile_path <- reactiveVal(NULL)
  ########## stuff for drop-down menu - end   ##########  
  
  # Reactive value to track the state of drawing tools
  drawing_enabled <- reactiveVal(TRUE)
  
  # Initialize selected_shapefile with NULL
  selected_shapefile(NULL)
  
  # Add a reactive value to store the selected polygon
  selected_polygon <- reactiveVal(NULL)
  selected_polygon(NULL)
  
  # observers for polygon name text field info buttons
  observeEvent(input$info_polygon_name, {
    showModal(modalDialog(
      title = "Polygon Name Information",
      "The full name of the location you have selected, sans Island name, ex: 'Waikiki Watershed'.  This will be in the title of the resulting powerpoint.  Do not include the island name, as that will already be automatically inserted into the title.  If you selected a pre-defined shape, this will be filled in for you, but you can edit it if you wish.",
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  # observer for polygon short name text field info button
  observeEvent(input$info_polygon_short_name, {
    showModal(modalDialog(
      title = "Polygon Short Name Information",
      "A short name of the location you have selected, sans Island name, ex: 'Waikiki' instead of 'Waikiki Watershed'.  This will be used in the text of many slides in the resulting powerpoint.  If you selected a pre-defined shape, this will be filled in for you, but you can edit it if you wish.",
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  # Reactive expression for email validation
  valid_text_inputs <- reactive({
    polygon_name <- trimws(gsub("[^a-zA-Z0-9 _-āēīōūʻ]", "", input$polygon_name))
    polygon_short_name <- trimws(gsub("[^a-zA-Z0-9 _-āēīōūʻ]", "", input$polygon_short_name))
    email <- trimws(input$email)
    
    # Check if the email is valid and if polygon names are not empty
    print(paste0("valid_text_inputs: ", email, ", ", polygon_name, ", ", polygon_short_name))
    validate_text_inputs(input$email, polygon_name, polygon_short_name)
  })
  
  # Reactive expression to check if a polygon is selected or drawn
  polygon_selected <- reactive({
    print(paste0("polygon_selected: ", 
      !is.null(drawnFeature()), ", ", 
      !is.null(selected_polygon()), ", ", 
      !is.null(selected_shapefile())))
    
    !is.null(drawnFeature()) || !is.null(selected_polygon()) || !is.null(selected_shapefile())
  })
  
  # Render a leaflet map with drawing tools, initialization, only runs once upon startup
  output$map <- renderLeaflet({
    cat("renderLeaflet\n")
    
    # made this not use pipes so I could add debugging statements between actions
    leaf <- leaflet()
    leaf <- addTiles(leaf)
    
    # Add satellite view
    leaf <- addProviderTiles(leaf, providers$Esri.WorldImagery)
    
    # Set initial map view to the Hawaiian Islands
    leaf <- setView(leaf, lng = -157.4983, lat = 20.2927, zoom = 7)
    
    # Add island boundaries as a layer
    leaf <- addPolygons(leaf, 
                        data = island_boundaries,
                        group = "Islands",
                        fillColor = "transparent",
                        color = "#B0C4DE", # Choose a light color for the boundaries
                        weight = 1, # Set the line weight (thickness) to a low value
                        opacity = 0.5 # Set line opacity to a lower value
    )
    
    # # Enable layers control for toggling drawn shapes and island boundaries
    # leaf <- addLayersControl(leaf,
    #                          overlayGroups = c("drawnFeatures", "Islands"),
    #                          options = layersControlOptions(collapsed = FALSE)
    # )

    # Check if drawing is enabled and add drawing tools if it is
    if (drawing_enabled()) {
      leaf <- addDrawToolbar(
        leaf,
        targetGroup = "drawnFeatures",
        polyline = FALSE, # Exclude polyline tool
        marker = FALSE,   # Exclude marker tool
        circleMarker = FALSE, # Exclude circlemarker tool
        singleFeature = TRUE,
        editOptions = editToolbarOptions(
          selectedPathOptions = selectedPathOptions()
        )
      )
    }
    # Disable the save button
    updateActionButton(session, "save_button", disabled = TRUE)

    leaf
  })

  ########## stuff for drop-down menu - start   ##########  
  # Create a dropdown menu of shapefile names
  output$shapefile_select <- renderUI({
    cat("renderUI\n")

    # Create another list with the paths and '.shp' extension stripped
    shapefile_names <- sub("\\.shp$", "", basename(shapefile_paths))
    cat("shapefile_names: ", shapefile_names, "\n")
    
    # Merge shapefile_names and shapefile_choices into one list
    shapefile_list <- setNames(c("Select a pre-defined location", shapefile_paths), 
      c("Select a pre-defined location", shapefile_names))
    cat("shapefile_list: ", shapefile_list, "\n")
    
    selectInput(
      inputId = "shapefile", 
      label = "Select a pre-defined location", 
      choices = shapefile_list)
  })

  # Observe changes in email input
  observe({
    update_submit_button(session, valid_text_inputs(), polygon_selected())
  })
  
  # Update the selected_shapefile reactive value when a shapefile is selected
  observe({
    cat("observe 1\n")
    req(input$shapefile)
    if (input$shapefile != "Select a pre-defined location") {
      cat(" shapefile selected from menu\n")
      
      #orig_selected_shapefile = selected_shapefile()
      #cat("orig_selected_shapefile: ", orig_selected_shapefile(), "\n")

      #shapefile <- selected_shapefile()
      #first_attribute <- colnames(shapefile)[1]
      # cat("first_attribute: ", shapefile[[first_attribute]], "\n")
      
      
      selected_shapefile_path(normalizePath(shapefile_paths[match(input$shapefile, shapefile_paths)]))
      selected_shapefile(sf::st_read(shapefile_paths[match(input$shapefile, shapefile_paths)]))
      #polygon_name <- input$polygon_name
      #if (orig_selected_shapefile != selected_shapefile) {
      #  cat(" no match, so wiping text fields\n")
      #  updateTextInput(session, "polygon_name", value = "")
      #  updateTextInput(session, "polygon_short_name", value = "")
      #}
      
      
      # Update button label
      #updateActionButton(session, "save_button", label = "Generate Report", disabled = TRUE)
      
      # Disable drawing when a shapefile is selected
      drawing_enabled(FALSE)
      
      # reset all values so user starts from scratch
      drawnFeature(NULL)

    } else {
      cat(" no menu selection made, pt1\n")
      # if there's a feature being drawn, don't clear all this out, but if we are 
      # switching from selecting a pre-defined shape to drawing our own, clear all
      # the cruft away so user starts fresh.
      if (is.null(drawnFeature())) {      
        # JEN HERE: this code is getting executed when anything happens on the map or inputs, even if something is already drawn.
        cat(" no menu selection made, pt2\n")
        selected_shapefile(NULL)
        # empty out labels
        updateTextInput(session, "polygon_name", value = "")
        updateTextInput(session, "polygon_short_name", value = "")
        # Update button label
        updateActionButton(session, "save_button", label = "Generate Report", disabled = TRUE)
        # Enable drawing when no shapefile is selected
        drawing_enabled(TRUE)
        
        # reset all values so user starts from scratch
        drawnFeature(NULL)
        selected_shapefile <- reactiveVal(NULL)
        selected_shapefile_path <- reactiveVal(NULL)
        selected_shapefile(NULL)
        selected_polygon <- reactiveVal(NULL)
        selected_polygon(NULL)
      }
    }
    update_submit_button(session, valid_text_inputs(), polygon_selected())
  })
  
  # Render the selected shapefile on the map
  observe({
    cat("observe 2\n")
    
    # clear any existing status message
    output$status_messages <- renderText({
      ""
    })

    shapefile <- selected_shapefile()
    if (!is.null(shapefile)) {
      proxy <- leafletProxy("map")
      proxy <- clearShapes(proxy)
      proxy <- addPolygons(proxy,
                data = sf::st_transform(shapefile, crs = 4326), # Reproject to WGS84
                group = "Selected Shapefile",
                fillColor = "transparent",
                color = "#B0C4DE", # Choose a light color for the boundaries
                weight = 2, # Set the line weight (thickness) to a low value
                opacity = 0.5, # Set line opacity to a lower value
                layerId = shapefile[[names(shapefile)[1]]]
      )
      # # Auto-select and zoom to a single polygon if it exists
      # if (nrow(shapefile) == 1) {
      #   polygon_id <- shapefile[[names(shapefile)[1]]][1]
      #   selected_polygon(shapefile$geometry)
      #   # Update the map to highlight the selected polygon
      #   proxy <- clearGroup(proxy, "highlightedPolygon")
      #   proxy <- addPolygons(proxy, 
      #                        data = shapefile, 
      #                        group = "highlightedPolygon",
      #                        fillColor = "blue",
      #                        color = "#B0C4DE", # Choose a light color for the boundaries
      #                        weight = 2, # Set the line weight (thickness) to a low value
      #                        opacity = 0.5, # Set line opacity to a lower value
      #                        fillOpacity = 0.5, # Set fill opacity to make the fill color visible
      #                        layerId = ~get(names(shapefile)[1])
      #   )
      #   
      #   # Zoom to the polygon's bounds
      #   current_bounds <- proxy %>% leaflet::getMapBounds()
      #   zoom_level <- current_bounds$zoom
      #   proxy %>% leaflet::fitBounds(
      #     lng1 = sf::st_bbox(shapefile)[1],
      #     lat1 = sf::st_bbox(shapefile)[2],
      #     lng2 = sf::st_bbox(shapefile)[3],
      #     lat2 = sf::st_bbox(shapefile)[4],
      #     zoom = zoom_level + 1  # Add 1 to zoom level for a closer view
      #   )
      # }
    } else {
      proxy <- leafletProxy("map")
      proxy <- clearShapes(proxy)
    }
  })
  
  
  ########## stuff for drop-down menu - end   ##########  
  
  # Observe event for when shapes are drawn on the map
  observeEvent(input$map_draw_new_feature, {
    cat("observeEvent: map_draw_new_feature\n")
    
    # Update the reactive value with the drawn feature
    drawnFeature(input$map_draw_new_feature)
    
    # enable the submit button
    update_submit_button(session, valid_text_inputs(), polygon_selected())
    print(paste0("just called update_submit_button: ", valid_text_inputs(), ", ", polygon_selected()))
    # reset all values so user starts from scratch
    selected_shapefile <- reactiveVal(NULL)
    selected_shapefile_path <- reactiveVal(NULL)
    selected_shapefile(NULL)
    selected_polygon <- reactiveVal(NULL)
    selected_polygon(NULL)
  })
  
  
  # Observe click event on polygons
  observeEvent(input$map_shape_click, {
    cat("observeEvent: map_shape_click\n")
    
    # reset anything that was drawn on the map
    drawnFeature(NULL)
    
    # Get the ID of the clicked polygon
    polygon_id <- input$map_shape_click$id
    #cat("polygon_id: ", polygon_id, "\n")
    # Find the corresponding polygon from the selected shapefile
    shapefile <- selected_shapefile()
    if (!is.null(shapefile)) {
      # Get the name of the first attribute
      first_attribute <- colnames(shapefile)[1]
      #print(paste0("first_attribute: ", first_attribute))
      
      # Extract the selected polygon
      selected_polygon_data <- shapefile[shapefile[[first_attribute]] == polygon_id, ]
      selected_polygon(selected_polygon_data$geometry)
      #selected_poly <- selected_polygon()
      
      test_polygon_name <- selected_polygon_data$name
      test_polygon_short_name <- selected_polygon_data$short_name
      # Check if test_polygon_short_name is not set to a string of length longer than one character
      if (is.null(test_polygon_short_name) || nchar(test_polygon_short_name) <= 1) {
        test_polygon_short_name <- test_polygon_name
      }
      
      print(paste0("test_polygon_name: ", test_polygon_name))
      print(paste0("test_polygon_short_name: ", test_polygon_short_name))
      # Update the text inputs with the selected polygon data
      updateTextInput(session, "polygon_name", value = test_polygon_name)
      updateTextInput(session, "polygon_short_name", value = test_polygon_short_name)
      
      #print(colnames(shapefile))
      
      # Update the map to highlight the selected polygon
      proxy <- leafletProxy("map")
      
      # Remove any previous highlights
      proxy <- clearGroup(proxy, "highlightedPolygon")
      
      # Add the selected polygon with blue fill color
      proxy <- addPolygons(proxy, 
                           data = selected_polygon_data, 
                           group = "highlightedPolygon",
                           fillColor = "blue",
                           color = "#B0C4DE", # Choose a light color for the boundaries
                           weight = 2, # Set the line weight (thickness) to a low value
                           opacity = 0.5, # Set line opacity to a lower value
                           fillOpacity = 0.5, # Set fill opacity to make the fill color visible
                           layerId = ~get(first_attribute)
      )
    }
    update_submit_button(session, valid_text_inputs(), polygon_selected())
  })
  
  # Observe event for the save button click
  observeEvent(input$save_button, {
    cat("observeEvent: save_button\n")
    email <- input$email

    # Get the drawn feature's geometry from the reactive value
    feature <- drawnFeature()
    # Check if a feature was drawn
    if (!is.null(feature)) {
      polygon_name <- input$polygon_name
      polygon_short_name <- input$polygon_short_name
      
      cat("feature not null\n")
      # Verify the feature is of type "Feature" and has geometry
      if (feature$type == "Feature" && !is.null(feature$geometry)) {
        cat("feature is of type feature and has geometry\n")
        # Extract the type and coordinates from the geometry
        geometry_type <- feature$geometry$type
        coordinates <- feature$geometry$coordinates
        
        # Convert the geometry type and coordinates into an `sf` object
        sf_object <- NULL
        
        # Handle different geometry types
        if (geometry_type == "Polygon") {
          # Convert coordinates to a list of numeric lists
          if (length(coordinates) == 1 && length(coordinates[[1]]) > 0) {
            polygon_coords <- lapply(coordinates[[1]], function(coord) {
              if (length(coord) == 2) as.numeric(coord) else stop("Invalid Polygon coordinates")
            })
            # Convert the list of coordinates to a matrix
            polygon_coords_matrix <- do.call(rbind, polygon_coords)
            sf_object <- sf::st_sfc(sf::st_polygon(list(polygon_coords_matrix)), crs = 4326)
          } else {
            showNotification("Invalid shape. Please draw a fully enclosed shape", type = "error")
            return
          }
        } else {
          showNotification("Unsupported geometry type. Please draw a fully enclosed shape", type = "error")
          return
        }
        
        # Ensure CRS of the drawn feature and island boundaries match
        if (sf::st_crs(sf_object) != sf::st_crs(island_boundaries)) {
          sf_object <- sf::st_transform(sf_object, sf::st_crs(island_boundaries))
        }

        # Generate the filename using the current date and time
        datetime_str <- format(Sys.time(), "%Y_%m_%d_%H_%M_%S")
        filename <- paste0(polygon_short_name, "_", datetime_str)
        full_filename <- paste0("Shapefiles/UserDefinedPolygon/", filename, ".shp")
        #full_filepath <- paste0(BASE_DIR, "PDKESite/", "Shapefiles/UserDefinedPolygon/", filename, ".shp")

        # Write the sf object to a shapefile
        sf::st_write(sf_object, full_filename)
        #showNotification(paste("The selected area has been saved as:", full_filename), type = "message")
        
        # call the other script asynchronously to do the processing
        run_ccvd(session, sf_object, island_boundaries, full_filename, polygon_name, polygon_short_name, email)
        # Render the confirmation message when the submit button is clicked
        output$status_messages <- renderText({
          status_message
        })
        
        # Disable the save button to prevent a re-submission.
        updateActionButton(session, "save_button", disabled = TRUE)
      } else {
        showNotification("You first need to select an area.", type = "warning")
      }
    # user selected a specific polygon from a pre-defined shapefile
    } else if (!is.null(selected_polygon())) {
      print("selected_polygon save")
      tryCatch({
        # Extract the geometry object from the reactive value
        sf_object <- selected_polygon()
        #print(paste0("sf_object: ", sf_object))
        polygon_name <- input$polygon_name
        polygon_short_name <- input$polygon_short_name
        print(paste0("polygon_name: ", polygon_name))
        print(paste0("polygon_short_name: ", polygon_short_name))

        # Generate the filename using the current date and time
        datetime_str <- format(Sys.time(), "%Y_%m_%d_%H_%M_%S")
        filename <- paste0(polygon_name, "_" , datetime_str)
        full_filename <- paste0("Shapefiles/SelectedPolygon/", filename, ".shp")

        # Save the selected polygon as a shapefile
        sf::st_write(sf_object, full_filename)
        #showNotification(paste("The selected polygon has been saved as:", full_filename), type = "message")

        # verify the CRS of the two spatial objects match
        if (st_crs(sf_object) != st_crs(island_boundaries)) {
          # reproject one of the objects so that it has the same CRS as the other.
          island_boundaries <- st_transform(island_boundaries, st_crs(sf_object)) # Transform island_boundaries to the CRS of sf_object
        }
        
        # need this because was getting "Error: Error in wk_handle.wk_wkb(wkb, s2_geography_writer(oriented = oriented, : Loop 0 is not valid: Edge 1865 is degenerate (duplicate vertex)\n"
        sf_object <- st_make_valid(sf_object)
        
        run_ccvd(session, sf_object, island_boundaries, full_filename, polygon_name, polygon_short_name, email)

        # Render the confirmation message when the submit button is clicked
        output$status_messages <- renderText({
          status_message
        })
        
        # Disable the save button to prevent a re-submission.
        updateActionButton(session, "save_button", disabled = TRUE)
        
      }, error = function(e) {
        print(paste("Error:", e))
        showNotification("An error occurred during submission, please contact hcdp@hawaii.edu.", type = "error")
      })
    # user selected a pre-defined shapefile, really, this was early on and should never happen now.  
    # All predefined shapefiles have multiple polygons now.
    } else if (!is.null(selected_shapefile())) {
      print("selected_shapefile")
      polygon_name <- input$polygon_name
      polygon_short_name <- input$polygon_short_name
      
      full_filename <- sub(".*PDKESite/", "", selected_shapefile_path())
      filename <- sub(".*/(.*)\\.shp$", "\\1", selected_shapefile_path())
      run_ccvd(session, selected_shapefile(), island_boundaries, full_filename, polygon_name, polygon_short_name, email)
      #cat(file=stderr(), "selected_shapefile2: ", selected_shapefile_path(), "\n")
    } 
    # reset all values so user starts from scratch
    drawnFeature(NULL)
    selected_shapefile <- reactiveVal(NULL)
    selected_shapefile_path <- reactiveVal(NULL)
    selected_shapefile(NULL)
    selected_polygon <- reactiveVal(NULL)
    selected_polygon(NULL)
  }, ignoreInit = TRUE)
 
  # JavaScript to make the info circles trigger Shiny events
  session$onFlushed(function() {
    session$sendCustomMessage(type = 'bindInfoCircleClick', message = list(
      id = "info_polygon_name",
      event = "info_polygon_name"
    ))
    session$sendCustomMessage(type = 'bindInfoCircleClick', message = list(
      id = "info_polygon_short_name",
      event = "info_polygon_short_name"
    ))
  })

}