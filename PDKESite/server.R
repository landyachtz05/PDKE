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


library(shiny)
library(leaflet)
library(leaflet.extras)
library(sf)

# TODO: 
# - for polygon selection, make the line start/end markers less obnoxious
# - test selected polygon and created polygon on the PDKE stuff

environ <- "dev"
#environ <- "prod"
# prod, default if environ is not dev
rscript_path = "/bin/Rscript" 
PDKE_dir = "/home/exouser/workflow/"
#myscript_path = "/srv/shiny-server/sample-apps/PDKESite/test.R"
#myscript_path = paste0(PDKE_dir, "test.R")
myscript_path = paste0(PDKE_dir, "CCVD_portfolio_content.R")
if (environ == "dev") {
  rscript_path = "/Rscript"
  PDKE_dir = "/Users/jgeis/Work/PDKE/"
  #myscript_path = paste0(PDKE_dir, "test.R")
  myscript_path = paste0(PDKE_dir, "CCVD_portfolio_content.R")
}
run_string = paste(rscript_path, myscript_path)


# Function to validate email
validate_email <- function(email) {
  if (grepl("^[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\\.[A-Za-z]{2,}$", email)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}


display_intersected_islands <- function(intersected_island_names) {
  island_names_string <- ""
  # Display intersected island names
  if (length(intersected_island_names) > 0) {
    island_names_string <- paste(intersected_island_names, collapse = ", ")
    cat("\nisland_names_string: ", island_names_string, "\n")
    island_names_text <- paste("The selected area intersects with the following island(s):", paste(intersected_island_names, collapse = ", "))
    showNotification(island_names_text, type = "message")
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

run_ccvd <- function(sf_object, island_boundaries, shapefile_full_path, name, short_name, email) {
  # work out short name from island name
  #NP_DIR <- paste0(INPUTS_FOLDER, "waikiki_watershed/")
  #NP_FILE <- paste0(NP_DIR, "waikiki_watershed.shp")
  #NM <- "Waikiki Watershed"
  #NM_s <- "Waikiki"
  #ILE <- "Oahu"
  #ILE_s <- "OA"
  
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
  # this works
  system(paste0(Sys.getenv("R_HOME"), run_string, " ", shQuote(email), " ", shQuote(paste0(PDKE_dir, "PDKESite/", shapefile_full_path)), " ", shQuote(name), " ", shQuote(short_name), " ", shQuote(island_full_name), " ", shQuote(island_short_name)), wait = FALSE)
  showNotification("Background R script has been initiated.", type = "message") 
  
  # # Construct the command string
  # run_command <- paste0(
  #   'nohup ',
  #   Sys.getenv("R_HOME"), 
  #   run_string
  # )
  # cat("run_command: ", run_command, "\n")
  # args <- c(
  #   shQuote(paste0(PDKE_dir, "PDKESite/", shapefile_full_path)),
  #   shQuote(name),
  #   shQuote(short_name),
  #   shQuote(island_full_name),
  #   shQuote(island_short_name)
  # )
  # cat("args: ", args, "\n")
  # 
  # output_file = paste0("../CCVD/CCVD_OUTPUTS/", name, ".txt")
  # error_file =paste0("../CCVD/CCVD_OUTPUTS/", name, "_err.txt")
  # print(paste0("output_file: ", output_file))
  # print(paste0("error_file: ", error_file))
  # 
  # system2(run_command, args, stdout=output_file, stderr=error_file)

  
}
#nohup ./hello.sh > myoutput.txt >2&1 
#nohup sh -c '/Library/Frameworks/R.framework/Resources/Rscript /Users/jgeis/Work/PDKE/CCVD_portfolio_content.R '/Users/jgeis/Work/PDKE/PDKESite/Shapefiles/waikiki_watershed.shp' 'waikiki_watershed' 'waikiki_watershed' 'Oahu' 'OA'' > /dev/null 2>&1 &

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
  
  # Reactive expression for email validation
  valid_email <- reactive({
    validate_email(input$email)
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
    
    # Enable layers control for toggling drawn shapes and island boundaries
    leaf <- addLayersControl(leaf,
                             overlayGroups = c("drawnFeatures", "Islands"),
                             options = layersControlOptions(collapsed = FALSE)
    )

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
    shapefile_choices <- c("Select a pre-defined shapefile", shapefile_paths)
    selectInput("shapefile", "Select a shapefile:", choices = shapefile_choices)
  })
  
  # Update the selected_shapefile reactive value when a shapefile is selected
  observe({
    cat("observe 1\n")
    req(input$shapefile)
    if (input$shapefile != "Select a pre-defined shapefile") {
      cat(" shapefile selected from menu\n")
      
      selected_shapefile_path(normalizePath(shapefile_paths[match(input$shapefile, shapefile_paths)]))
      selected_shapefile(sf::st_read(shapefile_paths[match(input$shapefile, shapefile_paths)]))
      
      # Update button label
      updateActionButton(session, "save_button", label = "Generate data", disabled = TRUE)
      
      # Disable drawing when a shapefile is selected
      drawing_enabled(FALSE)
      
      # reset all values so user starts from scratch
      drawnFeature(NULL)

    } else {
      cat(" no menu selection made\n")
      
      selected_shapefile(NULL)
      # Update button label
      updateActionButton(session, "save_button", label = "Generate data", disabled = TRUE)
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
  })
  
  # Render the selected shapefile on the map
  observe({
    cat("observe 2\n")

    shapefile <- selected_shapefile()
    if (!is.null(shapefile)) {
      proxy <- leafletProxy("map")
      proxy <- clearShapes(proxy)
      proxy <- addPolygons(proxy,
                           data = shapefile,
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
  
  # # Enable submit button based on email input
  # observeEvent(input$email_input, {
  #   if (nchar(input$email_input) > 0) {
  #     updateActionButton(session, "save_button", label = "Generate Data", disabled = FALSE)
  #   } else {
  #     updateActionButton(session, "save_button", disabled = TRUE)
  #   }
  # })
  
  # Observe event for when shapes are drawn on the map
  observeEvent(input$map_draw_new_feature, {
    cat("observeEvent: map_draw_new_feature\n")
    
    # Update the reactive value with the drawn feature
    drawnFeature(input$map_draw_new_feature)
    
    # enable the submit button
    #updateActionButton(session, "save_button", label = "Generate data using drawn area", disabled = FALSE)
    updateActionButton(session, "save_button", label = "Generate data using drawn area", disabled = !valid_email())
    
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
      
      # Extract the selected polygon
      selected_polygon_data <- shapefile[shapefile[[first_attribute]] == polygon_id, ]
      selected_polygon(selected_polygon_data$geometry)
      selected_poly <- selected_polygon()
      
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
    #updateActionButton(session, "save_button", label = "Generate data", disabled = FALSE)
    updateActionButton(session, "save_button", label = "Generate data using drawn area", disabled = !valid_email())
  })
  
  # Observe event for the save button click
  observeEvent(input$save_button, {
    cat("observeEvent: save_button\n")
    email <- input$email

    # Get the drawn feature's geometry from the reactive value
    feature <- drawnFeature()
    # Check if a feature was drawn
    if (!is.null(feature)) {
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
            showNotification("Invalid Polygon coordinates. Please draw a valid shape.", type = "error")
            return
          }
        } else {
          showNotification("Unsupported geometry type. Please draw a valid shape.", type = "error")
          return
        }
        
        # Ensure CRS of the drawn feature and island boundaries match
        if (sf::st_crs(sf_object) != sf::st_crs(island_boundaries)) {
          sf_object <- sf::st_transform(sf_object, sf::st_crs(island_boundaries))
        }

        # Generate the filename using the current date and time
        datetime_str <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
        #filename <- paste0("Shapefiles/selected_area_", datetime_str, ".shp")
        filename <- paste0("selected_area_", datetime_str)
        full_filename <- paste0("Shapefiles/", filename, ".shp")
        #full_filepath <- paste0(PDKE_dir, "PDKESite/", "Shapefiles/", filename, ".shp")
        
        # Write the sf object to a shapefile
        sf::st_write(sf_object, full_filename)
        showNotification(paste("The selected area has been saved as:", full_filename), type = "message")
        
        # call the other script asynchronously to do the processing
        #intersected_island_names <- get_intersected_islands(sf_object, island_boundaries)
        #system(paste0(Sys.getenv("R_HOME"), run_string, " ", shQuote(full_filepath), " ", shQuote(intersected_island_names)), wait = FALSE)
        run_ccvd(sf_object, island_boundaries, full_filename, filename, datetime_str, email)
        
        showNotification("Background R script has been initiated.", type = "message")
        
      } else {
        showNotification("No area selected to save as shapefile.", type = "warning")
      }
    # user selected a specific polygon from a pre-defined shapefile
    } else if (!is.null(selected_polygon())) {
      print("selected_polygon save")
      tryCatch({
        # Extract the geometry object from the reactive value
        sf_object <- selected_polygon()

        # Generate the filename using the current date and time
        datetime_str <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
        filename <- paste0("selected_polygon_", datetime_str)
        full_filename <- paste0("Shapefiles/", filename, ".shp")
        
        # Save the selected polygon as a shapefile
        sf::st_write(sf_object, full_filename)
        showNotification(paste("The selected polygon has been saved as:", full_filename), type = "message")

        run_ccvd(sf_object, island_boundaries, full_filename, filename, datetime_str, email)
        
      }, error = function(e) {
        print(paste("Error:", e))
        showNotification("An error occurred while saving the selected polygon.", type = "error")
      })
    # user selected a pre-defined shapefile
    } else if (!is.null(selected_shapefile())) {
      print("selected_shapefile")
      
      full_filename <- sub(".*PDKESite/", "", selected_shapefile_path())
      filename <- sub(".*/(.*)\\.shp$", "\\1", selected_shapefile_path())
      run_ccvd(selected_shapefile(), island_boundaries, full_filename, filename, filename, email)
        
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
  
}