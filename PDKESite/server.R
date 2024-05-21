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
# - don't enable the submit button until a valid polygon has been selected on the map
# - for polygon selection, make the line start/end markers less obnoxious
# - test selected polygon and created polygon on the PDKE stuff
# - add highlighting of selected polygon
# - when loading existing shapefile, make shape highlighting more obvious 
#   (actually, don't think I should do this as single-polygon shapefiles are unlikely to ever actually be used)

environ <- "dev"
#environ <- "prod"
# prod, default if environ is not dev
rscript_path = "/bin/Rscript" 
myscript_path = "/srv/shiny-server/sample-apps/PDKESite/test.R"
if (environ == "dev") {
  rscript_path = "/Rscript"
  myscript_path = "/Users/jgeis/Work/PDKE/test.R"
}
run_string = paste(rscript_path, myscript_path)

get_intersected_islands <- function(sf_object, island_boundaries) {
  # Check for intersection with island boundaries
  intersected_islands <- sf::st_intersects(sf_object, island_boundaries)
  
  # Extract island names or identifiers
  intersected_island_names <- unique(island_boundaries$isle[intersected_islands[[1]]])
  display_intersected_islands(intersected_island_names)
  return(intersected_island_names)
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
  
  # Render a leaflet map with drawing tools, initialization, only runs once upon startup
  output$map <- renderLeaflet({
    cat("renderLeaflet\n")
    
    # made this not use pipes so I could add debugging statements between actions
    leaf <- leaflet()
    leaf <- addTiles(leaf)
    
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
      selected_shapefile_path(normalizePath(shapefile_paths[match(input$shapefile, shapefile_paths)]))
      selected_shapefile(sf::st_read(shapefile_paths[match(input$shapefile, shapefile_paths)]))
      
      # Update button label
      updateActionButton(session, "save_button", label = "Generate data using selected shapefile")
      
      # Disable drawing when a shapefile is selected
      drawing_enabled(FALSE)
    } else {
      selected_shapefile(NULL)
      # Update button label
      updateActionButton(session, "save_button", label = "Generate data using selected area")
      # Enable drawing when no shapefile is selected
      drawing_enabled(TRUE)
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
  })
  
  
  # Observe click event on polygons
  observeEvent(input$map_shape_click, {
    cat("observeEvent: map_shape_click\n")
    
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
  })
  
  # Observe event for the save button click
  observeEvent(input$save_button, {
    cat("observeEvent: save_button\n")

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
        intersected_island_names <- get_intersected_islands(sf_object, island_boundaries)

        # Generate the filename using the current date and time
        datetime_str <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
        filename <- paste0("Shapefiles/selected_area_", datetime_str, ".shp")
        
        # Write the sf object to a shapefile
        sf::st_write(sf_object, filename)
        showNotification(paste("The selected area has been saved as:", filename), type = "message")
        
        # call the other script asynchronously to do the processing
        system(paste0(Sys.getenv("R_HOME"), run_string, " ", shQuote(filename), " ", shQuote(intersected_island_names)), wait = FALSE)
        
        showNotification("Background R script has been initiated.", type = "message")
        
      } else {
        showNotification("No area selected to save as shapefile.", type = "warning")
      }
    # user selected a specific polygon from a pre-defined shapefile
    } else if (!is.null(selected_polygon())) {
      print("selected_polygon save")
      # Generate the filename using the current date and time
      datetime_str <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
      filename <- paste0("Shapefiles/selected_polygon_", datetime_str, ".shp")
      tryCatch({
        # Extract the geometry object from the reactive value
        sf_object <- selected_polygon()

        # Save the selected polygon as a shapefile
        sf::st_write(sf_object, filename)
        showNotification(paste("The selected polygon has been saved as:", filename), type = "message")
        
        intersected_island_names <- get_intersected_islands(sf_object, island_boundaries)
        system(paste0(Sys.getenv("R_HOME"), run_string, " ", shQuote(filename), " ", shQuote(intersected_island_names)), wait = FALSE)
        
      }, error = function(e) {
        print(paste("Error:", e))
        showNotification("An error occurred while saving the selected polygon.", type = "error")
      })
    # user selected a pre-defined shapefile
    } else if (!is.null(selected_shapefile())) {
      # TODO: get isle from it's location, not the attribute as not all shapes will have isle attached
      print("selected_shapefile")
      sf_object<-selected_shapefile()
      intersected_island_names <- get_intersected_islands(sf_object, island_boundaries)
      system(paste0(Sys.getenv("R_HOME"), run_string, " ", shQuote(selected_shapefile_path()), " ", shQuote(intersected_island_names)), wait = FALSE)
      showNotification("Background R script has been initiated.", type = "message")
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