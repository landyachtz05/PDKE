#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#
# using https://maps.equatorstudios.com to test if saved shapefiles are legit

library(shiny)
library(leaflet)
library(leaflet.extras)
library(sf)

# dev
#rscript_path = "/Rscript"
#myscript_path = "/Users/jgeis/Work/PDKE/test.R"

# prod
rscript_path = "/bin/Rscript" 
myscript_path = "/srv/shiny-server/sample-apps/PDKESite/test.R"
run_string = paste(rscript_path, myscript_path)

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
    
    # Render a leaflet map with drawing tools, initialization, only runs once upon startup
    output$map <- renderLeaflet({
        leaf <- leaflet() %>%
            addTiles() %>%
            # Set initial map view to the Hawaiian Islands
            setView(lng = -157.4983, lat = 20.2927, zoom = 7) %>%
            # Add island boundaries as a layer
            addPolygons(data = island_boundaries, 
                        group = "Islands", 
                        fillColor = "transparent",
                        color = "#B0C4DE", # Choose a light color for the boundaries
                        weight = 1, # Set the line weight (thickness) to a low value
                        opacity = 0.5, # Set line opacity to a lower value
            ) %>%
            # Enable layers control for toggling drawn shapes and island boundaries
            addLayersControl(
                overlayGroups = c("drawnFeatures", "Islands"),
                options = layersControlOptions(collapsed = FALSE)
            )
            if (drawing_enabled())
              leaf <- leaf %>% 
              # Enable drawing tools if drawing is enabled
              addDrawToolbar(
                targetGroup = "drawnFeatures",
                singleFeature = TRUE,
                editOptions = editToolbarOptions(
                  selectedPathOptions = selectedPathOptions()
                )
              )
            leaf
    })
    
    ########## stuff for drop-down menu - start   ##########  
    # Create a dropdown menu of shapefile names
    output$shapefile_select <- renderUI({
        shapefile_choices <- c("Select a pre-defined shapefile", shapefile_paths)
        selectInput("shapefile", "Select a shapefile:", choices = shapefile_choices)
    })
    
    # Update the selected_shapefile reactive value when a shapefile is selected
    observe({
        req(input$shapefile)
        if (input$shapefile != "Select a pre-defined shapefile") {
            selected_shapefile_path(normalizePath(shapefile_paths[match(input$shapefile, shapefile_paths)]))
            selected_shapefile(sf::st_read(shapefile_paths[match(input$shapefile, shapefile_paths)]))
            
            # Update button label
            updateActionButton(session, "save_button", label = "Generate data using selected shapefile")
            
            #cat("Selected Shapefile1:", selected_shapefile_path(), "\n")
            #shape<-selected_shapefile()
            #cat("Selected info: ", names(shape), "\n")
            #cat("isle: ", shape$isle, "\n")

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
        shapefile <- selected_shapefile()
        if (!is.null(shapefile)) {
            leafletProxy("map") %>%
                clearShapes() %>%
                addPolygons(data = shapefile, group = "Selected Shapefile")
        } else {
          leafletProxy("map") %>%
            clearShapes()
        }
    })
    ########## stuff for drop-down menu - end   ##########  
    
    # Observe event for when shapes are drawn on the map
    observeEvent(input$map_draw_new_feature, {
        # Update the reactive value with the drawn feature
        drawnFeature(input$map_draw_new_feature)
    })
    
    # Observe event for the save button click
    observeEvent(input$save_button, {

        # Get the drawn feature's geometry from the reactive value
        feature <- drawnFeature()
        # Check if a feature was drawn
        if (!is.null(feature)) {
            # Verify the feature is of type "Feature" and has geometry
            if (feature$type == "Feature" && !is.null(feature$geometry)) {
                # Extract the type and coordinates from the geometry
                geometry_type <- feature$geometry$type
                coordinates <- feature$geometry$coordinates
                
                # Convert the geometry type and coordinates into an `sf` object
                sf_object <- NULL
                
                # Handle different geometry types
                if (geometry_type == "Point") {
                    # Convert coordinates to numeric and check length
                    if (is.numeric(coordinates) && length(coordinates) == 2) {
                        sf_object <- sf::st_sfc(sf::st_point(coordinates), crs = 4326)
                    } else {
                        showNotification("Invalid Point coordinates. Please draw a valid shape.", type = "error")
                        return
                    }
                } else if (geometry_type == "Polygon") {
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
                } else if (geometry_type == "LineString") {
                    # Convert coordinates to a list of numeric lists
                    if (length(coordinates) > 0) {
                        linestring_coords <- lapply(coordinates, function(coord) {
                            if (length(coord) == 2) as.numeric(coord) else stop("Invalid LineString coordinates")
                        })
                        # Convert the list of coordinates to a matrix
                        linestring_coords_matrix <- do.call(rbind, linestring_coords)
                        sf_object <- sf::st_sfc(sf::st_linestring(linestring_coords_matrix), crs = 4326)
                    } else {
                        showNotification("Invalid LineString coordinates. Please draw a valid shape.", type = "error")
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
                
                # Check for intersection with island boundaries
                intersected_islands <- sf::st_intersects(sf_object, island_boundaries)
                
                # Extract island names or identifiers
                intersected_island_names <- unique(island_boundaries$isle[intersected_islands[[1]]])
                island_names_string <- ""
                # Display intersected island names
                if (length(intersected_island_names) > 0) {
                    island_names_string <- paste(intersected_island_names, collapse = ", ")
                    cat("\nisland_names_string: ", island_names_string, "\n")
                    #sf_object$island <- island_names_string # this runs, but the program crashes when it tries to write out the shapefile
                    island_names_text <- paste("The selected area intersects with the following island(s):", paste(intersected_island_names, collapse = ", "))
                    showNotification(island_names_text, type = "message")
                } else {
                    showNotification("The selected area does not intersect with any islands.", type = "warning")
                }
                
                # Generate the filename using the current date and time
                datetime_str <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
                filename <- paste0("Shapefiles/selected_area_", datetime_str, ".shp")
                
                # Save the sf_object as a shapefile
                #sf::st_write(sf_object, filename)
                #showNotification(paste("The selected area has been saved as:", filename), type = "message")
                
                # Create a data frame with island name and geometry
                shapefile_data <- data.frame(isle = island_names_string, geometry = sf_object)
                
                # Convert the data frame to an sf object
                sf_object <- sf::st_as_sf(shapefile_data)
                
                # Write the sf object to a shapefile
                sf::st_write(sf_object, filename)
                showNotification(paste("The selected area has been saved as:", filename), type = "message")
                
                # Check if sf_object is not NULL and has geometry
                #cat("is.null: ", !is.null(sf_object), "\n")
                #cat("geometry: ", !is.null(sf_object$geometry), "\n")
                #print(paste("sf_object: ", sf_object))
                
                #print(sf_object)
                #print(paste("filename: ", filename))
                #print(str(sf_object))
                #print(str(sf_object[[1]]))
                #print(paste("postcall"))
                
                #sf::st_write(do.call(sf::st_sfc, sf_object), filename)
                #showNotification(paste("The selected area has been saved as:", filename), type = "message")
                
                # call the other script asynchronously to do the processing
                system(paste0(Sys.getenv("R_HOME"), run_string, " ", shQuote(filename), " ", shQuote(intersected_island_names)), wait = FALSE)

                showNotification("Background R script has been initiated.", type = "message")
                
            } else {
                showNotification("No area selected to save as shapefile.", type = "warning")
            }
        } else if (!is.null(selected_shapefile())) {
          shape<-selected_shapefile()
          #cat(file=stderr(), "isle: ", shape$isle, "\n")
          system(paste0(Sys.getenv("R_HOME"), run_string, " ", shQuote(selected_shapefile_path()), " ", shape$isle), wait = FALSE)
          showNotification("Background R script has been initiated.", type = "message")
          #cat(file=stderr(), "selected_shapefile2: ", selected_shapefile_path(), "\n")
        }
    })
}