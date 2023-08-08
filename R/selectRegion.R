#' Select region of interest from plot
#'
#' @param data A data.frame object.
#' @param xcol Column name of the x coordinates.
#' @param ycol Column name of the y coordinates.
#'
#' @return A matrix.
#' @export
#'
#' @examples 
#' 
#' data("xenium_bc_spe")
#' 
#' spe_b <- spe[,colData(spe)$cell_type == "B cells"]
#' 
#' dat <- as.data.frame(spatialCoords(spe_b))
#' 
#' selectRegion(dat, x_col = "x_centroid", y_col = "y_centroid")
#' 
selectRegion <- function(data, x_col = "x", y_col = "y") {
  
  ui <- fluidPage(
    sidebarLayout(
      sidebarPanel(
        sliderInput("point_size", "Point Size:", min = 1, max = 10, value = 5),
        selectInput("color_by", "Color Points by:", choices = names(data), selected = NULL),
        actionButton("export_region", "Export Selected Points")
      ),
      mainPanel(
        style = "display: flex; flex-direction: column; align-items: center;",
        plotlyOutput("scatterplot", height = "60vh"),
        verbatimTextOutput("sel_points")
      )
    )
  )
  
  server <- function(input, output) {
    x <- reactiveVal(NULL)
    
    output$scatterplot <- renderPlotly({
      color_var <- input$color_by
      
      p <- plot_ly(data, x = ~get(x_col), y = ~get(y_col), type = "scatter", 
                   mode = "markers", marker = list(size = input$point_size))
      
      if (!is.null(color_var)) {
        if (is.numeric(data[[color_var]])) {
          p <- add_markers(p, color = ~get(color_var))
        } else {
          color_palette <- scales::hue_pal()(nlevels(factor(data[[color_var]])))
          p <- add_markers(p, color = ~factor(data[[color_var]]), colors = color_palette)
        }
      }
      
      p <- layout(p, dragmode = "select", 
                       xaxis = list(title = "X"), 
                       yaxis = list(title = "Y"),
                       showlegend = FALSE)
    })
    
    observeEvent(event_data("plotly_selected"), {
      sel_indices <- event_data("plotly_selected")$pointNumber
      sel_points <- data[sel_indices, ]
      x(sel_points)
    })
    
    output$sel_points <- renderPrint({
      x()
    })
    
    observeEvent(input$export_region, {
      sel_points <- x()
      if (!is.null(sel_points)) {
        sel_region <- as.matrix(sel_points)
        assign("sel_region", sel_region, envir = .GlobalEnv)
        cat("Selected region exported as 'sel_region' in the global environment.\n")
      }
    })
  }
  
  shinyApp(ui, server)
}
