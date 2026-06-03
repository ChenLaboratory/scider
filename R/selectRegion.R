#' Select region of interest from plot
#'
#' @param data A data.frame object.
#' @param x.col Column name of the x coordinates.
#' @param y.col Column name of the y coordinates.
#' @param save.region Whether to also export the region defined by box/lasso 
#' select as an sf polygon.
#'
#' @return A data.frame object in the global environment. If save.region is TRUE,
#' output is a list with a data.frame of selected points and a sf polygon instead.
#' @export
#'
#' @examples
#'
#' data("xenium_bc_spe")
#'
#' spe_b <- spe[, SummarizedExperiment::colData(spe)$cell_type == "B cells"]
#'
#' dat <- as.data.frame(SpatialExperiment::spatialCoords(spe_b))
#'
#' # selectRegion(dat, x.col = "x_centroid", y.col = "y_centroid")
#'
selectRegion <- function(data, x.col = "x", y.col = "y",save.region=FALSE) {
    data <- as.data.frame(data)
    
    ui <- fluidPage(
        sidebarLayout(
            sidebarPanel(
                sliderInput("point_size", "Point Size:",
                    min = 1,
                    max = 10, value = 5
                ),
                selectInput("color_by", "Color Points by:",
                    choices = names(data)[!names(data) %in%c(x.col,y.col)],
                    selected = NULL
                ),
                actionButton("export_region", "Export Selected Points")
            ),
            mainPanel(
                style = "display: flex; flex-direction: column;
                align-items: center;",
                plotly::plotlyOutput("scatterplot", height = "60vh"),
                verbatimTextOutput("sel_points")
            )
        )
    )

    server <- function(input, output) {
        x <- reactiveVal(NULL) # points
        y <- reactiveVal(NULL) # region

        output$scatterplot <- renderPlotly({
            color <- data[[input$color_by]]
            colors=NULL
            if (is.null(color)) {
                legend = FALSE
            } else {
                legend = TRUE
                if (!is.numeric(color)) {
                    color = factor(color)
                    colors = grDevices::colorRampPalette(col.spec)(
                        nlevels(color))
                }
            }

            p <- plot_ly(data,
                x = ~ get(x.col), y = ~ get(y.col), type = "scatter",
                mode = "markers", marker = list(size = input$point_size),
                color = color,
                colors = colors
            )

            p <- layout(p,
                dragmode = "select",
                xaxis = list(title = "X"),
                yaxis = list(title = "Y"),
                showlegend = legend
            )
            
            # Speed up. https://plotly-r.com/performance 
            p <- plotly::toWebGL(p)
        })

        observeEvent(event_data("plotly_selected"), {
            i <- event_data("plotly_selected")$pointNumber + 1 # JS is 0-indexed
            sel_points <- data[i, ]
            x(sel_points)
        })
        
        if (save.region) {
            observeEvent(event_data("plotly_brushed"), {
                y(event_data("plotly_brushed"))
            })
        }

        output$sel_points <- renderPrint({
            x()
        })

        observeEvent(input$export_region, {
            sel_points <- x()
            region_coords <- y()
            if (!is.null(sel_points)) {
                # Reformat the selected points & region 
                sel_region <- as.data.frame(sel_points)
                if (save.region) {
                    sel_region <- list(points=sel_region)
                    sel_region$region <- plotly2sfpolygon(region_coords)
                }
                
                pos <- 1
                assign("sel_region", sel_region, envir = as.environment(pos))
                message("Selected region exported as 'sel_region'
                in the global environment.\n")
            }
        })
    }

    shinyApp(ui, server)
}

# Take coords returned by plotly_brushed event and convert them into sf polygon
plotly2sfpolygon <- function(coords) {
    if(is.null(coords)) return(NULL)
    
    n <- length(coords$x)
    if(n==2) { # box select (bounding box)
        res <- sf::st_polygon(list(
            cbind(coords$x[c(1,2,2,1,1)],coords$y[c(1,1,2,2,1)])
        ))
    } else { # lasso select
        # Closing polygon
        if (coords$x[1]!=coords$x[n]) {
            coords$x <- c(coords$x,coords$x[1])
            coords$y <- c(coords$y,coords$y[1])
        }
        res <- sf::st_polygon(list(
            cbind(coords$x,coords$y)
        ))
    }
    return (res)
}
