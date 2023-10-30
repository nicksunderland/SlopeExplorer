#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#' @import shiny
#' @importFrom shinyWidgets prettyRadioButtons sliderTextInput
#' @noRd
app_ui <- function(request) {

  # Create the page content
  fluidPage(
    # custom css
    tags$head(
      tags$style(HTML("hr {border-top: 1px solid #000000;}"))
    ),
    # Leave this function for adding external resources
    golem_add_external_resources(),
    # A title
    titlePanel("SlopeExplorer"),
    # Horizontal line
    tags$hr(),
    # sidebar layout
    sidebarLayout(position = "right",
                  sidebarPanel("Controls",
                               width=3,
                               prettyRadioButtons(inputId = "surface_type",
                                                  label = "Select plot:",
                                                  choiceNames = c("2D", "3D"),
                                                  choiceValues = c("2d", "3d"),
                                                  selected = "2d",
                                                  inline = TRUE,
                                                  animation = "smooth"),
                               prettyRadioButtons(inputId = "data_source",
                                                  label = "Select data:",
                                                  choiceNames = c("file", "example_1", "composite_1"),
                                                  choiceValues = c("file", "example_1", "composite_1"),
                                                  selected = "example_1",
                                                  inline = TRUE,
                                                  animation = "smooth"),
                               fileInput(inputId = "file_input", "Choose a file:"),
                               # distribution guesses / inits
                               hr(),
                               fluidRow(
                                 column(6, p(strong("SlopeHunter:"))),
                                 column(6, actionButton(inputId = "run_button", "Run"))
                               ),
                               fluidRow(
                                 column(4, numericInput(inputId="xp_thresh", HTML(paste0("p", tags$sub("i"))), value = 0.1, step=0.001)),
                                 column(4, numericInput(inputId="pi0", "\u03C00", value = 0.6, step=0.01)),
                                 column(4, numericInput(inputId="sxy1_slope_init", "Cov(i,p)", value = 0.00001, step=0.00001))
                               ),
                               sliderInput(inputId = "iter",
                                           "SlopeHunter step:",
                                           min = 0,
                                           max = 0,
                                           value = 0),
                               # distribution values
                               hr(),
                               p(strong("Distribution parameters:")),
                               p(HTML(paste0("G", tags$sub("i"), " - incidence only SNPs"))),
                               fluidRow(
                                 column(3, numericInput(inputId="ux0", HTML(paste0("\u03BC", tags$sub("i"))), value = 0, step=0.01)),
                                 column(3, numericInput(inputId="uy0", HTML(paste0("\u03BC", tags$sub("p"))), value = 0, step=0.01)),
                                 column(3, numericInput(inputId="sx0", HTML(paste0("\u03C3", tags$sub("i"))), value = 1, step=0.01, min=1e-7)),
                                 column(3, numericInput(inputId="sy0", HTML(paste0("\u03C3", tags$sub("p"))), value = 1, step=0.01, min=1e-7))
                               ),
                               fluidRow(
                                 column(3, numericInput(inputId="sxy0", "Cov(i,p)", value = 0, step=0.01)),
                                 column(3, numericInput(inputId="dir0", "Sign(cov)", value = 0, step=0.01))
                               ),
                               p(HTML(paste0("G", tags$sub("p"), " - incidence & progression SNPs"))),
                               fluidRow(
                                 column(3, numericInput(inputId="ux1", HTML(paste0("\u03BC", tags$sub("i"))), value = 0, step=0.01)),
                                 column(3, numericInput(inputId="uy1", HTML(paste0("\u03BC", tags$sub("p"))), value = 0, step=0.01)),
                                 column(3, numericInput(inputId="sx1", HTML(paste0("\u03C3", tags$sub("i"))), value = 1, step=0.01, min=1e-7)),
                                 column(3, numericInput(inputId="sy1", HTML(paste0("\u03C3", tags$sub("p"))), value = 1, step=0.01, min=1e-7))
                               ),
                               fluidRow(
                                 column(3, numericInput(inputId="sxy1", "Cov(i,p)", value = 0, step=0.01)),
                                 column(3, numericInput(inputId="dir1", "Sign(cov)", value = 0, step=0.01))
                               )

                  ),
                  mainPanel(p(strong("Display:")),
                            width=9,
                            plotlyOutput("main_plot", height = "800px")
                  )
     )
  )
}




textInputRow<-function (inputId, label, value = "")
{
  div(style="display:inline-block",
      tags$label(label, `for` = inputId),
      tags$input(id = inputId, type = "text", value = value,class="input-small"))
}












#' Add external Resources to the Application
#'
#' This function is internally used to add external
#' resources inside the Shiny application.
#'
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function() {
  add_resource_path(
    "www",
    app_sys("app/www")
  )

  tags$head(
    favicon(),
    bundle_resources(
      path = app_sys("app/www"),
      app_title = "SlopeExplorer"
    )
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert()
  )
}
