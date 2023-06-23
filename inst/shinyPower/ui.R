library(shiny)

inline_wrap <- function(f, ...){
  out <- f(...)
  div(style='display:inline-block; width: 100px; vertical-align:top', out)
}

ui <- fluidPage(
  tags$head(
    tags$style(HTML('#run{background-color:orange}'))
  ),
  titlePanel("Power Analysis"),
  sidebarLayout(
    sidebarPanel(width=4,
                 htmlOutput("mod"),
                 htmlOutput("class"),
                 htmlOutput("sites"),
                 br(),
                 inline_wrap(numericInput, inputId="alpha", label="Type I error (alpha)",
                             value=0.05, min=0.001, max=1),
                 inline_wrap(numericInput, inputId="nsims", label="Number of simulations",
                             value=10, min=1, max=300, step=1),
                 br(), br(),
                 #h3("Site scenarios"),
                 HTML("<b>Number of site (M) scenarios:</b>"),
                 inline_wrap(numericInput, inputId="ndesign_sites", label=NULL,
                             min=1, max=10, value=1, step=1),
                 uiOutput("design_sites"),

                 HTML("<b>Number of obs (J) scenarios:</b>"),
                 inline_wrap(numericInput, inputId="ndesign_obs", label=NULL,
                             min=1, max=10, value=1, step=1),
                 uiOutput("design_obs"),
                 #uiOutput("scenarios"),
                 br(),
                 uiOutput("coef_ui"),
                 br(),
                 actionButton("run", "Run analysis")
                 ),
    mainPanel(width=8,
              tabsetPanel(
                tabPanel("Summary", tableOutput("summary")),
                tabPanel("Plot",
                         uiOutput("param_selector"),
                         plotOutput("plot"))
              )
    )
  )
)

ui
