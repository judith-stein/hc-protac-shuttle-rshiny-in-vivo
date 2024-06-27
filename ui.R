# Define UI for the app ----
ui <- shinyUI(dashboardPage(
  dashboardHeader(title = "SimADC"),
  dashboardSidebar(
    actionButton(
      inputId = "linkBack",
      label = "Back to portal",
      onclick = "location.href='../index.html'"
    ),
    tags$h3("Databook"),
    helpText("Select components below to change their respective values in the fields on the right.
                    Double-check that values make sense. If a component has no value in the databook, no changes are made."),
    helpText(
      "Databook:",
      tags$a(href = substring(databookFileName, 5), target = "blank", databookTitle)
    ),
    div(
      style = "padding: 0px 0px; margin-bottom: -2em",
      selectInput("presetAb", "Antibody:", abList, selected = NULL)
    ),
    div(
      style = "padding: 0px 0px; margin-bottom: -2em",
      selectInput("presetCell", "Cell line:", cellList, selected = NULL)
    ),
    div(
      style = "padding: 0px 0px; margin-bottom: -2em",
      selectInput("presetAg", "Antigen:", agList, selected = NULL)
    ),
    div(
      style = "padding: 0px 0px; margin-bottom: -2em",
      selectInput("presetLinker", "Linker:", linkList, selected = NULL)
    ),
    selectInput("presetPL", "Payload:", payloadList, selected = NULL),
    helpText("Rows that turn red mean that the value is an assumption.", style = "color:red;"),
    tags$h3("Save or Load File"),
    helpText("Download current parameters or set new ones from a CSV file."),
    downloadButton("download_button", label = "Save to CSV"),
    tags$style(type = "text/css", "#download_button { margin-left: 15px; color: #444;}"),
    fileInput("savedSettingsFile", "Load from CSV",
              accept = c(
                "text/csv",
                "text/comma-separated-values,text/plain",
                ".csv"
              )
    ),
    helpText("This R Shiny App source code:"),
    downloadButton("downloadZip", label = "Download ZIP File"),
    tags$style(type = "text/css", "#downloadZip { margin-left: 15px; color: #444;}"),
    br(), br(),
    helpText("This R Shiny App is version:"),
    helpText(versionString),
    collapsed = TRUE
  ),
  dashboardBody(
    tags$script(HTML('$(document).ready(function() {
                      $("header").find("nav").append(\'<div class="myTopBar">\\
                        Protac Shuttle </div>\');
                      })')),
    
    # Include custom CSS to add top bar title and to change background color
    includeCSS("www/custom_style_to_dashboard.css"),
    fluidPage(
      
      # To higlight changes, needs shinyjs, also maybe V8
      useShinyjs(),
      extendShinyjs(text = jsCode, functions = c("backgroundCol")),
      
      # Sidebar layout with input and output definitions ----
      fluidRow(
        # Sidebar panel for inputs ----
        # width of this column is the first argument, sum of all has to be 12
        column(
          5,
          helpText(
            "Model equations: ",
            HTML("<a href=\"src/Model.R\" target=\"_blank\">Code</a>
                          <a href=\"logistic_pk_pd_equations_explained.pptx\" target=\"_blank\">PPTX</a>
                          <a href=\"logistic_pk_pd_equations_explained.pdf\" target=\"_blank\">PDF</a><br>")
          ),
          splitLayout(
            cellWidths = c("50%", "50%"),
            div(
          numericInput("maxTime", "Max t (days)", value = 52),
          selectInput("doseTo", "Dose of which PROxAb:",
                      list(
                        "Ab_C1_f", "Ab_C1_b1", "Ab_C1_b2", "Ab_C1_b3", "Ab_C1_b4"
                      ),
                      selected = "Ab_C1_b1"
          ),
          "bx <= max. number of bound protacs per Ab",
          numericInput("dose", "Dose (mg/kg)", value = 2),
          numericInput("doseDrug", "Extra dose of free Protac Drug_C1_f (mg/kg)", value = 0),
          "Dose converted to nmol/kg or nM with corresponding MW", br(),
          "Same dosing scheme for free Protac and Shuttle", 
            ),
          div(
          numericInput("doseT", "Dose given every x day", value = 0),
          numericInput("doseMaxT", "For a total number of times of (1 + x)", value = 0),
          numericInput("startT", "Dosing start time (h)", value = 0), 
          numericInput("max", "Max. number of bound protacs per antibody", value = 4),
          "min=1 and max=4"
          )
          ),
          br(),
          actionButton("update", "Calculate and plot"), br(), br(),
          tabsetPanel(
            id = "inputTablePanel",
            tabPanel(
              "Conditions",
              div(DTOutput("conditionDT"),
                  style = "font-size: 75%"
              )
            ),
            tabPanel(
              "Parameters",
              div(DTOutput("parameterDT"),
                  style = "font-size: 75%;
                                    margin-bottom: 50px"
              )
            ),
            selected = "Parameters"
          )
        ),
        
        column(
          7,
          fillRow(
            width = 800,
            tabsetPanel(
              tabPanel(
                "Plots",
                plotlyOutput(outputId = "plotPlasma", width = 600 * 1, height = 450 * 0.8),
                plotlyOutput(outputId = "plotCell", width = 600 * 1, height = 450 * 0.8),
                plotlyOutput(outputId = "plotTV", width = 600 * 1, height = 450 * 0.8)
              ),
              tabPanel("Points", dataTableOutput("resultTable")),
              tabPanel("Fitting Experimental Data",
                       checkboxGroupInput("fitTheseParameters",
                                          "Fit these parameters:",
                                          choiceNames =
                                            list(
                                              HTML(paste0("EC50")),
                                              HTML(paste0("Ag", tags$sup("cell"), tags$sub("t"))),
                                              HTML(paste0("tau"))
                                            ),
                                          inline = TRUE,
                                          choiceValues =
                                            list("EC50", "Ag_cell_t", "tau"),
                                          selected = c()
                       ),
                       checkboxGroupInput("fitToTheseColumns",
                                          "Fit model to these data:",
                                          choiceNames =
                                            list(
                                              HTML(paste0("V_tumor"))
                                            ),
                                          choiceValues =
                                            list("V_tumor_mm3"),
                                          selected = c("V_tumor_mm3")
                       ),
                       radioButtons("fitting",
                                    "Select type of parameter fitting:",
                                    choices = list("Local (fast)", "Global"),
                                    selected = c("Local (fast)")
                       ),
                       conditionalPanel(condition = "input.fitting == 'Global'",
                                        "The implementation searches between lower and upper for the global optimum", br(),
                                        "The optimization process will stop if either the maximum number of iterations or the error is reached",
                                        numericInput("lower", "Lower bound", value = 10^(-9)),
                                        numericInput("upper", "Upper bound", value = 1),
                                        numericInput("itermax", "Maximum number of iterations", value = 10),
                                        numericInput("VTR", "Error to be reached", value = 1e-02),
                       ),
                       # Output for parameter fitting, if fitting was done
                       verbatimTextOutput(outputId = "fittingResult"), 
                       textAreaInput("expDataAsString", 
                                     HTML("Input: time (days) &lt;space&gt; V_tumor (mm^3) &lt;space&gt; deviance <br/> (space delimited and comma as decimals)"),
                                     value = defaultExpDataAsString,
                                     rows = 6
                       ),
                       plotlyOutput(outputId = "plotExpData", width = 800 * 1, height = 600 * 0.8)
              )
            )
          )
        )
      )
    )
  )
))
