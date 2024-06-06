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
                        ADC In-Vivo: Cytosol Effect: Cytotoxic </div>\');
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
          ),
          numericInput("thresholdTumorVolume", "Threshold tumor volume (mm^3, 0 = no limit)", value = 0),
          numericInput("maxTime", "Max t (days)", value = 52),
          numericInput("doseT", "Dose given every x day", value = 0),
          numericInput("doseMaxT", "For a total number of times of (1 + x)", value = 0),
          numericInput("startT", "Dosing start time (h)", value = 0),
          checkboxInput("enableTmdd", "Add TMDD cell population")
        ),
        
        column(
          7,
          fillRow(
            width = 800,
            tabsetPanel(
              tabPanel(
                "Plots",
                plotlyOutput(outputId = "plot1", width = 600 * 1, height = 450 * 0.8),
                plotlyOutput(outputId = "plot2", width = 600 * 1, height = 450 * 0.8),
                plotlyOutput(outputId = "plot3", width = 600 * 1, height = 450 * 0.8),
                plotlyOutput(outputId = "plot4", width = 600 * 1, height = 450 * 0.8),
                br(),
                verbatimTextOutput("AUC") 
              ),
              tabPanel("Points", dataTableOutput("resultTable")),
              tabPanel("Fitting Experimental Data",
                       checkboxGroupInput("fitTheseParameters",
                                          "Fit these parameters:",
                                          choiceNames =
                                            list(
                                              HTML(paste0("EC50")),
                                              HTML(paste0("Ag", tags$sup("Ex"))),
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
              ),
              tabPanel("Parameter Variability", 
                       helpText("Both standard deviations must be set!"),
                       numericInput("std.EC50", "Standard deviation for EC50 (nM)", value = 0),
                       numericInput("std.Ag_cell_t", "Standard deviation for receptor count on each tumor cell Ag_cell_t (1 per cell)", value = 0),
                       numericInput("nSub", "nSub: Number between subject variabilities", value = 1),
                       helpText("Note that nSub < 2500 is needed to create plot!"),
                       numericInput("level", "level: % for confidence interval", value = 0.90),
                       helpText("Probability that the value of a parameter lies within this interval.
                                Thus, 100 - level % is the probalility that the true value lies outside the interval.
                                The interval ranges from min = (1-level)/2 to max = 1 - min."),
                       textAreaInput("expDataAsString2", 
                                     HTML("Input: time (days) &lt;space&gt; V_tumor (mm^3) &lt;space&gt; deviance <br/> (space delimited and comma as decimals)"),
                                     value = defaultExpDataAsString,
                                     rows = 6
                       ),
                       helpText("Simulation is calculated without threshold tumor volume."),
                       plotlyOutput(outputId = "plotVar", width = 600 * 1, height = 450 * 0.8)
              ),
              tabPanel("Waterfall Plots", 
                       wellPanel(
                         style = "background: #3c8dbc",
                         span("3D-plot options:", style = "color:white"),
                         textAreaInput(
                           "xAxisInstructions",
                           "Instructions for X-axis",
                           LoadDataDefault("default/defaultXAxisInstructions.json"),
                           rows = 5
                         ),
                         selectInput("xAxisItem", "Choose X-axis:",
                                     list(
                                       `Initial conditions` = LoadData("default/defaultConditions.tsv")[[1]],
                                       `Parameters` = LoadData("default/defaultParameters.tsv")[[1]]
                                     ),
                                     selected = "Ag_cell_t"
                         ),
                         textAreaInput(
                           "yAxisInstructions",
                           "Instructions for Y-axis",
                           LoadDataDefault("default/defaultYAxisInstructions.json"),
                           rows = 6
                         ),
                         selectInput("yAxisItem", "Choose Y-axis:",
                                     list(
                                       `Initial conditions` = LoadData("default/defaultConditions.tsv")[[1]],
                                       `Parameters` = LoadData("default/defaultParameters.tsv")[[1]]
                                     ),
                                     selected = "DAR"
                         ),
                         selectInput(
                           "zAxisItem", "Choose Z-axis:",
                           list("Tumor size (waterfall plot)", "Sum of longest diameter (waterfall plot)"),
                           "MED"
                         )
                       ),
                       tabsetPanel(
                         tabPanel(
                           "Graph",
                           plotlyOutput(outputId = "plot3d"),
                           textOutput("figureCaption")
                         ),
                         tabPanel("Points", dataTableOutput("resultTable3D"))
                       )
              )
            )
          )
        )
      )
    )
  )
))
