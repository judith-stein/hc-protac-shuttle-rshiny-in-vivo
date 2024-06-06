# Define server logic required to solve ODE's ----
server <- function(input, output, session) {
  
  #
  # Input DT
  #
  
  options(DT.options = list(pageLength = 100, dom = "t"))
  
  # Editable table for conditions:
  conditionTable <- LoadData("default/defaultConditions.tsv")
  names(conditionTable) <- c("Initial", "Value", "Unit", "Definition")
  output$conditionDT <- renderDT(conditionTable,
                                 selection = "none",
                                 rownames = FALSE, editable = TRUE
  )
  
  proxyConditionDT <- dataTableProxy("conditionDT")
  
  observeEvent(input$conditionDT_cell_edit, {
    info <- input$conditionDT_cell_edit
    str(info)
    i <- info$row
    j <- info$col + 1 # column index offset by 1
    v <- info$value
    if (j == 2) {
      conditionTable[i, j] <<- DT::coerceValue(v, conditionTable[i, j])
      replaceData(proxyConditionDT, conditionTable, resetPaging = FALSE, rownames = FALSE)
    }
  })
  
  # Editable table for parameters:
  parameterTable <- LoadData("default/defaultParameters.tsv")
  names(parameterTable) <- c("Parameter", "Value", "Unit", "Definition")
  
  # To highlight certain rows in the table, we add a new column with 0 or 1,
  # that is hidden from the view but denotes if a row should be highlighted or not
  parameterTable$background_color <- 0
  hideCols <- which(names(parameterTable) %in% c("background_color")) - 1
  parameterDatatable <- datatable(parameterTable,
                                  selection = "none",
                                  rownames = FALSE, editable = TRUE,
                                  options = list(columnDefs = list(list(visible = FALSE, targets = hideCols)))
  ) %>%
    formatStyle("background_color",
                target = "row",
                backgroundColor = styleEqual(c(0, 1), c("", "red"))
    )
  
  output$parameterDT <- renderDT(parameterDatatable)
  
  proxyParameterDT <- dataTableProxy("parameterDT")
  
  observeEvent(input$parameterDT_cell_edit, {
    info <- input$parameterDT_cell_edit
    str(info)
    i <- info$row
    j <- info$col + 1 # column index offset by 1
    v <- info$value
    if (j == 2) {
      parameterTable[i, j] <<- DT::coerceValue(v, parameterTable[i, j])
      parameterTable[i, 5] <<- 0 # Set red color back to normal
      replaceData(proxyParameterDT, parameterTable,
                  resetPaging = FALSE, rownames = FALSE
      )
    }
  })
  
  # -----------------------------
  # Preset feature from SIDEBAR
  
  UpdateDTTable <- function(name, newValue, isRed = FALSE) {
    if (name %in% conditionTable$Initial) {
      conditionTable[which(conditionTable$Initial == name), "Value"] <<- DT::coerceValue(
        newValue,
        conditionTable[which(conditionTable$Initial == name), "Value"]
      )
      replaceData(proxyConditionDT, conditionTable, resetPaging = FALSE, rownames = FALSE)
    } else if (name %in% parameterTable$Parameter) {
      parameterTable[which(parameterTable$Parameter == name), "Value"] <<- DT::coerceValue(
        newValue,
        parameterTable[which(parameterTable$Parameter == name), "Value"]
      )
      # Set row color to red if toggled
      indexOfRowToBeRed <- match(name, parameterTable$Parameter)
      parameterTable[indexOfRowToBeRed, "background_color"] <<- if (isRed) 1 else 0
      
      replaceData(proxyParameterDT, parameterTable, resetPaging = FALSE, rownames = FALSE)
    }
  }
  
  observeEvent(
    {
      input$presetAb
      input$presetCell
      input$presetAg
      input$presetLinker
      input$presetPL
    },
    {
      # When preset is chosen, update the relevant fields
      
      UpdateFromSheet <- function(idSheet, idCol, newId, varRow, listOfVars, foreignKey = NULL) {
        if (newId == "") {
          return()
        }
        
        rowforThisId <- tx %>%
          filter(sheet == idSheet, col == idCol, character == newId) %>%
          pull(row)
        for (var in listOfVars) {
          sheetVar <- var
          
          if (var == "Ag_cell_t") {
            if (foreignKey != "") {
              sheetVar <- foreignKey
            } else {
              next
            }
          }
          
          colforThisVar <- tx %>%
            filter(sheet == idSheet, row == varRow, character == sheetVar) %>%
            pull(col)
          varCell <- tx %>%
            filter(sheet == idSheet, row == rowforThisId, col == colforThisVar)
          
          if (var %in% c("K_ADC_cell_int", "EC50") && foreignKey != "") {
            colForThisCellLine <- tx %>%
              filter(sheet == idSheet, row == varRow, character == foreignKey) %>%
              pull(col)
            specificVarCell <- tx %>%
              filter(sheet == idSheet, row == rowforThisId, col == colForThisCellLine)
            specificVar <- pull(specificVarCell, numeric)
            if (!identical(specificVar, numeric(0))) {
              varCell <- specificVarCell
            }
          }
          
          # Check that excel cell value is either numeric (preferred)
          # or something like "3.12*10^-12" or, if empty, then do nothing
          valueNumeric <- pull(varCell, numeric)
          valueCharacter <- pull(varCell, character)
          if (!is.na(valueNumeric) && (length(valueNumeric) != 0)) {
            value <- valueNumeric
          } else if (!is.na(valueCharacter) && (length(valueCharacter) != 0)) {
            value <- valueCharacter
          } else {
            next
          }
          
          # Is the cell red in the excel file? Can be also NA if cell font color is "Automatic"
          isRed <- (txFormat$local$font$color$"rgb"[varCell$local_format_id] == "FFFF0000")
          isRed <- if (is.na(isRed)) FALSE else isRed
          
          UpdateDTTable(var, value, isRed)
        }
      }
      
      newAbId <- input$presetAb
      newCellLineId <- input$presetCell
      newAgId <- input$presetAg
      newLinkerId <- input$presetLinker
      newPL <- input$presetPL
      
      UpdateFromSheet("Antibodies", 2, newAbId, 15, c(
        "K_ADC_cell_ag_on", "K_ADC_cell_ag_off", "CL_ADC", "CLD_ADC",
        "V_C1_ADC", "V_C2_ADC", "K_ADC_cell_int", "K_ADC_cell_lyso_pino", "K_all_cell_lyso_exo"
      ), newCellLineId)
      UpdateFromSheet("Cell Lines", 2, newCellLineId, 3, c("V_cell", "V_cell_lyso", "V_cell_nuc", "DT_tumor", "Ag_cell_t"), newAgId)
      UpdateFromSheet("Linkers", 2, newLinkerId, 3, c("K_ADC_deg", "K_ADC_dec"))
      UpdateFromSheet("Payloads", 2, newPL, 3, c(
        "K_Drug_ex_in", "K_Drug_ex_out",
        "K_Drug_lyso_in", "K_Drug_lyso_out", "K_Nuc_in", "K_Nuc_out", "Drug_Target_cell_cyto_t",
        "Binding_sites", "K_Drug_cyto_dt_on", "K_Drug_cyto_dt_off", "CL_Drug",
        "CLD_Drug", "V_C1_Drug", "V_C2_Drug", "P_Drug", "D_Drug", "f_ex_ub",
        "tau", "k_kill_max", "k_g", "k_r", "k_z", "EC50", "K_Drug_all_charged", "f_Drug_charged_ph4", "f_Drug_charged_ph7"
      ), newCellLineId)
    },
    ignoreInit = TRUE
  )
  
  
  # Preset feature end
  # -----------------------------
  
  #
  # Save State to file feature
  #
  
  output$download_button <- downloadHandler(
    # Save tables and other relevant (visible) input parameters
    # to a csv file. The file has #-comments and empty lines and
    # structure of: name,value
    filename = function() {
      # File name eg. savefile_adc-pk-pd_2020-03-27.csv
      paste0(
        "savefile_", basename(getwd()), "_", Sys.Date(), "_",
        versionNumbers[1], "-", versionNumbers[2], ".csv"
      )
    },
    content = function(file) {
      # Add app name, explanation and download date as comments to first lines
      write(paste("# Save file for", basename(getwd())), file, append = TRUE)
      write(paste("# Version:", versionString), file, append = TRUE)
      write("# Each non-commented line: name,value", file, append = TRUE)
      write(paste("# Downloaded", Sys.Date()), file, append = TRUE)
      
      # Add conditionTable, parameterTable and inputs as lines: name,value
      write("\n# Conditions", file, append = TRUE)
      write.table(conditionTable %>% select(Initial, Value), file, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",", append = T)
      
      write("\n# Parameters", file, append = TRUE)
      write.table(parameterTable %>% select(Parameter, Value), file, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",", append = T)
      
      write("\n# Other inputs", file, append = TRUE)
      for (name in c("thresholdTumorVolume", "maxTime", "doseT", "doseMaxT", "startT", "enableTmdd", "fitting",
                     "std.EC50", "std.Ag_cell_t", "nSub", "level")) {
        # Hardcoded because no distinction between important and unimportant shiny-inputs
        write(paste(name, input[[name]], sep = ","), file, append = TRUE)
      }
      if (input$fitting == "Global") {
        for (name in c("lower", "upper", "itermax", "VTR")) {
          # Hardcoded because no distinction between important and unimportant shiny-inputs
          write(paste(name, input[[name]], sep = ","), file, append = TRUE)
        }
      }
      
      for (name in c("expDataAsString", "expDataAsString2")) {
        write(paste0("\n# ", name), file, append = TRUE)
        write(paste(input[[name]], sep = ","), file, append = TRUE)
      }
    }
  )
  
  testUpload <- observeEvent(input$savedSettingsFile, {
    # When user uploads a file, read it as csv to df, loop through the rows,
    # turn safely to numeric values and update fields in the UI accordingly
    upfile <- input$savedSettingsFile
    if (!is.null(upfile)) {
      
      
      fileAsString <- readChar(upfile$datapath, upfile$size)
      inputSetStrings <- strsplit(fileAsString, "#")[[1]]
      
      for (s in inputSetStrings) {
        if (startsWith(s, " expDataAsString")) { # || startsWith(s, " expDataAsString2")
          test <- strsplit(s, "\n")
          content <- test[[1]][-1]
          text <- content[1]
          if (length(content) > 1){
            for (i in 2:length(content)) {
              text <- paste(text, content[i], sep = "\n")
            }
          }
          if (startsWith(s, " expDataAsString2")) {
            updateTextAreaInput(session, "expDataAsString2", value = text)
          } else {
            updateTextAreaInput(session, "expDataAsString", value = text)
          }
        } else {
          upDataFrame <- read.csv(
            text = s,
            header = FALSE,
            comment.char = "#",
            blank.lines.skip = TRUE,
            col.names = c("name", "rawValue"),
            colClasses = c("character", "character")
          )
          
          if (dim(upDataFrame)[1] > 0) {
            for (rowIndex in 1:nrow(upDataFrame)) {
              name <- upDataFrame[rowIndex, "name"]
              rawValue <- upDataFrame[rowIndex, "rawValue"]
              # eg. "10^9/(6.023*10^23)" -> 1.660302e-15
              numericValue <- StringToNumeric(rawValue)
              UpdateDTTable(name, numericValue) # Updates tables
              if (name %in% c("fitting")) {
                updateRadioButtons(session, inputId = name, selected = c(rawValue))
              } else if (name %in% names(input)) {
                # For other inputs (not in tables but eg. numericInput shiny elements)
                updateTextInput(session, name, value = numericValue)
              }
            }
          }
        }
      }
    }
  })
  
  
  output$downloadZip <- downloadHandler(
    # This button downloads all the files in this repository/app as a zip
    filename = function() {
      paste0(
        "source_", basename(getwd()), "_", Sys.Date(), "_",
        versionNumbers[1], "-", versionNumbers[2], ".zip"
      )
    },
    content = function(file) {
      # Download all source files as zip, excluding .git
      zip(zipfile = file, files = ".", extras = "-x '*.git*'")
    },
    contentType = "application/zip"
  )
  
  # Save State to file feature end
  # -----------------------------
  
  # Panel tab loaded second tab first to put it also to memory, set it back to first tab here
  updateTabsetPanel(session, "inputTablePanel", selected = "Conditions")
  
  #
  # Simulation
  #
  
  RunSim <- reactive({
    input$update
    isolate({
      withProgress({
        setProgress(message = "Calculating...")
        
        results <- SimThis(input, conditionTable, parameterTable)
        
        # Update the parameter fields with the new results
        if (input$fitting == "Local (fast)") {
          parametersToUpdate <- coef(results$fittingResults)
        } else { # input$fitting == "Global"
          parametersToUpdate <- results$fittingResults$optim$bestmem
        }
        
        for (name in names(parametersToUpdate)) {
          UpdateDTTable(name, newValue = parametersToUpdate[[name]])
        }
        
        results
      })
    })
  })
  
  RunSimVar <- reactive({
    input$update
    isolate({
      withProgress({
        setProgress(message = "Calculating parameter variability...")
        
        results <- SimThisVar(input, conditionTable, parameterTable)
        
        results
      })
    })
  })
  
  MatrixSim <- reactive({
    input$update
    isolate({
      if (input$update != 0) {
        withProgress({
          setProgress(value = 0, message = "Calculating waterfall plot...")
          validate(
            need(
              input$xAxisItem != input$yAxisItem,
              "Same variables selected for both X- and Y-axis, please change one of them."
            )
          )
          MatrixSimResults <- MultipleSim(input, conditionTable, parameterTable)
        })
      }
    })
  })
  
  #
  # Output
  #
  
  output$resultTable <- renderDataTable({
    result <- RunSim()$numericalSolution
    result <- result[ , !names(result) %in% c("Ab_C1_f", "Ab_C2_f", "EC50", "Ag_cell_t", "ADC_ex_f_E_ADC")]
    custom_formatted_results <- FormatSimRawDataForOutput(result)
    custom_formatted_results
  })
  
  output$AUC <- renderPrint({
    data <- RunSim()$numericalSolution
    AUCdegADC <- auc(data$time, data$degADC)[1]
    AUCdegADCperTime <- auc(data$time, data$degADCperTime)[1]
    str1 <- paste("AUC of sum of degradated ADC in tumor =", AUCdegADC, "nmol*h")
    str2 <- paste("AUC of degradated ADC in tumor =", AUCdegADCperTime, "nmol*h")
    cat(paste(str1, str2, sep = "\n"))
  })
  
  output$plot1 <- renderPlotly({
    show <- c("legendonly", rep(T, 5), rep("legendonly", 16))
    MakePlot1(RunSim()$numericalSolution, show, "Plasma PK")
  })
  
  output$plot2 <- renderPlotly({
    show <- c(rep("legendonly", 9), T, T, T, T, T, rep("legendonly", 8))
    MakePlot1(RunSim()$numericalSolution, show, "Number of molecules in a single tumor cell")
  })
  
  output$plot3 <- renderPlotly({
    show <- c(rep("legendonly", 18), T, rep("legendonly", 3))
    MakePlot1(RunSim()$numericalSolution, show, "Occupancy of drug target")
  })
  
  output$plot4 <- renderPlotly({
    show <- c(T, rep("legendonly", 13), rep(T, 4), rep("legendonly", 4))
    MakePlot1(RunSim()$numericalSolution, show, "Tumor volume")
  })
  
  output$plotExpData <- renderPlotly({
    ExpPlot(RunSim()$numericalSolution, input, "Comparison of simulation with experimental data")
  })
  
  output$plotVar<- renderPlotly({
    VarPlot(RunSimVar(), input, "Parameter variability")
  })
  
  output$fittingResult <- renderPrint({
    if (input$fitting == "Local (fast)") {
      PrintSummaryIfFittingExists(RunSim()$fittingResults)
    } else { # input$fitting == "Global"
      if (!is.null(RunSim()$fittingResults$optim$bestmem)) {
        data <- paste(input$fitTheseParameters, RunSim()$fittingResults$optim$bestmem)
        cat("Estimated parameters:", data)
      } 
    }
  })
  
  output$plot3d <- renderPlotly({
    validate(need(input$update != 0, "Press 'Calculate and plot' button first and progress bar should appear."))
    validate(need(MatrixSim(), "Wait for the calculations to happen."))
    
    isolate({
      withProgress({
        setProgress(message = "Building waterfall plot ...")
        MakeWaterfallPlot(input, MatrixSim())
      })
    })
  })
  
  output$figureCaption <- renderText({
    input$update
    isolate({
      if (input$zAxisItem == "Tumor size (waterfall plot)") {
        validate(need(MatrixSim(), ""))
        initialV_tumor <- sum(conditionTable[conditionTable$Initial %in% c("V_tumor_pro_mm3", "V_tumor_dyi_1_mm3", "V_tumor_dyi_2_mm3", "V_tumor_dyi_3_mm3"), "Value"])
        paste(
          "Figure 1. Waterfall plot that shows change in tumor volume (%) from initial", initialV_tumor,
          "cubic millimeters after", input$maxTime, "days for different parameter combinations. Parameter X is",
          input$xAxisItem, " (blue hue) and Y is", input$yAxisItem, "(asterisks)."
        )
      } else if (input$zAxisItem == "Sum of longest diameter (waterfall plot)") {
        validate(need(MatrixSim(), ""))
        initialV_tumor <- sum(conditionTable[conditionTable$Initial %in% c("V_tumor_pro_mm3", "V_tumor_dyi_1_mm3", "V_tumor_dyi_2_mm3", "V_tumor_dyi_3_mm3"), "Value"])
        initialD <- (sign(initialV_tumor)*initialV_tumor * 6 / pi)^(1/3)
        paste(
          "Figure 1. Waterfall plot that shows change in sum of longest diameter (%) from initial", round(initialD, digits = 2),
          "mm after", input$maxTime, "days for different parameter combinations. Parameter X is",
          input$xAxisItem, " (blue hue) and Y is", input$yAxisItem, "(asterisks)."
        )
      }
    })
  })
  
  output$resultTable3D <- renderDataTable({
    validate(need(MatrixSim(), "Wait for the calculations to happen."))
    result <- MatrixSim()
    custom_formatted_results <- FormatSimRawDataForOutput3D(input, result)
    custom_formatted_results
  })
}
