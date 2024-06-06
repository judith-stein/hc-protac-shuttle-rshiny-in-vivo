# Web GUI to plot a system of differential equations
# How to run:
# > runApp("app-dir/")

library(shiny) # Web GUI
library(plotly) # Interactive plots
library(deSolve) # For solving differential equations
library(DT) # Interactive DataTables in the UI
library(shinydashboard) # For UI layout with side and top bars
library(tibble)
library(tidyxl) # ALSO NEEDS package "tibble"! For excel file opening
library(dplyr) # For dataframe manipulation and %>%
library(V8)
library(shinyjs) # ALSO NEEDS package "V8"! For UI style
library(flux) # Calculation of AUC
library(DEoptim) # for global parameter fitting using DEoptim
library(minpack.lm) # For least squares fit using levenberg-marquart algorithm
library(ggplot2)
library(rxode2)

source("Model.R")
source("functions.R")

# Make it possible to download source files
addResourcePath(prefix = "src", directoryPath = ".")

# Location of the default data
databookFileName <- "www/simADC_databook_in_vivo_logistic.xlsx"

# Load default experimental data
defaultExpDataAsString <- readChar(
  "default/defaultData.tsv",
  file.info("default/defaultData.tsv")$size
)


SimThis <- function(input, conditions, parameters,
                    paramXname = NULL, paramXvalue = NULL,
                    paramYname = NULL, paramYvalue = NULL) {
  
  experimentalData <- LoadExpData(input$expDataAsString)
  
  initialConditions <- PickValuesToNumeric(conditions)
  initialParameters <- PickValuesToNumeric(parameters)
  
  initialParameters["thresholdV_tumor"] <- input$thresholdTumorVolume
  initialParameters["enableTmdd"] <- input$enableTmdd
  
  initialConditions["degADC"] <- 0
  initialConditions["A"] <- 0
  initialConditions["D"] <- 0
  initialConditions["Ab_C1_f"] <- initialConditions["Ab_C1_f"] * 0.001 / 150000 * 1e9 # convert from mg/kg to nmol/kg
  
  # Parameter variability
  initialParameters["eta.EC50"] <- 0
  initialParameters["eta.Ag_cell_t"] <- 0
  
  
  # Switch waterfall plot variables in-place
  if (!is.null(paramXname)) {
    # For possible x-values:
    initialParameters <- ReplaceIfThere(initialParameters, paramXname, paramXvalue)
    initialConditions <- ReplaceIfThere(initialConditions, paramXname, paramXvalue)
  }
  if (!is.null(paramYname)) {
  # For possible y-values:
  initialParameters <- ReplaceIfThere(initialParameters, paramYname, paramYvalue)
  initialConditions <- ReplaceIfThere(initialConditions, paramYname, paramYvalue)
  }
  
  
  # EventTable as Input for RxODE: dosing and observation (sampling) events
  et <- eventTable(amount.units = "1", time.units = "hours")
  timeSteps <- seq(from = 0, to = input$maxTime * 24, by = 1/2) 
  et$add.sampling(timeSteps)
  et$add.dosing(dosing.to = "ADC_C1_f",
                dose = initialConditions["Dose"] * 0.001 / 150000 * 1e9,
                nbr.doses = input$doseMaxT + 1,
                dosing.interval=input$doseT*24,
                start.time = input$startT)
  # Dosing for DAR
  et$add.dosing(dosing.to = "D",
                dose = initialConditions["Dose"] * 0.001 / 150000 * 1e9 * initialConditions["DAR"],
                nbr.doses = input$doseMaxT + 1,
                dosing.interval=input$doseT*24,
                start.time = input$startT)
  et$add.dosing(dosing.to = "A",
                dose = initialConditions["Dose"] * 0.001 / 150000 * 1e9,
                nbr.doses = input$doseMaxT + 1,
                dosing.interval=input$doseT*24,
                start.time = input$startT) 
  
  #################################### Fitting
  
  fitval <- NULL # To have it available for output even without fitting
  if (!is.null(input$fitTheseParameters)) {
    # Fit if toggled
    
    fitTheseParameters <- ParametersToCharacter(
      input$fitTheseParameters,
      initialParameters
    )
    
    fitTheseParameters <- sapply(fitTheseParameters, as.numeric)
    
    if (input$fitting == "Local (fast)") {
      fitval <- nls.lm(
        par = fitTheseParameters, 
        lower = rep(10^(-9), length(fitTheseParameters)), 
        fn = TakeVectorOfResiduals,
        fitTheseParameters = fitTheseParameters,
        fitting = input$fitting,
        currentExperimentalData = experimentalData,
        initialConditions = initialConditions, 
        initialParameters = initialParameters, 
        dataToFit = input$fitToTheseColumns,
        timeSteps = timeSteps, et = et,
        control = nls.lm.control(epsfcn = 1e-04)
      )
      
      estimatedParameters <- coef(fitval)
      
    } else { # input$fitting == "Global"
      # Use a genetic algorithm to fit the parameters
      lower <- vector(mode="integer", length=length(fitTheseParameters))
      lower[] <- input$lower
      upper <- vector(mode="integer", length=length(fitTheseParameters))
      upper[] <- input$upper
      fitval <- DEoptim(TakeVectorOfResiduals,
                        fitTheseParameters,
                        fitting = input$fitting,
                        currentExperimentalData = experimentalData,
                        initialConditions = initialConditions, 
                        initialParameters = initialParameters, 
                        dataToFit = input$fitToTheseColumns,
                        timeSteps = timeSteps, et = et,
                        lower = lower, upper = upper, 
                        control = list(NP = 50, VTR = input$VTR, itermax = input$itermax, trace = FALSE))
      
      names(fitval$optim$bestmem) <- names(fitTheseParameters)
      estimatedParameters <- fitval$optim$bestmem
    }
    
    
    initialParameters <- replace(
      initialParameters,
      names(estimatedParameters),
      estimatedParameters
    )
  } else {
    # initialParameters <- initialParameters
  }
  
  #################################### Fitting end
  
  # Solve equation with given parameters with RxODE
  numericalSolution <- as.data.frame(Model$solve(initialParameters,
                                                 et,
                                                 initialConditions)
  )
  
  # Convert ADC_C1_f and ADC_C2_f from nmol/Kg to nM for plotting: (nmol/kg) / (L/kg) = nM
  numericalSolution["ADC_C1_f_nM"] <- sapply(
    numericalSolution["ADC_C1_f"],
    function(x) x / initialParameters["V_C1_ADC"]
  )
  numericalSolution["ADC_C2_f_nM"] <- sapply(
    numericalSolution["ADC_C2_f"],
    function(x) x / initialParameters["V_C2_ADC"]
  )
  # Convert Ab_C1_f and Ab_C2_f from nmol/Kg to nM for plotting: (nmol/kg) / (L/kg) = nM
  numericalSolution["Ab_C1_f_nM"] <- sapply(
    numericalSolution["Ab_C1_f"],
    function(x) x / initialParameters["V_C1_Ab"]
  )
  numericalSolution["Ab_C2_f_nM"] <- sapply(
    numericalSolution["Ab_C2_f"],
    function(x) x / initialParameters["V_C2_Ab"]
  )
  
  # Intracellular occupancy of drug target
  SF <- 10^9 / (6.023 * 10^23)
  Drug_Target_cell_cyto_t_num <- initialParameters[["Drug_Target_cell_cyto_t"]] * initialParameters[["V_cell"]] / SF
  numericalSolution["Occupancy"] <- numericalSolution["Drug_cell_cyto_b_dt"] / Drug_Target_cell_cyto_t_num * 100
  
  numericalSolution["degADCperTime"] <- c()
  numericalSolution$degADCperTime[1] <- numericalSolution$degADC[1]
  if (length(numericalSolution$degADC) > 1) {
    for (i in 2:length(numericalSolution$degADC)) {
      numericalSolution$degADCperTime[i] <- numericalSolution$degADC[i] - numericalSolution$degADC[i-1]}
  }
  numericalSolution["ADC_in_tumor"] <- numericalSolution$ADC_ex_f / initialParameters[["E_ADC"]] * numericalSolution$V_tumor_mm3 * 10^-6 +                      
    numericalSolution$ADC_cell_b_ag * initialParameters[["NCL_tumor"]] * numericalSolution$V_tumor_mm3 * 10^-6 * SF +
    numericalSolution$ADC_cell_lyso_f * initialParameters[["NCL_tumor"]] * numericalSolution$V_tumor_mm3 * 10^-6 * SF +
    numericalSolution["degADCperTime"]
  
  numericalSolution["ADC_ex_f_E_ADC"] <- numericalSolution$ADC_ex_f / initialParameters[["E_ADC"]]
  
  # waterfall plot
  dataAtGoal <- numericalSolution[numericalSolution$time == input$maxTime * 24, ]
  comparableV_tumor <- as.numeric(dataAtGoal["V_tumor_pro_mm3"] +
                                    dataAtGoal["V_tumor_dyi_1_mm3"] +
                                    dataAtGoal["V_tumor_dyi_2_mm3"] +
                                    dataAtGoal["V_tumor_dyi_3_mm3"])
  comparableD <- (sign(comparableV_tumor)*comparableV_tumor * 6 / pi)^(1/3)
  
  originalV_tumor <- as.numeric(initialConditions["V_tumor_pro_mm3"] +
                                  initialConditions["V_tumor_dyi_1_mm3"] +
                                  initialConditions["V_tumor_dyi_2_mm3"] +
                                  initialConditions["V_tumor_dyi_3_mm3"])
  originalD <- (sign(originalV_tumor)*originalV_tumor * 6 / pi)^(1/3)
  
  percentageChangeInTumorVolume <- (comparableV_tumor - originalV_tumor) / originalV_tumor * 100
  percentageChangeInDiameter <- (comparableD - originalD) / originalD * 100

  
  return(list(numericalSolution = numericalSolution, fittingResults = fitval,
              V_tumor = percentageChangeInTumorVolume, diameter = percentageChangeInDiameter))
}


SimThisVar <- function(input, conditions, parameters) {
  
  initialConditions <- PickValuesToNumeric(conditions)
  initialParameters <- PickValuesToNumeric(parameters)
  
  initialParameters["thresholdV_tumor"] <- input$thresholdTumorVolume
  initialParameters["enableTmdd"] <- input$enableTmdd
  
  initialConditions["A"] <- 0
  initialConditions["D"] <- 0
  initialConditions["Ab_C1_f"] <- initialConditions["Ab_C1_f"] * 0.001 / 150000 * 1e9 # convert from mg/kg to nmol/kg
  
  # Parameter variability
  std.EC50 <- input$std.EC50
  std.Ag_cell_t <- input$std.Ag_cell_t
  omega <- lotri(eta.EC50 ~ std.EC50^2,
                 eta.Ag_cell_t ~ std.Ag_cell_t^2
  )
  
  # EventTable as Input for RxODE: dosing and observation (sampling) events
  et <- eventTable(amount.units = "1", time.units = "hours")
  timeSteps <- seq(from = 0, to = input$maxTime * 24, by = 1/2) 
  et$add.sampling(timeSteps)
  et$add.dosing(dosing.to = "ADC_C1_f",
                dose = initialConditions["Dose"] * 0.001 / 150000 * 1e9,
                nbr.doses = input$doseMaxT + 1,
                dosing.interval=input$doseT*24,
                start.time = input$startT)
  # Dosing for DAR
  et$add.dosing(dosing.to = "D",
                dose = initialConditions["Dose"] * 0.001 / 150000 * 1e9 * initialConditions["DAR"],
                nbr.doses = input$doseMaxT + 1,
                dosing.interval=input$doseT*24,
                start.time = input$startT)
  et$add.dosing(dosing.to = "A",
                dose = initialConditions["Dose"] * 0.001 / 150000 * 1e9,
                nbr.doses = input$doseMaxT + 1,
                dosing.interval=input$doseT*24,
                start.time = input$startT) 
  
  # Solve equation with given parameters with RxODE
  numericalSolution <- rxSolve(Model,
                               initialParameters,
                               et,
                               initialConditions,
                               omega = omega,
                               nSub = input$nSub
  )
  
  return(numericalSolution = numericalSolution)
}

MultipleSim <- function(input, conditionTable, parameterTable) {
  xn <- ReadStepsFromJson(input$xAxisInstructions)
  yn <- ReadStepsFromJson(input$yAxisInstructions)
  
  x <- ReadRangeFromJson(input$xAxisInstructions)
  y <- ReadRangeFromJson(input$yAxisInstructions)
  
  # Initialize result matrix that will be z-axis
  # Prepare for three surfaces, even though we use only one at first
  matrixTumorSize <- matrix(data = NA, nrow = yn, ncol = xn)
  matrixDiameter <- matrix(data = NA, nrow = yn, ncol = xn)
  
  # Initialize messages for progress bar shown when calculations run
  previousResult <- ""
  progressCounter <- 0
  maxProg <- 1
  stepProg <- ""
  
  for (i in 1:yn) {
    for (j in 1:xn) {
      
      if (input$zAxisItem == "Tumor size (waterfall plot)") {
        # Tumor volume after x days with these variables
        
        # Update progress bar in the UI to show how many data points to be calculated
        progressCounter <- progressCounter + 1
        nowCalculating <- paste0(
          "Now calculating V_tumor % with x: ", x[j], " and y: ", y[i],
          " ... (", progressCounter, stepProg, "/", xn * yn, ")"
        )
        incProgress(amount = 1 / (xn * yn * maxProg + 1), message = paste(previousResult, nowCalculating))
        
        tv_after_x <- SimThis(input, conditionTable, parameterTable,
                              paramXname = input$xAxisItem,
                              paramXvalue = x[j],
                              paramYname = input$yAxisItem,
                              paramYvalue = y[i]
        )$V_tumor
        if (!is.null(tv_after_x)) {
          matrixTumorSize[i, j] <- tv_after_x
        } else {
          matrixTumorSize[i, j] <- NA
        }
        
        previousResult <- paste("Previously, V_tumor % was", round(tv_after_x), "(when x:", x[j], "and y:", y[i], ").")
      }
      
      if (input$zAxisItem == "Sum of longest diameter (waterfall plot)") {
        # Longest diameter after x days with these variables
        
        # Update progress bar in the UI to show how many data points to be calculated
        progressCounter <- progressCounter + 1
        nowCalculating <- paste0(
          "Now calculating diameter % with x: ", x[j], " and y: ", y[i],
          " ... (", progressCounter, stepProg, "/", xn * yn, ")"
        )
        incProgress(amount = 1 / (xn * yn * maxProg + 1), message = paste(previousResult, nowCalculating))
        
        d_after_x <- SimThis(input, conditionTable, parameterTable,
                             paramXname = input$xAxisItem,
                             paramXvalue = x[j],
                             paramYname = input$yAxisItem,
                             paramYvalue = y[i]
        )$diameter
        if (!is.null(d_after_x)) {
          matrixDiameter[i, j] <- d_after_x
        } else {
          matrixDiameter[i, j] <- NA
        }
        
        previousResult <- paste("Previously, diameter % was", round(d_after_x), "(when x:", x[j], "and y:", y[i], ").")
      }
      
    }
  }
  
  return(list(
    "x" = x, "y" = y,
    "Tumor size (waterfall plot)" = matrixTumorSize,
    "Sum of longest diameter (waterfall plot)" = matrixDiameter
  ))
}



MakePlot1 <- function(simData, vis, title) {
  mfigure <- plot_ly(simData, x = ~ time / 24, y = ~V_tumor_mm3, name = "V_tumor (mm^3)", type = "scatter", mode = "lines", visible = vis[1]) %>%
    add_trace(y = ~ADC_C1_f_nM, name = "ADC_C1_f (nM)", mode = "lines", visible = vis[2]) %>%
    add_trace(y = ~ADC_C2_f_nM, name = "ADC_C2_f (nM)", mode = "lines", visible = vis[3]) %>%
    add_trace(y = ~Drug_C1_f, name = "Drug_C1_f (nM)", mode = "lines", visible = vis[4]) %>%
    add_trace(y = ~Drug_C2_f, name = "Drug_C2_f (nM)", mode = "lines", visible = vis[5]) %>%
    add_trace(y = ~Drug_C1_b_ntp, name = "Drug_C1_b_ntp (nM)", mode = "lines", visible = vis[6]) %>%
    add_trace(y = ~DAR, name = "DAR", mode = "lines", visible = vis[7]) %>%
    add_trace(y = ~ADC_ex_f, name = "ADC_ex_f (nM)", mode = "lines", visible = vis[8]) %>%
    add_trace(y = ~ADC_ex_f_E_ADC, name = "ADC_ex_f/E_ADC (nM)", mode = "lines", visible = "legendonly") %>%
    add_trace(y = ~Drug_ex_f, name = "Drug_ex_f (nmol)", mode = "lines", visible = vis[9]) %>%
    add_trace(y = ~ADC_cell_b_ag, name = "ADC_cell_b_ag", mode = "lines", visible = vis[10]) %>%
    add_trace(y = ~ADC_cell_lyso_f, name = "ADC_cell_lyso_f", mode = "lines", visible = vis[11]) %>%
    add_trace(y = ~Drug_cell_lyso_f, name = "Drug_cell_lyso_f", mode = "lines", visible = vis[12]) %>%
    add_trace(y = ~Drug_cell_cyto_f, name = "Drug_cell_cyto_f", mode = "lines", visible = vis[13]) %>%
    add_trace(y = ~Drug_cell_cyto_b_dt, name = "Drug_cell_cyto_b_dt", mode = "lines", visible = vis[14]) %>%
    add_trace(y = ~ADC_TMDD_cell_b_ag, name = "ADC_TMDD_cell_b_ag", mode = "lines", visible = "legendonly") %>%
    add_trace(y = ~ADC_TMDD_cell_lyso_f, name = "ADC_TMDD_cell_lyso_f", mode = "lines", visible = "legendonly") %>%
    add_trace(y = ~Drug_TMDD_cell_lyso_f, name = "Drug_TMDD_cell_lyso_f", mode = "lines", visible = "legendonly") %>%
    add_trace(y = ~Drug_TMDD_cell_cyto_f, name = "Drug_TMDD_cell_cyto_f", mode = "lines", visible = "legendonly") %>%
    add_trace(y = ~V_tumor_pro_mm3, name = "V_tumor_pro_mm3", mode = "lines", visible = vis[15]) %>%
    add_trace(y = ~V_tumor_dyi_1_mm3, name = "V_tumor_dyi_1_mm3", mode = "lines", visible = vis[16]) %>%
    add_trace(y = ~V_tumor_dyi_2_mm3, name = "V_tumor_dyi_2_mm3", mode = "lines", visible = vis[17]) %>%
    add_trace(y = ~V_tumor_dyi_3_mm3, name = "V_tumor_dyi_3_mm3", mode = "lines", visible = vis[18]) %>%
    add_trace(y = ~Occupancy, name = "Target_occup.", mode = "lines", visible = vis[19]) %>%
    add_trace(y = ~degADC, name = "Sum of deg. ADC (nmol)", mode = "lines", visible = vis[20]) %>%
    add_trace(y = ~degADCperTime, name = "Deg. ADC (nmol)", mode = "lines", visible = vis[21]) %>%
    add_trace(y = ~ADC_in_tumor, name = "ADC in tumor (nmol)", mode = "lines", visible = vis[22]) %>%
    add_trace(y = ~Ab_C1_f, name = "Ab_C1_f (nmol/kg)", mode = "lines", visible = "legendonly") %>%
    add_trace(y = ~Ab_C2_f, name = "Ab_C2_f (nmol/kg)", mode = "lines", visible = "legendonly") %>%
    add_trace(y = ~Ab_C1_f_nM, name = "Ab_C1_f (nM)", mode = "lines", visible = "legendonly") %>%
    add_trace(y = ~Ab_C2_f_nM, name = "Ab_C2_f (nM)", mode = "lines", visible = "legendonly") %>%
    add_trace(y = ~Ab_ex_f, name = "Ab_ex_f (nM)", mode = "lines", visible = "legendonly") %>%
    add_trace(y = ~Ab_cell_b_ag, name = "Ab_cell_b_ag", mode = "lines", visible = "legendonly") %>%
    add_trace(y = ~Ab_TMDD_cell_b_ag, name = "Ab_TMDD_cell_b_ag", mode = "lines", visible = "legendonly") %>%
    add_trace(y = ~Ag_C1_f, name = "Ag_C1_f (nmol/kg)", mode = "lines", visible = "legendonly") %>%
    add_trace(y = ~Ag_C1_b, name = "Ag_C1_b (nmol/kg)", mode = "lines", visible = "legendonly") %>%
    add_trace(y = ~Ag_ex_f, name = "Ag_ex_f (nM)", mode = "lines", visible = "legendonly") %>%
    add_trace(y = ~Ag_ex_b, name = "Ag_ex_b (nM)", mode = "lines", visible = "legendonly") %>%
    layout(
      margin = list(t = 50), title = title,
      xaxis = list(title = "Time (days)"),
      yaxis = list(title = "Simulation")
    ) %>%
    config(displaylogo = FALSE, toImageButtonOptions = list(format = "png")) # To hide the link
  
  return(mfigure)
}

# Plot with experimental data
ExpPlot <- function(simData, input, title) {
  
  experimentalData <- LoadExpData(input$expDataAsString)
  
  mfigure <- ggplot() +
    geom_line(data=simData, aes(x=time/24, y=V_tumor_mm3, color="Sim")) +
    geom_point(data=experimentalData, aes(x=time, y=V_tumor_mm3, color="Exp")) +
    geom_ribbon(data=experimentalData, aes(x=time, y=V_tumor_mm3, ymin=V_tumor_mm3-dev,ymax=V_tumor_mm3+dev), fill="lightblue", alpha=0.5) +
    theme(legend.position="right") +
    scale_color_manual(name  = "", values = c("blue", "red"), labels = c("Sim", "Exp")) +
    ggtitle(title) +
    labs(x = "Time (days)", y = paste0("Tumor volume (mm", tags$sup("3"), ")")) 
  
  return(mfigure)
}

VarPlot <- function(simData, input, title) {
  experimentalData <- LoadExpData(input$expDataAsString2)

  level <- input$level
  ymin <- (1-level)/2
  ymax <- 1 - (1-level)/2
  parm <- "V_tumor_mm3"
  simLegendTitle <- "Sim"
  data <- confint(simData, parm, level)
  data$time <- data$time/24

  if ("eff" %in% names(data)) {
    data <- data %>%
      tidyr::pivot_wider(
        id_cols = time,
        names_from = p1,
        values_from = eff
      )
    p <- plot_ly(data)
    p <- p %>%
      add_trace(
        x=~time,
        y=~`0.5`,
        name=simLegendTitle,
        type = "scatter",
        mode = "lines"
      ) %>%
      add_ribbons(
        data=data,
        x=~time,
        ymin=as.formula(paste0("~`", ymin, "`")),
        ymax=as.formula(paste0("~`", ymax, "`")),
        name=paste("Confidence interval range:", ymin, "-", ymax),
        type = "scatter",
        mode = "lines"
      )
  } else if ("p50" %in% names(data)) {
    cols <- colnames(data)
    pcols <- sort(cols[startsWith(cols, "p") & !(cols == "p1")])
    data <- data %>%
      tidyr::pivot_longer(
          cols=pcols,
          names_to="p2",
          values_to="tmp"
        ) %>%
      tidyr::pivot_wider(
        id_cols=time,
        names_from=c("p2", "p1"),
        values_from="tmp"
      )
    p <- plot_ly(data)
    p <- p %>%
      add_trace(
        x=~time,
        y=~`p50_0.5`,
        name=paste(simLegendTitle, "0.5"),
        type = "scatter",
        mode = "lines"
      ) %>%
      add_ribbons(
        data=data,
        x=~time,
        ymin=as.formula(paste0("~`", pcols[1], "_0.5`")),
        ymax=as.formula(paste0("~`", pcols[length(pcols)], "_0.5`")),
        name=paste("0.5 CL:", ymin, ymax),
        type = "scatter",
        mode = "lines"
      ) %>%
      add_trace(
        x=~time,
        y=as.formula(paste0("~`p50_", ymin, "`")),
        name=paste(simLegendTitle, ymin),
        type = "scatter",
        mode = "lines"
      ) %>%
      add_ribbons(
        data=data,
        x=~time,
        ymin=as.formula(paste0("~`", pcols[1], "_", ymin, "`")),
        ymax=as.formula(paste0("~`", pcols[length(pcols)], "_", ymin, "`")),
        name=paste(ymin, "CL:", ymin, ymax),
        type = "scatter",
        mode = "lines"
      ) %>%
      add_trace(
        x=~time,
        y=as.formula(paste0("~`p50_", ymax, "`")),
        name=paste(simLegendTitle, ymax),
        type = "scatter",
        mode = "lines"
      ) %>%
      add_ribbons(
        data=data,
        x=~time,
        ymin=as.formula(paste0("~`", pcols[1], "_", ymax, "`")),
        ymax=as.formula(paste0("~`", pcols[length(pcols)], "_", ymax, "`")),
        name=paste(ymax, "CL:", ymin, ymax),
        type = "scatter",
        mode = "lines"
      )
  }

  if (!is.null(experimentalData)) {
    experimentalData$ymin <- experimentalData$V_tumor_mm3 - experimentalData$dev
    experimentalData$ymax <- experimentalData$V_tumor_mm3 + experimentalData$dev
    p <- p %>%
      add_trace(
        data=experimentalData,
        x=~time,
        y=~V_tumor_mm3,
        name="Exp",
        type = "scatter",
        mode = "lines+markers"
      ) %>%
      add_ribbons(
        data=experimentalData,
        x=~time,
        ymin=~ymin,
        ymax=~ymax,
        name=paste("Exp deviance"),
        type = "scatter",
        mode = "lines"
      )
  }

  p <- p %>%
    layout(
      margin = list(t = 50),
      title = title,
      xaxis = list(title = "Time (days)"),
      yaxis = list(title = paste0("Tumor volume (mm", tags$sup("3"), ")"))
    ) %>%
    config(
      displaylogo = FALSE,
      toImageButtonOptions = list(format = "png")
    )
  return(p)
}

MakeWaterfallPlot <- function(input, xyzData) {
  # Waterfall plot
  # Note, that NA's do not show on the plot at all
  
  x <- xyzData[["x"]]
  y <- xyzData[["y"]]
  d <- xyzData[[input$zAxisItem]]
  
  # Create y-markers because x gives color and we want to
  # distinquish y somehow as well
  leny <- length(y)
  ym <- rep("", leny)
  for (i in 1:leny) {
    # eg. for four values, make "*<br>" "**<br>" "***<br>" "***<br>* "
    # which would look like: * ** *** ***
    #                                 *
    k <- i
    layers <- 0
    while (k > 0) {
      ym[i] <- paste0(c(ym[i], rep("*", min(3, k)), "<br>"), collapse = "")
      layers <- layers + 1
      k <- k - 3
    }
    if (layers > 1) ym[i] <- paste0(ym[i], " ")
  }
  
  df <- data.frame(
    col = rep(x, each = nrow(d)),
    row = rep(y, ncol(d)),
    ym = rep(ym, ncol(d)),
    value = as.vector(d)
  )
  
  df2 <- df[order(-df$value), ]
  
  xInOrder <- paste0("X: ", df2$col, " and Y: ", df2$row)
  yInOrder <- df2$value
  
  xWithColors <- df2$col
  yWithStars <- df2$ym
  
  if (input$zAxisItem == "Tumor size (waterfall plot)") {
    title = list(title = "Tumor size change (%) from initial")
  } else {title = list(title = "Sum of longest diameter change (%) from initial")}
  
  plotWaterfall <- plot_ly(
    x = xInOrder,
    y = yInOrder,
    name = "Waterfall Plot",
    type = "bar",
    marker = list(
      color = xWithColors,
      colorscale = list( # From light blue:
        list(0, "rgb(198,219,239)"),
        # To dark blue:
        list(1, "rgb(33,113,181)")
      )
    )
  ) %>%
    layout(
      xaxis = list(
        categoryorder = "array",
        # Only way to order bars:
        categoryarray = yInOrder
      ),
      yaxis = title
    ) %>%
    add_annotations(
      text = yWithStars,
      x = xInOrder,
      y = yInOrder,
      xref = "x",
      yref = "y",
      showarrow = FALSE
    )
  return(plotWaterfall)
}


# Databook feature
tx <- xlsx_cells(databookFileName)
txFormat <- xlsx_formats(databookFileName)

GetInitialListFromVertical <- function(idSheet, labelStartRow, labelColumn, idStartRow, idColumn) {
  labels <- tx %>%
    filter(sheet == idSheet, row >= labelStartRow, col == labelColumn) %>%
    pull(character)
  
  idList <- tx %>%
    filter(sheet == idSheet, row >= idStartRow, col == idColumn) %>%
    pull(character)
  
  names(idList) <- labels
  idList <- c("", idList) # Add empty item to the front for initial empty selection
  return(idList)
}

GetInitialListFromHorizontal <- function(idSheet, labelRow, labelStartColumn, idRow, idStartColumn) {
  labels <- tx %>%
    filter(sheet == idSheet, row == labelRow, col >= labelStartColumn) %>%
    pull(character)
  
  idList <- tx %>%
    filter(sheet == idSheet, row == idRow, col >= idStartColumn) %>%
    pull(character)
  
  names(idList) <- labels
  idList <- c("", idList) # Add empty item to the front for initial empty selection
  return(idList)
}

abList <- GetInitialListFromVertical("Antibodies", 16, 1, 16, 2)
agList <- GetInitialListFromHorizontal("Cell Lines", 2, 7, 3, 7)
cellList <- GetInitialListFromVertical("Cell Lines", 4, 1, 4, 2)
linkList <- GetInitialListFromVertical("Linkers", 4, 1, 4, 2)
payloadList <- GetInitialListFromVertical("Payloads", 4, 1, 4, 2)

# Databook title
databookTitle <- tx %>%
  filter(sheet == "Antibodies", row == 2, col == 1) %>%
  pull(character)

# Read git version number
if (file.exists("www/version.txt")) {
  versionString <- read.table(file = "www/version.txt", header = F, nrows = 1)[[1, 1]]
  versionNumbers <- tail(strsplit(versionString, "-")[[1]], n = 2)
} else {
  versionString <- "unknown"
  versionNumbers <- c("unknown", "unknown")
}

# Shinyjs code to highlight if a field is updated
jsCode <- 'shinyjs.backgroundCol = function(params) {
                var defaultParams = {
                    id : null,
                    col : "red"
                };
                params = shinyjs.getParams(params, defaultParams);
                var el = $("#" + params.id);
                el.css("background-color", params.col);
            }'
