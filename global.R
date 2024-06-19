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
library(glue)
library(readr)

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


#
############### Simulation ----------------
#

SimThis <- function(input, conditions, parameters) {

  experimentalData <- LoadExpData(input$expDataAsString)
  
  initialConditions <- PickValuesToNumeric(conditions)
  initialParameters <- PickValuesToNumeric(parameters)
  initialParameters["max"] <- input$max
  
  # EventTable as Input for RxODE: dosing and observation (sampling) events
  et <- eventTable(amount.units = "1", time.units = "hours")
  timeSteps <- seq(from = 0, to = input$maxTime * 24, by = 1/2) 
  et$add.sampling(timeSteps)
  et$add.dosing(dosing.to = input$doseTo,
                dose = input$dose * 0.001 / initialParameters["MW_Ab"] * 1e9, # convert from mg/kg to nmol/kg
                nbr.doses = input$doseMaxT + 1,
                dosing.interval=input$doseT*24,
                start.time = input$startT)
  et$add.dosing(dosing.to = "Drug_C1_f",
                dose = input$doseDrug * 1e6 / initialParameters["MW_Drug"] / initialParameters["V_C1_Drug"], # From mg/kg to nmol/L
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

  # writeLines(BuildModelForN(max = input$max), "model.txt")
  Model <<- rxode2(BuildModelForN(max = input$max))
  
  # Solve equation with given parameters with RxODE
  numericalSolution <- as.data.frame(Model$solve(initialParameters,
                                                 et,
                                                 initialConditions)
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
  
  
  return(list(numericalSolution = numericalSolution, fittingResults = fitval))
}

#
#################### Plots ---------------------
#

# Plasma and rest
MakePlot1 <- function(simData, vis, title, max) {
  mfigure <- plot_ly(simData, x = ~ time / 24, y = ~V_tumor_mm3, name = "V_tumor (mm^3)", type = "scatter", mode = "lines", visible = vis[1]) %>%
    add_trace(y = ~Ab_C1_f_nM, name = "Ab_C1_f (nM)", mode = "lines", visible = vis[2]) %>%
  add_trace(y = ~Ab_C1_f, name = "Ab_C1_f (nmol/kg)", mode = "lines", visible = "legendonly") 
  for (i in 1:max) {
    mfigure <- mfigure %>%
      add_trace(y = as.formula(paste0("~Ab_C1_b", i)), name = paste0("Ab_C1_b", i, " (nmol/kg)"), mode = "lines", visible = "legendonly") %>%
      add_trace(y = as.formula(paste0("~Ab_C1_m", i)), name = paste0("Ab_C1_m", i, " (nmol/kg)"), mode = "lines", visible = "legendonly") 
    if (max > 1 && i < max) {
      for (j in 1:(max-i)) {
        mfigure <- mfigure %>%
          add_trace(y = as.formula(paste0("~Ab_C1_b", i, "_m", j)), name = paste0("Ab_C1_b", i, "_m", j, " (nmol/kg)"), mode = "lines", visible = "legendonly") 
      }
    }
  }
  mfigure <- mfigure %>%
    add_trace(y = ~Ab_C2_f_nM, name = "Ab_C2_f (nM)", mode = "lines", visible = vis[3]) 
  for (i in 1:max) {
    mfigure <- mfigure %>%
      add_trace(y = as.formula(paste0("~Ab_C2_b", i)), name = paste0("Ab_C2_b", i, " (nmol/kg)"), mode = "lines", visible = "legendonly") %>%
      add_trace(y = as.formula(paste0("~Ab_C2_m", i)), name = paste0("Ab_C2_m", i, " (nmol/kg)"), mode = "lines", visible = "legendonly") 
    if (max > 1 && i < max) {
      for (j in 1:(max-i)) {
        mfigure <- mfigure %>%
          add_trace(y = as.formula(paste0("~Ab_C2_b", i, "_m", j)), name = paste0("Ab_C2_b", i, "_m", j, " (nmol/kg)"), mode = "lines", visible = "legendonly") 
      }
    }
  }
  mfigure <- mfigure %>%
    add_trace(y = ~Drug_C1_f, name = "Drug_C1_f (nM)", mode = "lines", visible = vis[4]) %>%
    add_trace(y = ~Meta1_C1_f, name = "Meta1_C1_f (nM)", mode = "lines", visible = vis[5]) %>%
    add_trace(y = ~Meta2_C1_f, name = "Meta2_C1_f (nM)", mode = "lines", visible = vis[6]) %>%
    add_trace(y = ~Drug_C2_f, name = "Drug_C2_f (nM)", mode = "lines", visible = vis[7]) %>%
    add_trace(y = ~Drug_C1_b_ntp, name = "Drug_C1_b_ntp (nM)", mode = "lines", visible = vis[8]) %>%
    add_trace(y = ~Ab_ex_f, name = "Ab_ex_f (nM)", mode = "lines", visible = vis[9]) 
  for (i in 1:max) {
    mfigure <- mfigure %>%
      add_trace(y = as.formula(paste0("~Ab_ex_b", i)), name = paste0("Ab_ex_b", i, " (nM)"), mode = "lines", visible = "legendonly") %>%
      add_trace(y = as.formula(paste0("~Ab_ex_m", i)), name = paste0("Ab_ex_m", i, " (nM)"), mode = "lines", visible = "legendonly") 
    if (max > 1 && i < max) {
      for (j in 1:(max-i)) {
        mfigure <- mfigure %>%
          add_trace(y = as.formula(paste0("~Ab_ex_b", i, "_m", j)), name = paste0("Ab_ex_b", i, "_m", j, " (nM)"), mode = "lines", visible = "legendonly") 
      }
    }
  }
  mfigure <- mfigure %>%
    add_trace(y = ~Drug_ex_f, name = "Drug_ex_f (nmol)", mode = "lines", visible = vis[10]) %>%
    add_trace(y = ~Meta1_ex_f, name = "Meta1_ex_f (nmol)", mode = "lines", visible = vis[11]) %>%
    add_trace(y = ~Meta2_ex_f, name = "Meta2_ex_f (nmol)", mode = "lines", visible = vis[12]) %>%
    add_trace(y = ~Ab_cell_f_b_ag, name = "Ab_cell_f_b_ag", mode = "lines", visible = vis[13]) 
  for (i in 1:max) {
    mfigure <- mfigure %>%
      add_trace(y = as.formula(paste0("~Ab_cell_b", i, "_b_ag")), name = paste0("Ab_cell_b", i, "_b_ag"), mode = "lines", visible = "legendonly") %>%
      add_trace(y = as.formula(paste0("~Ab_cell_m", i, "_b_ag")), name = paste0("Ab_cell_m", i, "_b_ag"), mode = "lines", visible = "legendonly") 
    if (max > 1 && i < max) {
      for (j in 1:(max-i)) {
        mfigure <- mfigure %>%
          add_trace(y = as.formula(paste0("~Ab_cell_b", i, "_m", j, "_b_ag")), name = paste0("Ab_cell_b", i, "_m", j, "_b_ag"), mode = "lines", visible = "legendonly") 
      }
    }
  }
  for (i in 1:max) {
    mfigure <- mfigure %>%
      add_trace(y = as.formula(paste0("~Ab_cell_lyso_b", i)), name = paste0("Ab_cell_lyso_b", i), mode = "lines", visible = "legendonly") %>%
      add_trace(y = as.formula(paste0("~Ab_cell_lyso_m", i)), name = paste0("Ab_cell_lyso_m", i), mode = "lines", visible = "legendonly") 
    if (max > 1 && i < max) {
      for (j in 1:(max-i)) {
        mfigure <- mfigure %>%
          add_trace(y = as.formula(paste0("~Ab_cell_lyso_b", i, "_m", j)), name = paste0("Ab_cell_lyso_b", i, "_m", j), mode = "lines", visible = "legendonly") 
      }
    }
  }
  mfigure <- mfigure %>%
    add_trace(y = ~Drug_cell_lyso_f, name = "Drug_cell_lyso_f", mode = "lines", visible = vis[14]) %>%
    add_trace(y = ~Meta1_cell_lyso_f, name = "Meta1_cell_lyso_f", mode = "lines", visible = vis[15]) %>%
    add_trace(y = ~Meta2_cell_lyso_f, name = "Meta2_cell_lyso_f", mode = "lines", visible = vis[16]) %>%
    add_trace(y = ~Drug_cell_cyto_f, name = "Drug_cell_cyto_f", mode = "lines", visible = vis[17]) %>%
    add_trace(y = ~Meta1_cell_cyto_f, name = "Meta1_cell_cyto_f", mode = "lines", visible = vis[18]) %>%
    add_trace(y = ~Meta2_cell_cyto_f, name = "Meta2_cell_cyto_f", mode = "lines", visible = vis[10]) %>%
    add_trace(y = ~Drug_cell_cyto_b_dt, name = "Drug_cell_cyto_b_dt", mode = "lines", visible = vis[20]) %>%
    add_trace(y = ~Meta2_cell_cyto_b_dt, name = "Meta2_cell_cyto_b_dt", mode = "lines", visible = vis[21]) %>%
    add_trace(y = ~V_tumor_pro_mm3, name = "V_tumor_pro_mm3", mode = "lines", visible = vis[22]) %>%
    add_trace(y = ~V_tumor_dyi_1_mm3, name = "V_tumor_dyi_1_mm3", mode = "lines", visible = vis[23]) %>%
    add_trace(y = ~V_tumor_dyi_2_mm3, name = "V_tumor_dyi_2_mm3", mode = "lines", visible = vis[24]) %>%
    add_trace(y = ~V_tumor_dyi_3_mm3, name = "V_tumor_dyi_3_mm3", mode = "lines", visible = vis[25]) %>%
    layout(
      margin = list(t = 50), title = title,
      xaxis = list(title = "Time (days)"),
      yaxis = list(title = "Simulation")
    ) %>%
    config(displaylogo = FALSE, toImageButtonOptions = list(format = "png")) # To hide the link
  
  return(mfigure)
}

# Tumor cell
MakePlot2 <- function(simData, title, max) {
  mfigure <- plot_ly(simData, x = ~ time / 24, y = ~Ab_cell_f_b_ag, name = "Ab_cell_f_b_ag", type = "scatter", mode = "lines") 
  for (i in 1:max) {
    mfigure <- mfigure %>%
      add_trace(y = as.formula(paste0("~Ab_cell_b", i, "_b_ag")), name = paste0("Ab_cell_b", i, "_b_ag"), mode = "lines", visible = "legendonly") %>%
      add_trace(y = as.formula(paste0("~Ab_cell_m", i, "_b_ag")), name = paste0("Ab_cell_m", i, "_b_ag"), mode = "lines", visible = "legendonly") 
    if (max > 1 && i < max) {
      for (j in 1:(max-i)) {
        mfigure <- mfigure %>%
          add_trace(y = as.formula(paste0("~Ab_cell_b", i, "_m", j, "_b_ag")), name = paste0("Ab_cell_b", i, "_m", j, "_b_ag"), mode = "lines", visible = "legendonly") 
      }
    }
  }
  for (i in 1:max) {
    mfigure <- mfigure %>%
      add_trace(y = as.formula(paste0("~Ab_cell_lyso_b", i)), name = paste0("Ab_cell_lyso_b", i), mode = "lines", visible = "legendonly") %>%
      add_trace(y = as.formula(paste0("~Ab_cell_lyso_m", i)), name = paste0("Ab_cell_lyso_m", i), mode = "lines", visible = "legendonly") 
    if (max > 1 && i < max) {
      for (j in 1:(max-i)) {
        mfigure <- mfigure %>%
          add_trace(y = as.formula(paste0("~Ab_cell_lyso_b", i, "_m", j)), name = paste0("Ab_cell_lyso_b", i, "_m", j), mode = "lines", visible = "legendonly") 
      }
    }
  }
  mfigure <- mfigure %>%
    add_trace(y = ~Drug_cell_lyso_f, name = "Drug_cell_lyso_f", mode = "lines") %>%
    add_trace(y = ~Meta1_cell_lyso_f, name = "Meta1_cell_lyso_f", mode = "lines") %>%
    add_trace(y = ~Meta2_cell_lyso_f, name = "Meta2_cell_lyso_f", mode = "lines") %>%
    add_trace(y = ~Drug_cell_cyto_f, name = "Drug_cell_cyto_f", mode = "lines") %>%
    add_trace(y = ~Meta1_cell_cyto_f, name = "Meta1_cell_cyto_f", mode = "lines") %>%
    add_trace(y = ~Meta2_cell_cyto_f, name = "Meta2_cell_cyto_f", mode = "lines") %>%
    add_trace(y = ~Drug_cell_cyto_b_dt, name = "Drug_cell_cyto_b_dt", mode = "lines") %>%
    add_trace(y = ~Meta2_cell_cyto_b_dt, name = "Meta2_cell_cyto_b_dt", mode = "lines") %>%
    layout(
      margin = list(t = 50), title = title,
      xaxis = list(title = "Time (days)"),
      yaxis = list(title = "Simulation")
    ) %>%
    config(displaylogo = FALSE, toImageButtonOptions = list(format = "png")) # To hide the link
  
  return(mfigure)
}

# Tumor volume
MakePlot3 <- function(simData, vis, title, max) {
  mfigure <- plot_ly(simData, x = ~ time / 24, y = ~V_tumor_mm3, name = "V_tumor (mm^3)", type = "scatter", mode = "lines", visible = vis[1]) %>%
    add_trace(y = ~Ab_ex_f, name = "Ab_ex_f (nM)", mode = "lines", visible = vis[2]) 
  for (i in 1:max) {
    mfigure <- mfigure %>%
      add_trace(y = as.formula(paste0("~Ab_ex_b", i)), name = paste0("Ab_ex_b", i, " (nM)"), mode = "lines", visible = "legendonly") %>%
      add_trace(y = as.formula(paste0("~Ab_ex_m", i)), name = paste0("Ab_ex_m", i, " (nM)"), mode = "lines", visible = "legendonly") 
    if (max > 1 && i < max) {
      for (j in 1:(max-i)) {
        mfigure <- mfigure %>%
          add_trace(y = as.formula(paste0("~Ab_ex_b", i, "_m", j)), name = paste0("Ab_ex_b", i, "_m", j, " (nM)"), mode = "lines", visible = "legendonly") 
      }
    }
  }
  mfigure <- mfigure %>%
    add_trace(y = ~Drug_ex_f, name = "Drug_ex_f (nmol)", mode = "lines", visible = vis[3]) %>%
    add_trace(y = ~Meta1_ex_f, name = "Meta1_ex_f (nmol)", mode = "lines", visible = vis[4]) %>%
    add_trace(y = ~Meta2_ex_f, name = "Meta2_ex_f (nmol)", mode = "lines", visible = vis[5]) %>%
    add_trace(y = ~V_tumor_pro_mm3, name = "V_tumor_pro_mm3", mode = "lines", visible = vis[6]) %>%
    add_trace(y = ~V_tumor_dyi_1_mm3, name = "V_tumor_dyi_1_mm3", mode = "lines", visible = vis[7]) %>%
    add_trace(y = ~V_tumor_dyi_2_mm3, name = "V_tumor_dyi_2_mm3", mode = "lines", visible = vis[8]) %>%
    add_trace(y = ~V_tumor_dyi_3_mm3, name = "V_tumor_dyi_3_mm3", mode = "lines", visible = vis[9]) %>%
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


#
#################### Databook feature ---------------------
#

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
