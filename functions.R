IsNaturalOrZero <- function(x, tol = .Machine$double.eps^0.5) {
  x == 0 | (x > tol & abs(x - round(x)) < tol)
}


SetColumnNameswithUnits <- function(df) {
  stopifnot(
    identical(
      colnames(df),
      c(
        "time_d",
        "time_h",
        "DAR",
        "V_tumor_mm3",
        "Ag_C1_b",
        "Ag_C1_f",
        "Ag_ex_b",
        "Ag_ex_f",
        "Ab_ex_f",
        "Ab_cell_b_ag",
        "Ab_TMDD_cell_b_ag",
        "ADC_C1_f",
        "ADC_C2_f",
        "Drug_C1_f",
        "Drug_C2_f",
        "Drug_C1_b_ntp",
        "ADC_ex_f",
        "Drug_ex_f",
        "ADC_cell_b_ag",
        "ADC_cell_lyso_f",
        "Drug_cell_lyso_f",
        "Drug_cell_cyto_f",
        "Drug_cell_cyto_b_dt",
        "ADC_TMDD_cell_b_ag",
        "ADC_TMDD_cell_lyso_f",
        "Drug_TMDD_cell_lyso_f",
        "Drug_TMDD_cell_cyto_f",
        "V_tumor_pro_mm3",
        "V_tumor_dyi_1_mm3",
        "V_tumor_dyi_2_mm3",
        "V_tumor_dyi_3_mm3",
        "degADC",
        "ADC_C1_f_nM",
        "ADC_C2_f_nM",
        "Ab_C1_f_nM",
        "Ab_C2_f_nM",
        "Occupancy",
        "degADCperTime",
        "ADC_in_tumor"
      )
    )
  )
  
  colnames(df) <- c(
    "time d",
    "time h",
    "DAR 1",
    "V_tumor mm^3",
    "Ag_C1_b nmol/kg",
    "Ag_C1_f nmol/kg",
    "Ag_ex_b nM",
    "Ag_ex_f nM",
    "Ab_ex_f nM",
    "Ab_cell_b_ag",
    "Ab_TMDD_cell_b_ag",
    "ADC_C1_f nmol/kg",
    "ADC_C2_f nmol/kg",
    "Drug_C1_f nM",
    "Drug_C2_f nM",
    "Drug_C1_b_ntp nM",
    "ADC_ex_f nM",
    "Drug_ex_f nmol",
    "ADC_cell_b_ag 1",
    "ADC_cell_lyso_f 1",
    "Drug_cell_lyso_f 1",
    "Drug_cell_cyto_f 1",
    "Drug_cell_cyto_b_dt 1",
    "ADC_TMDD_cell_b_ag 1",
    "ADC_TMDD_cell_lyso_f 1",
    "Drug_TMDD_cell_lyso_f 1",
    "Drug_TMDD_cell_cyto_f 1",
    "V_tumor_pro_mm3",
    "V_tumor_dyi_1_mm3",
    "V_tumor_dyi_2_mm3",
    "V_tumor_dyi_3_mm3",
    "degADC nmol",
    "ADC_C1_f nM",
    "ADC_C2_f nM",
    "Ab_C1_f nM",
    "Ab_C2_f nM",
    "Occupancy %",
    "degADCperTime nmol",
    "ADC_in_tumor nmol"
  )
  
  return(df)
}


FormatSimRawDataForOutput <- function(raw) {
  names(raw)[names(raw) == "time"] <- "time_h"
  # Add days column as well
  formatted <- tibble::add_column(raw, time_d = unlist(raw["time_h"]) / 24, .before = 1)
  # Shorten the output for other than time columns
  formatted[-c(1, 2)] <- format(formatted[-c(1, 2)], scientific = TRUE)
  # Only whole days
  filtered <- formatted[unlist(lapply(formatted["time_d"], IsNaturalOrZero)), ]
  filtered <- SetColumnNameswithUnits(filtered)
  return(filtered)
}

FormatSimRawDataForOutput3D <- function(input, raw) {
  points <- expand.grid(raw$x, raw$y)
  names(points) <- c(paste(" X:", input$xAxisItem), paste(" Y:", input$yAxisItem))
  points[[paste(" Z:", input$zAxisItem)]] <- as.vector(t(raw[[input$zAxisItem]]))
  
  return(points)
}

LoadData <- function(inputString) {
  # Load tab-limited input string into a table:
  outputTable <- read.table(
    inputString,
    row.names = NULL,
    header = FALSE,
    stringsAsFactors = FALSE
  )
  return(outputTable)
}

LoadDataDefault <- function(parameterDataFileName) {
  defaultAsString <- readChar(
    parameterDataFileName,
    file.info(parameterDataFileName)$size
  )
  return(defaultAsString)
}

# Eingabe von Experimental Data
LoadExpData <- function(inputString) {
  # Load space-delimited input string into a dataframe:
  #   Time (days) <space> ADC_C1_f_ug_ml <space> deviation
  expDataInput <- read.table(text = inputString, dec = ",")
  names(expDataInput) <- c("time", "V_tumor_mm3", "dev")
  return(expDataInput)
}


StringToNumeric <- function(inputString) {
  # Safe(er) way to evaluate string that includes arithmetic expressions
  # as numeric. Eg. "10^9/6.023*10^-23" --> 1.66030217499585e-15
  
  # Try first if string is just a normal number (eg. 21 or 1e2)
  output <- tryCatch(
    as.numeric(inputString),
    warning = function(w) w
  )
  
  # If as.numeric raised an error because there were
  # arithmetic operators like * or ^, get only these functions
  # from base R to an empty environment and evaluate the string there
  if (inherits(output, "simpleWarning")) {
    safeFunctions <- c(
      getGroupMembers("Arith"), # include * and ^ etc.
      "(" # also include () if input is eg. 10^(-2)
    )
    
    safeEnv <- new.env(parent = emptyenv())
    
    for (f in safeFunctions) {
      safeEnv[[f]] <- get(f, "package:base")
    }
    
    # If the input had something else than valid arithmetic operations,
    # we return NA
    output <- tryCatch(
      {
        evaluatedExpression <- eval(parse(text = inputString), env = safeEnv)
        output <- as.numeric(evaluatedExpression)
      },
      error = function(e) e
    )
    
    if (inherits(output, "simpleError")) {
      output <- NA
    }
  }
  
  return(output)
}


PickValuesToNumeric <- function(inputDataFrame) {
  # Take Values column from input table and return
  # it as named numeric vector
  
  outputVector <- as.character(inputDataFrame[, 2])
  names(outputVector) <- inputDataFrame[, 1]
  
  # Character to numeric, also translates limited math operations
  # eg. "10^9/(6.023*10^23)" -> 1.660302e-15
  outputVector <- sapply(
    outputVector,
    function(x) StringToNumeric(x)
  )
}


numericInput <- function(inputId, label, value, min = NA, max = NA, step = "any", width = NULL) {
  # Wrapper for numericInput to change the default arguments
  shiny::numericInput(inputId, label, value, min, max, step, width)
}


ParametersToCharacter <- function(mfitTheseParameters, minitialParameters) {
  # Take only the initial parameters that are to be fitted
  # and turn to character
  lengthOfFitParameterList <- length(mfitTheseParameters)
  fitTheseParameters <- character(lengthOfFitParameterList)
  
  for (i in 1:lengthOfFitParameterList) {
    name <- mfitTheseParameters[i]
    
    if (name %in% names(minitialParameters)) {
      fitTheseParameters[i] <- minitialParameters[name]
      names(fitTheseParameters)[i] <- name
    }
  }
  
  return(fitTheseParameters)
}


TakeVectorOfResiduals <- function(parametersToEstimate, fitTheseParameters, fitting, 
                                  currentExperimentalData,
                                  initialConditions, initialParameters,
                                  dataToFit, timeSteps, et) {
  # returns residuals, from which lm will minimize the sum of squares
  
  # time points for which conc. is reported
  # include the points where data is available
  times <- c(timeSteps, currentExperimentalData$time*24)
  times <- sort(unique(times))
  
  # EventTable  
  et <- etRep(et, samples="clear") # clear sampling times, but keep dosing
  et$add.sampling(times)
  
  names(parametersToEstimate) <- names(fitTheseParameters)
  
  newParameters <- replace(
    initialParameters,
    names(parametersToEstimate),
    parametersToEstimate
  )
  
  # solve ODE for a given set of parameters 
  outdf <- as.data.frame(Model$solve(newParameters,
                                     et,
                                     initialConditions)
  )
  
  # Filter data that contains time points where data is available
  outdf <- outdf[outdf$time %in% (currentExperimentalData$time*24), ]
  
  # Evaluate predicted vs experimental residual
  
  # Both amount and conc. by stacking scaled totalDRUG and concDRUG to one column:
  stacked_outdf <- stack(outdf, select = dataToFit) # Time column is not needed
  stacked_experimentalData <- stack(currentExperimentalData, select = dataToFit)
  
  if (fitting == "Local (fast)") {
    # Stack has columns 'values' and 'ind'
    ssqres <- abs(stacked_outdf$values - stacked_experimentalData$values)/abs(stacked_experimentalData$values)
  } else { # fitting == "Global"
    ssqres <- sum((abs(stacked_outdf$values - stacked_experimentalData$values)/abs(stacked_experimentalData$values))^2)
  }
  
  # return predicted vs experimental residual
  return(ssqres)
}


PrintSummaryIfFittingExists <- function(results) {
  # Print summary if fitting exists, otherwise return nothing
  
  if (!is.null(results)) {
    printResult <- summary(results)
    return(printResult)
  }
}


ParseJson <- function(input) {
  output <- tryCatch(
    {
      jsonlite::fromJSON(input)
    },
    error = function(e) stop(safeError(e))
  )
  return(output)
}


ReadStepsFromJson <- function(possibleJson) {
  parsedJson <- ParseJson(possibleJson)
  fields <- names(parsedJson)
  steps <- NA
  if ("values" %in% fields) {
    values <- parsedJson[["values"]]
    if (is.numeric(values)) {
      steps <- length(values)
    }
  } else if (
    all(c("xMin", "xMax", "xSteps") %in% fields)
  ) {
    steps <- parsedJson[["xSteps"]]
  }
  if (any(is.na(steps))) {
    stop(
      safeError(
        paste(
          "Wrong kind of text input. See examples here and try again:\n\n",
          defaultXAxisInstructions,
          "\n\nor\n\n",
          defaultYAxisInstructions
        )
      )
    )
  }
  return(steps)
}


ReadRangeFromJson <- function(possibleJson) {
  parsedJson <- ParseJson(possibleJson)
  fields <- names(parsedJson)
  range <- NA
  if ("values" %in% fields) {
    values <- parsedJson[["values"]]
    if (is.numeric(values)) {
      range <- sort(values)
    }
  } else if (
    all(c("xMin", "xMax", "xSteps") %in% fields)
  ) {
    range <- seq(
      from = parsedJson[["xMin"]],
      to = parsedJson[["xMax"]],
      length.out = parsedJson[["xSteps"]]
    )
  }
  if (any(is.na(range))) {
    stop(
      safeError(
        paste(
          "Wrong kind of text input. See examples here and try again:\n\n",
          defaultXAxisInstructions,
          "\n\nor\n\n",
          defaultYAxisInstructions
        )
      )
    )
  }
  return(range)
}


ReplaceIfThere <- function(namedList, name, newValue) {
  # Replaces entry in a named list if there is already
  # an entry with that name
  if (name %in% names(namedList)) {
    namedList[[name]] <- newValue
  }
  return(namedList)
}

BuildModelForN <- function(max = 1) {
  
  Ab_C1_b1KD <- ifelse(max >= 2, 
                     paste0("- K_Ab_Drug_on * (max - 1) * Ab_C1_b1 * Drug_C1_f + K_Ab_Drug_off * Ab_C1_b2"), 
                      paste0(""))
  
  Ab_C2_b1KD <- ifelse(max >= 2, 
                       paste0("- K_Ab_Drug_on * (max - 1) * Ab_C2_b1 * Drug_C2_f + K_Ab_Drug_off * Ab_C2_b2"), 
                       paste0(""))
  
  Ab_C1_b2KD <- ifelse(max > 2, 
                       paste0("- K_Ab_Drug_on * (max - 2) * Ab_C1_b2 * Drug_C1_f + K_Ab_Drug_off * Ab_C1_b3
                              "), 
                       paste0(""))
  Ab_C1_b2 <- ifelse(max >= 2, 
                     paste0("# Amount (nmol/kg) of Antibody bound to 2 Protacs in central compartment/plasma
d/dt(Ab_C1_b2) <- - (CL_Ab / V_C1_Ab) * Ab_C1_b2 
                  - (CLD_Ab / V_C1_Ab) * Ab_C1_b2 
                  + (CLD_Ab / V_C2_Ab) * Ab_C2_b2 
                  - (Ab_C1_b2 / V_C1_Ab - Ab_ex_b2 / E_Ab) * V_tumor / BW * (2 * P_Ab * R_Cap / R_Krogh^2 + 6 * D_Ab / R_Tumor^2)
                  + K_Ab_Drug_on * (max - 1) * Ab_C1_b1 * Drug_C1_f
                  - K_Ab_Drug_off * Ab_C1_b2
                  ", Ab_C1_b2KD
), 
                     paste0(""))

  Ab_C1_b3KD <- ifelse(max > 3, 
                       paste0(" - K_Ab_Drug_on * (max - 3) * Ab_C1_b3 * Drug_C1_f + K_Ab_Drug_off * Ab_C1_b4
                              "), 
                       paste0(""))
  Ab_C1_b3 <- ifelse(max >= 3, 
                     paste0("# Amount (nmol/kg) of Antibody bound to 3 Protacs in central compartment/plasma
d/dt(Ab_C1_b3) <- - (CL_Ab / V_C1_Ab) * Ab_C1_b3 
                  - (CLD_Ab / V_C1_Ab) * Ab_C1_b3 
                  + (CLD_Ab / V_C2_Ab) * Ab_C2_b3 
                  - (Ab_C1_b3 / V_C1_Ab - Ab_ex_b3 / E_Ab) * V_tumor / BW * (2 * P_Ab * R_Cap / R_Krogh^2 + 6 * D_Ab / R_Tumor^2)
                  + K_Ab_Drug_on * (max - 2) * Ab_C1_b2 * Drug_C1_f
                  - K_Ab_Drug_off * Ab_C1_b3
                  ", Ab_C1_b3KD
), 
                     paste0(""))
  
  Ab_C2_b2KD <- ifelse(max > 2, 
                       paste0(" - K_Ab_Drug_on * (max - 2) * Ab_C2_b2 * Drug_C2_f + K_Ab_Drug_off * Ab_C2_b3
                              "), 
                       paste0(""))
  Ab_C2_b2 <- ifelse(max >= 2, 
                     paste0("# Amount (nmol/kg) of Antibody bound to 2 Protacs in peripheral compartment
d/dt(Ab_C2_b2) <- (CLD_Ab / V_C1_Ab) * Ab_C1_b2 
                  - (CLD_Ab / V_C2_Ab) * Ab_C2_b2 
                  + K_Ab_Drug_on * (max - 1) * Ab_C2_b1 * Drug_C2_f
                  - K_Ab_Drug_off * Ab_C2_b2
                  ", Ab_C2_b2KD
), 
                     paste0(""))

  Ab_C2_b3KD <- ifelse(max > 3, 
                       paste0(" - K_Ab_Drug_on * (max - 3) * Ab_C2_b3 * Drug_C2_f + K_Ab_Drug_off * Ab_C2_b4
                              "), 
                       paste0(""))
  Ab_C2_b3 <- ifelse(max >= 3, 
                     paste0("# Amount (nmol/kg) of Antibody bound to 3 Protacs in peripheral compartment
d/dt(Ab_C2_b3) <- (CLD_Ab / V_C1_Ab) * Ab_C1_b3 
                  - (CLD_Ab / V_C2_Ab) * Ab_C2_b3 
                  + K_Ab_Drug_on * (max - 2) * Ab_C2_b2 * Drug_C2_f
                  - K_Ab_Drug_off * Ab_C2_b3
                  ", Ab_C2_b3KD
), 
                     paste0(""))

  Ab_C1_b4 <- ifelse(max == 4, 
                     paste0("# Amount (nmol/kg) of Antibody bound to 4 Protacs in central compartment/plasma
d/dt(Ab_C1_b4) <- - (CL_Ab / V_C1_Ab) * Ab_C1_b4 
                  - (CLD_Ab / V_C1_Ab) * Ab_C1_b4 
                  + (CLD_Ab / V_C2_Ab) * Ab_C2_b4 
                  - (Ab_C1_b4 / V_C1_Ab - Ab_ex_b4 / E_Ab) * V_tumor / BW * (2 * P_Ab * R_Cap / R_Krogh^2 + 6 * D_Ab / R_Tumor^2)
                  + K_Ab_Drug_on * Ab_C1_b3 * Drug_C1_f
                  - K_Ab_Drug_off * Ab_C1_b4"), 
                     paste0(""))
  
  Ab_C2_b4 <- ifelse(max == 4, 
                     paste0("# Amount (nmol/kg) of Antibody bound to 4 Protacs in peripheral compartment
d/dt(Ab_C2_b4) <- (CLD_Ab / V_C1_Ab) * Ab_C1_b4 
                  - (CLD_Ab / V_C2_Ab) * Ab_C2_b4 
                  + K_Ab_Drug_on * Ab_C2_b3 * Drug_C2_f
                  - K_Ab_Drug_off * Ab_C2_b4"), 
                     paste0(""))
  
  Drug_C1_f <- ifelse(max >= 1, 
                       paste0("CL_Ab * ", 1:max," * (Ab_C1_b", 1:max," / V_C1_Ab) / V_C1_Drug", collapse=" + "), 
                       paste0(""))
  
  Ab_ex_f <- ifelse(max >= 2, 
                    paste0("Ab_cell_b", 2:max,"_b_ag", collapse=" + "), 
                    paste0("0"))
  
  Ab_ex_b1KD <- ifelse(max >= 2, 
                       paste0(" - K_Ab_Drug_on * (max - 1) * Ab_ex_b1 / E_Ab * Drug_ex_f / V_tumor + K_Ab_Drug_off * Ab_ex_b2 / E_Ab
                              "), 
                       paste0(""))

  Ab_ex_b2KD <- ifelse(max > 2, 
                       paste0(" - K_Ab_Drug_on * (max - 2) * Ab_ex_b2 / E_Ab * Drug_ex_f / V_tumor + K_Ab_Drug_off * Ab_ex_b3 / E_Ab
                              "), 
                       paste0(""))
  Ab_ex_b2 <- ifelse(max >= 2, 
                     paste0("# Concentration (nM) of Antibody bound to 2 Protacs in tumor extracellular space
d/dt(Ab_ex_b2) <- (Ab_C1_b2 / V_C1_Ab - Ab_ex_b2 / E_Ab) * (2 * P_Ab * R_Cap / R_Krogh^2 + 6 * D_Ab / R_Tumor^2)
                  + (- K_Ab_cell_ag_on * Ab_ex_b2 / E_Ab * (Ag_cell_t - Ab_cell_f_b_ag - Ab_cell_b1_b_ag - (", Ab_ex_f, ")) 
                  + K_Ab_cell_ag_off * Ab_cell_b2_b_ag) * NC_Tumor * SF / V_tumor 
                  + 1 / tau * V_tumor_dyi_3_mm3 * 10^5 * (Ab_cell_b2_b_ag + Ab_cell_lyso_b2) * SF / V_tumor 
                  - K_Ab_cell_lyso_pino * NCL_tumor / E_Ab * Ab_ex_b2 
                  + K_Ab_Drug_on * (max - 1) * Ab_ex_b1 / E_Ab * Drug_ex_f / V_tumor 
                  - K_Ab_Drug_off * Ab_ex_b2 / E_Ab
                  ", Ab_ex_b2KD), 
                     paste0(""))

  Ab_ex_b3KD <- ifelse(max > 3, 
                       paste0(" - K_Ab_Drug_on * (max - 3) * Ab_ex_b3 / E_Ab * Drug_ex_f / V_tumor + K_Ab_Drug_off * Ab_ex_b4 / E_Ab
                              "), 
                       paste0(""))
  Ab_ex_b3 <- ifelse(max >= 3, 
                     paste0("# Concentration (nM) of Antibody bound to 3 Protacs in tumor extracellular space
d/dt(Ab_ex_b3) <- (Ab_C1_b3 / V_C1_Ab - Ab_ex_b3 / E_Ab) * (2 * P_Ab * R_Cap / R_Krogh^2 + 6 * D_Ab / R_Tumor^2)
                  + (- K_Ab_cell_ag_on * Ab_ex_b3 / E_Ab * (Ag_cell_t - Ab_cell_f_b_ag - Ab_cell_b1_b_ag - (", Ab_ex_f, ")) 
                  + K_Ab_cell_ag_off * Ab_cell_b3_b_ag) * NC_Tumor * SF / V_tumor 
                  + 1 / tau * V_tumor_dyi_3_mm3 * 10^5 * (Ab_cell_b3_b_ag + Ab_cell_lyso_b3) * SF / V_tumor 
                  - K_Ab_cell_lyso_pino * NCL_tumor / E_Ab * Ab_ex_b3 
                  + K_Ab_Drug_on * (max - 2) * Ab_ex_b2 / E_Ab * Drug_ex_f / V_tumor 
                  - K_Ab_Drug_off * Ab_ex_b3 / E_Ab
                  ", Ab_ex_b3KD), 
                     paste0(""))

  Ab_ex_b4 <- ifelse(max == 4, 
                     paste0("# Concentration (nM) of Antibody bound to 4 Protacs in tumor extracellular space
d/dt(Ab_ex_b4) <- (Ab_C1_b4 / V_C1_Ab - Ab_ex_b4 / E_Ab) * (2 * P_Ab * R_Cap / R_Krogh^2 + 6 * D_Ab / R_Tumor^2)
                  + (- K_Ab_cell_ag_on * Ab_ex_b4 / E_Ab * (Ag_cell_t - Ab_cell_f_b_ag - Ab_cell_b1_b_ag - (", Ab_ex_f, ")) 
                  + K_Ab_cell_ag_off * Ab_cell_b4_b_ag) * NC_Tumor * SF / V_tumor 
                  + 1 / tau * V_tumor_dyi_3_mm3 * 10^5 * (Ab_cell_b4_b_ag + Ab_cell_lyso_b4) * SF / V_tumor 
                  - K_Ab_cell_lyso_pino * NCL_tumor / E_Ab * Ab_ex_b4 
                  + K_Ab_Drug_on * Ab_ex_b3 / E_Ab * Drug_ex_f / V_tumor 
                  - K_Ab_Drug_off * Ab_ex_b4 / E_Ab"), 
                     paste0(""))
  
  Ab_cell_b1KD <- ifelse(max >= 2, 
                       paste0(" - K_Ab_Drug_on * (max - 1) * Ab_cell_b1_b_ag * Drug_ex_f / V_tumor + K_Ab_Drug_off * Ab_cell_b2_b_ag"), 
                       paste0(""))
  
  Ab_cell_b2KD <- ifelse(max > 2, 
                         paste0(" - K_Ab_Drug_on * (max - 2) * Ab_cell_b2_b_ag * Drug_ex_f / V_tumor + K_Ab_Drug_off * Ab_cell_b3_b_ag
                                "), 
                         paste0(""))
  Ab_cell_b2_b_ag <- ifelse(max >= 2, 
                     paste0("# Number of Antibody molecules bound to 2 Protacs bound to binding target on a single tumor cell
d/dt(Ab_cell_b2_b_ag) <- K_Ab_cell_ag_on * Ab_ex_b2 / E_Ab * (Ag_cell_t - Ab_cell_f_b_ag - Ab_cell_b1_b_ag - (", Ab_ex_f, ")) 
                         - K_Ab_cell_ag_off * Ab_cell_b2_b_ag 
                         - K_Ab_cell_int * Ab_cell_b2_b_ag 
                         - log(2) / DT_tumor * Ab_cell_b2_b_ag  
                         + K_Ab_Drug_on * (max - 1) * Ab_cell_b1_b_ag * Drug_ex_f / V_tumor 
                         - K_Ab_Drug_off * Ab_cell_b2_b_ag
                         ", Ab_cell_b2KD), 
                     paste0(""))

  Ab_cell_b3KD <- ifelse(max > 3, 
                         paste0(" - K_Ab_Drug_on * (max - 3) * Ab_cell_b3_b_ag * Drug_ex_f / V_tumor + K_Ab_Drug_off * Ab_cell_b4_b_ag
                                "), 
                         paste0(""))
  Ab_cell_b3_b_ag <- ifelse(max >= 3, 
                            paste0("# Number of Antibody molecules bound to 3 Protacs bound to binding target on a single tumor cell
d/dt(Ab_cell_b3_b_ag) <- K_Ab_cell_ag_on * Ab_ex_b3 / E_Ab * (Ag_cell_t - Ab_cell_f_b_ag - Ab_cell_b1_b_ag - (", Ab_ex_f, ")) 
                         - K_Ab_cell_ag_off * Ab_cell_b3_b_ag 
                         - K_Ab_cell_int * Ab_cell_b3_b_ag 
                         - log(2) / DT_tumor * Ab_cell_b3_b_ag  
                         + K_Ab_Drug_on * (max - 2) * Ab_cell_b2_b_ag * Drug_ex_f / V_tumor 
                         - K_Ab_Drug_off * Ab_cell_b3_b_ag
                         ", Ab_cell_b3KD), 
                            paste0(""))
  
  Ab_cell_b4_b_ag <- ifelse(max == 4, 
                            paste0("# Number of Antibody molecules bound to 4 Protacs bound to binding target on a single tumor cell
d/dt(Ab_cell_b4_b_ag) <- K_Ab_cell_ag_on * Ab_ex_b4 / E_Ab * (Ag_cell_t - Ab_cell_f_b_ag - Ab_cell_b1_b_ag - (", Ab_ex_f, ")) 
                         - K_Ab_cell_ag_off * Ab_cell_b4_b_ag 
                         - K_Ab_cell_int * Ab_cell_b4_b_ag 
                         - log(2) / DT_tumor * Ab_cell_b4_b_ag  
                         + K_Ab_Drug_on * Ab_cell_b3_b_ag * Drug_ex_f / V_tumor 
                         - K_Ab_Drug_off * Ab_cell_b4_b_ag"), 
                            paste0(""))
  
  Ab_cell_lyso_bi <- paste0("# Number of Antibody molecules bound to ", 1:max, " Protacs in lysosomal space on a single tumor cell
d/dt(Ab_cell_lyso_b", 1:max, ") <- K_Ab_cell_int * Ab_cell_b", 1:max, "_b_ag 
                              - K_Ab_deg * Ab_cell_lyso_b", 1:max, " 
                              + K_Ab_cell_lyso_pino * Ab_ex_b", 1:max, " / (E_Ab * SF) 
                              - log(2) / DT_tumor * Ab_cell_lyso_b", 1:max, 
                              "
                              ", collapse="\n")
  
  Drug_cell_lyso_f <- paste0("K_Ab_deg * Ab_cell_lyso_b", 1:max, " * ", 1:max, collapse=" + ")

  
  Model <- paste0("
 # Modeling assignemnts -----------------------------------

  SF ~ 10^9 * 1 / 6.023 * 10^-23

  # Add total tumor volume in mm³
  V_tumor_mm3 <- V_tumor_pro_mm3 + V_tumor_dyi_1_mm3 + V_tumor_dyi_2_mm3 + V_tumor_dyi_3_mm3
 
  # Tumor volume in liter
  V_tumor ~ V_tumor_mm3 * 10^-6
 
  # Number of tumor cells 
  NC_Tumor ~ NCL_tumor * V_tumor
  
  # Radius tumor in cm
  R_Tumor ~ (3 * V_tumor / (4 * pi))^(1/3) * 10
  
  
  # Ordinary differential equations ------------------------------
  ########################
  
  # Amount (nmol/kg) of free Antibody (bound to 0 Protacs) in central compartment/plasma
  d/dt(Ab_C1_f) <- - (CL_Ab / V_C1_Ab) * Ab_C1_f 
                    - (CLD_Ab / V_C1_Ab) * Ab_C1_f 
                    + (CLD_Ab / V_C2_Ab) * Ab_C2_f 
                    - (Ab_C1_f / V_C1_Ab - Ab_ex_f / E_Ab) * V_tumor / BW * (2 * P_Ab * R_Cap / R_Krogh^2 + 6 * D_Ab / R_Tumor^2)
                    - K_Ab_Drug_on * max * Ab_C1_f * Drug_C1_f
                    + K_Ab_Drug_off * Ab_C1_b1
  
  # Amount (nmol/kg) of Antibody bound to 1 Protac in central compartment/plasma
  d/dt(Ab_C1_b1) <- - (CL_Ab / V_C1_Ab) * Ab_C1_b1 
                    - (CLD_Ab / V_C1_Ab) * Ab_C1_b1 
                    + (CLD_Ab / V_C2_Ab) * Ab_C2_b1 
                    - (Ab_C1_b1 / V_C1_Ab - Ab_ex_b1 / E_Ab) * V_tumor / BW * (2 * P_Ab * R_Cap / R_Krogh^2 + 6 * D_Ab / R_Tumor^2)
                    + K_Ab_Drug_on * max * Ab_C1_f * Drug_C1_f
                    - K_Ab_Drug_off * Ab_C1_b1
                    {{Ab_C1_b1KD}}
 
 {{Ab_C1_b2}}
 {{Ab_C1_b3}}
 {{Ab_C1_b4}}

  # Concentration (nmol/l) of free (unbound) Drug in central compartment/plasma
  d/dt(Drug_C1_f) <- - CL_Drug / V_C1_Drug * Drug_C1_f 
                     - CLD_Drug / V_C1_Drug * Drug_C1_f 
                     + CLD_Drug / V_C1_Drug * Drug_C2_f 
                     - (Drug_C1_f - Drug_ex_f / (V_tumor * E_Drug)) * V_tumor / (V_C1_Drug * BW) * (2 * P_Drug * R_Cap / R_Krogh^2 + 6 * D_Drug / R_Tumor^2)
                     - K_Drug_ex_ntp_on_off * (1 - f_ex_ub) * Drug_C1_f 
                     + K_Drug_ex_ntp_on_off * f_ex_ub * Drug_C1_b_ntp
                     + {{Drug_C1_f}}
 
  # Concentration (nmol/l) of drug bound to unspecific protein in central compartment/plasma
  d/dt(Drug_C1_b_ntp) <- K_Drug_ex_ntp_on_off * (1 - f_ex_ub) * Drug_C1_f 
                        - K_Drug_ex_ntp_on_off * f_ex_ub * Drug_C1_b_ntp
 
  # Amount (nmol/kg) of free Antibody (bound to 0 Protacs) in peripheral compartment
  d/dt(Ab_C2_f) <- CLD_Ab / V_C1_Ab * Ab_C1_f 
                  - CLD_Ab / V_C2_Ab * Ab_C2_f
                   - K_Ab_Drug_on * max * Ab_C2_f * Drug_C2_f
                   + K_Ab_Drug_off * Ab_C2_b1
 
  # Amount (nmol/kg) of Antibody bound to 1 Protac in peripheral compartment
  d/dt(Ab_C2_b1) <- (CLD_Ab / V_C1_Ab) * Ab_C1_b1 
                    - (CLD_Ab / V_C2_Ab) * Ab_C2_b1 
                    + K_Ab_Drug_on * max * Ab_C2_f * Drug_C2_f
                    - K_Ab_Drug_off * Ab_C2_b1
                    {{Ab_C2_b1KD}}
 
 {{Ab_C2_b2}}
 {{Ab_C2_b3}}
 {{Ab_C2_b4}}
 
  # Concentration (nmol/l) of drug in peripheral compartment
  d/dt(Drug_C2_f) <- CLD_Drug / V_C2_Drug * Drug_C1_f
                    - CLD_Drug / V_C2_Drug * Drug_C2_f
 
  # Concentration (nM) of free Antibody (bound to 0 Protacs) in tumor extracellular space
  d/dt(Ab_ex_f) <- (Ab_C1_f / V_C1_Ab - Ab_ex_f / E_Ab) * (2 * P_Ab * R_Cap / R_Krogh^2 + 6 * D_Ab / R_Tumor^2)
                    + (- K_Ab_cell_ag_on * Ab_ex_f / E_Ab * (Ag_cell_t - Ab_cell_f_b_ag - Ab_cell_b1_b_ag - ({{Ab_ex_f}})) 
                    + K_Ab_cell_ag_off * Ab_cell_f_b_ag) * NC_Tumor * SF / V_tumor 
                    + 1 / tau * V_tumor_dyi_3_mm3 * 10^5 * (Ab_cell_f_b_ag) * SF / V_tumor 
                    - K_Ab_cell_lyso_pino * NCL_tumor / E_Ab * Ab_ex_f
                    - K_Ab_Drug_on * max * Ab_ex_f / E_Ab * Drug_ex_f / V_tumor
                    + K_Ab_Drug_off * Ab_ex_b1 / E_Ab
 
  # Concentration (nM) of Antibody bound to 1 Protac in tumor extracellular space
  d/dt(Ab_ex_b1) <- (Ab_C1_b1 / V_C1_Ab - Ab_ex_b1 / E_Ab) * (2 * P_Ab * R_Cap / R_Krogh^2 + 6 * D_Ab / R_Tumor^2)
                    + (- K_Ab_cell_ag_on * Ab_ex_b1 / E_Ab * (Ag_cell_t - Ab_cell_f_b_ag - Ab_cell_b1_b_ag - ({{Ab_ex_f}})) 
                    + K_Ab_cell_ag_off * Ab_cell_b1_b_ag) * NC_Tumor * SF / V_tumor 
                    + 1 / tau * V_tumor_dyi_3_mm3 * 10^5 * (Ab_cell_b1_b_ag + Ab_cell_lyso_b1) * SF / V_tumor 
                    - K_Ab_cell_lyso_pino * NCL_tumor / E_Ab * Ab_ex_b1
                    + K_Ab_Drug_on * max * Ab_ex_f / E_Ab * Drug_ex_f / V_tumor
                    - K_Ab_Drug_off * Ab_ex_b1 / E_Ab
                    {{Ab_ex_b1KD}}
 
 {{Ab_ex_b2}}
 {{Ab_ex_b3}}
 {{Ab_ex_b4}}
 
  # Amount (nmol) of drug in tumor extracellular space
  d/dt(Drug_ex_f) <- (Drug_C1_f - Drug_ex_f / (V_tumor * E_Drug)) * V_tumor * (2 * P_Drug * R_Cap / R_Krogh^2 + 6 * D_Drug / R_Tumor^2) 
                    + K_Drug_ex_out * Drug_cell_cyto_f * NC_Tumor * SF 
                    - K_Drug_ex_in * NC_Tumor * (V_cell / (V_tumor * E_Drug)) * Drug_ex_f 
                    + 1 / tau * V_tumor_dyi_3_mm3 * 10^5 * (Drug_cell_cyto_f + Drug_cell_cyto_b_dt + Drug_cell_lyso_f) * SF
 
  # Number of free Antibody (bound to 0 Protacs) molecules bound to binding target on a single tumor cell
  d/dt(Ab_cell_f_b_ag) <- K_Ab_cell_ag_on * Ab_ex_f / E_Ab * (Ag_cell_t - Ab_cell_f_b_ag - Ab_cell_b1_b_ag - ({{Ab_ex_f}})) 
                          - K_Ab_cell_ag_off * Ab_cell_f_b_ag 
                          - K_Ab_cell_int * Ab_cell_f_b_ag 
                          - log(2) / DT_tumor * Ab_cell_f_b_ag 
                          - K_Ab_Drug_on * max * Ab_cell_f_b_ag * Drug_ex_f / V_tumor
                          + K_Ab_Drug_off * Ab_cell_b1_b_ag 
 
  # Number of Antibody molecules bound to 1 Protac bound to binding target on a single tumor cell
  d/dt(Ab_cell_b1_b_ag) <- K_Ab_cell_ag_on * Ab_ex_b1 / E_Ab * (Ag_cell_t - Ab_cell_f_b_ag - Ab_cell_b1_b_ag - ({{Ab_ex_f}})) 
                          - K_Ab_cell_ag_off * Ab_cell_b1_b_ag 
                          - K_Ab_cell_int * Ab_cell_b1_b_ag 
                          - log(2) / DT_tumor * Ab_cell_b1_b_ag 
                          + K_Ab_Drug_on * max * Ab_cell_f_b_ag * Drug_ex_f / V_tumor
                          - K_Ab_Drug_off * Ab_cell_b1_b_ag 
                          {{Ab_cell_b1KD}}
 
 {{Ab_cell_b2_b_ag}}
 {{Ab_cell_b3_b_ag}}
 {{Ab_cell_b4_b_ag}}
 
 {{Ab_cell_lyso_bi}}

  # Number of unbound (free) drug molecules in lysosomal space on a single tumor cell
  d/dt(Drug_cell_lyso_f) <- - K_Drug_lyso_out * (V_cell / V_cell_lyso) * Drug_cell_lyso_f 
                            + K_Drug_lyso_in * Drug_cell_cyto_f 
                            - log(2) / DT_tumor * Drug_cell_lyso_f
                            + {{Drug_cell_lyso_f}}
 
  # Number of unbound (free) drug molecules in cytosol on a single tumor cell
  d/dt(Drug_cell_cyto_f) <- K_Drug_lyso_out * (V_cell / V_cell_lyso) * Drug_cell_lyso_f 
                            - K_Drug_lyso_in * Drug_cell_cyto_f 
                            - K_Drug_ex_out * Drug_cell_cyto_f 
                            - K_Drug_cyto_dt_on * SF / V_cell * Drug_cell_cyto_f * (DrugTarget_cell_cyto_t * V_cell / SF - Drug_cell_cyto_b_dt) 
                            - K_Drug_met * Drug_cell_cyto_f 
                            + K_Drug_cyto_dt_off * Drug_cell_cyto_b_dt 
                            + K_Drug_ex_in * (V_cell / (V_tumor * E_Drug)) * (Drug_ex_f / SF) 
                            - log(2) / DT_tumor * Drug_cell_cyto_f
 
  # Number of target-bound drug molecules in cytosol on a single tumor cell
  d/dt(Drug_cell_cyto_b_dt) <- K_Drug_cyto_dt_on * SF / V_cell * Drug_cell_cyto_f * (DrugTarget_cell_cyto_t * V_cell / SF - Drug_cell_cyto_b_dt) 
                              - K_Drug_cyto_dt_off * Drug_cell_cyto_b_dt 
                              - log(2) / DT_tumor * Drug_cell_cyto_b_dt
 
  # Logistic killing function from Thomas Rysiok ------------------------------
  ########################
 
  if (!is.na(Drug_cell_cyto_f) && Drug_cell_cyto_f > 0) {
    nM_Drug_cell_cyto_f ~ Drug_cell_cyto_f * SF / V_cell
  } else {nM_Drug_cell_cyto_f ~ 0}  # Otherwise log of -1 gets NA
 
  tl ~ log(nM_Drug_cell_cyto_f) - log(EC50)
  LOGI ~ k_g / (1 + (k_g / k_z - 1) * exp(-k_r * k_g * tl))
  R_Kill ~ k_kill_max * (log(2) / DT_tumor)^f_DT_kill * LOGI
 
  Kg ~ log(2) / DT_tumor * (1 - V_tumor_pro_mm3 / V_tumor_max) / (1 + (log(2) / DT_tumor * V_tumor_pro_mm3 / k_lin)^psi)^(1/psi)
  
  d/dt(V_tumor_pro_mm3) <- (Kg - R_Kill) * V_tumor_pro_mm3
  d/dt(V_tumor_dyi_1_mm3) <- R_Kill * V_tumor_pro_mm3 - 1 / tau * V_tumor_dyi_1_mm3
  d/dt(V_tumor_dyi_2_mm3) <- 1 / tau * (V_tumor_dyi_1_mm3 - V_tumor_dyi_2_mm3)
  d/dt(V_tumor_dyi_3_mm3) <- 1 / tau * (V_tumor_dyi_2_mm3 - V_tumor_dyi_3_mm3)
  ")
  Model <- glue(Model, .open = "{{", .close = "}}")
  return(Model)
}
