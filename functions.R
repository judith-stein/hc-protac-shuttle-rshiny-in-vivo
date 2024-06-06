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
