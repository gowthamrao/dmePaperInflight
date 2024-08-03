#' Impute Missing Years in a Data Frame
#'
#' This function imputes missing years in a data frame. It groups the data by cohortId, databaseId, gender, and ageGroup, and fills in missing years with specified values.
#' The function also adds a flag to indicate whether the data was imputed.
#'
#' @param df A data frame that contains the columns cohortId, databaseId, gender, ageGroup, and calendarYear.
#' @return A data frame with missing years imputed.
#' @export
imputeMissingYears <- function(df) {
  df_with_flag <- df |>
    dplyr::mutate(imputedData = as.integer(0))
  
  output <- df_with_flag |>
    dplyr::group_by(cohortId, databaseId, gender, ageGroup) |>
    tidyr::complete(
      calendarYear = seq(min(calendarYear), max(calendarYear), by = 1),
      fill = list(
        cohortCount = 0,
        personYears = 0,
        incidenceRate = 0,
        imputedData = 1
      )
    ) |>
    dplyr::mutate(imputedData = ifelse(is.na(imputedData), 0, imputedData)) |>
    dplyr::ungroup() |>
    dplyr::arrange(cohortId, databaseId, gender, ageGroup, calendarYear)
  
  return(output)
}


#' Get Predicted Count in a Data Frame
#'
#' This function predicts the count of events in a data frame using a Poisson regression model with natural splines. It groups the data by timeSequenceField and fills in missing counts with predicted values.
#' The function also adds a flag to indicate whether the data was imputed.
#'
#' @param data A data frame that contains the columns specified in the function parameters.
#' @param timeSequenceField The field in the data frame that represents the time sequence.
#' @param countField The field in the data frame that represents the count of events.
#' @param personTimeField The field in the data frame that represents the person time.
#' @param maxNumberOfSplines The maximum number of splines to use in the model. If NULL, the default is 5.
#' @param splineTickInterval The interval at which to place the splines. The default is 3.
#' @param zeroCountAdjustment A logical value indicating whether to adjust for zero counts. If TRUE, a small constant is added to avoid log(0).
#' @return A data frame with the observed and expected counts.
#' @export
getPredictedCount <- function(data,
                              timeSequenceField,
                              countField,
                              personTimeField,
                              maxNumberOfSplines = NULL,
                              splineTickInterval = 3,
                              zeroCountAdjustment = TRUE) {
  # Check if there are duplicate records for the same timeSequenceField
  if (length(data[[timeSequenceField]]) != length(data[[timeSequenceField]] |> unique())) {
    stop("Cant have more than one record per ", timeSequenceField)
  }
  
  # Create a new column 'observed' based on the countField
  data$observed <- data[[countField]]
  
  # Arrange the data based on timeSequenceField and create a new column 'timeId'.
  # timeSequenceField is the natural interval order in the source data. timeId is the index of the record in the ordered data.
  data <- data |>
    dplyr::arrange(timeSequenceField) |>
    dplyr::mutate(timeId = dplyr::row_number())
  
  # If all observed values are 0, return the data with 'expected' column filled with 0s
  if (all(data$observed == 0)) {
    data$expected <- 0
    return(data)
  }
  
  # Set the default value for maxNumberOfSplines if it's NULL
  if (is.null(maxNumberOfSplines)) {
    maxNumberOfSplines <- 5
  }
  
  # Determine the number of splines based on the number of unique timeId values
  numberOfSplines <-
    max(maxNumberOfSplines, ceiling(length(unique(data$timeId)) / splineTickInterval)) # we will put a spline for every 3 ticks
  
  # Create a Cyclops data object using the observed counts, timeId, and personTimeField
  # The formula used here is a Poisson regression model with natural splines
  offsetTerm <- if (zeroCountAdjustment) {
    paste0("offset(log(", personTimeField, " + 1e-8))") # adding small constant to avoid log(0)
  } else {
    paste0("offset(log(", personTimeField, "))")
  }
  model <- Cyclops::createCyclopsData(formula =
                                        as.formula(
                                          paste(
                                            "observed ~ splines::ns(x = timeId, df =",
                                            numberOfSplines,
                                            ") + ",
                                            offsetTerm
                                          )
                                        ),
                                      data = data,
                                      modelType = "pr") |>
    Cyclops::fitCyclopsModel() # Fit the Cyclops model
  
  # Predict the counts using the fitted model
  predictions <- stats::predict(model, newdata = data)
  
  # Add the predicted counts to the data
  data$expected <- round(as.double(predictions))
  
  # Arrange the data based on timeId and select the 'observed' and 'expected' columns
  data <- data |>
    dplyr::arrange(timeId) |>
    dplyr::select(-timeId) |>
    dplyr::tibble() |>
    dplyr::select("observed", "expected")
  
  #replace NA in expected with observed values
  data <- data |>
    dplyr::mutate(expected = ifelse(is.na(expected), observed, expected))
  
  return(data)
}




#' Compare Likelihoods in a Data Frame
#'
#' This function compares the likelihoods of observed and expected counts in a data frame. It uses a Poisson regression model with natural splines and calculates the log-likelihood ratio. The function also adds a flag to indicate whether the data was imputed.
#'
#' @param data A data frame that contains the columns 'observed' and 'expected'.
#' @param maxRatio The maximum ratio for the likelihood comparison. If NULL, the default is 1.25.
#' @param alpha The significance level for the likelihood ratio test. If NULL, the default is 0.05.
#' @param zeroCountAdjustment A logical value indicating whether to adjust for zero counts. If TRUE, a small constant is added to avoid log(0).
#' @return A data frame with the ratio, p-value, and stability indicator.
#' @export
likelihoodComparison <- function(data,
                                 maxRatio = 1.25,
                                 alpha = 0.05,
                                 zeroCountAdjustment = TRUE) {
  observed <- data$observed
  expected <- data$expected
  # From https://cdsmithus.medium.com/the-logarithm-of-a-sum-69dd76199790
  smoothMax <- function(x, y, zeroCountAdjustment) {
    smoothMax <- ifelse(
      test = abs(x - y) > 100,
      yes = pmax(x, y),
      no = if (zeroCountAdjustment) {
        x + log(1 + exp(y - x + 1e-8)) # adding a small constant
      } else {
        x + log(1 + exp(y - x))
      }
    )
    return(smoothMax)
  }
  
  logLikelihood <- function(x, zeroCountAdjustment) {
    xTerm <- if (zeroCountAdjustment) {
      dpois(observed, expected * x + 1e-8, log = TRUE) # adding a small constant
    } else {
      dpois(observed, expected * x, log = TRUE)
    }
    yTerm <- if (zeroCountAdjustment) {
      dpois(observed, expected / x + 1e-8, log = TRUE) # adding a small constant
    } else {
      dpois(observed, expected / x, log = TRUE)
    }
    return(-sum(
      smoothMax(
        x = xTerm,
        y = yTerm,
        zeroCountAdjustment = zeroCountAdjustment
      )
    ))
  }
  
  likelihood <- function(x, zeroCountAdjustment) {
    return(exp(-logLikelihood(x, zeroCountAdjustment)))
  }
  
  vectorLikelihood <- function(x, zeroCountAdjustment) {
    return(sapply(x, likelihood, zeroCountAdjustment = zeroCountAdjustment))
  }
  
  x <- seq(1, 10, by = 0.1)
  ll <- sapply(x, logLikelihood, zeroCountAdjustment = zeroCountAdjustment)
  
  indices <- which(!is.na(ll) & !is.infinite(ll))
  if (length(indices) > 0) {
    maxX <- x[max(indices)]
    minX <- x[min(indices)]
  } else {
    # Handle the case when indices is empty by setting to NA
    maxX <- NA
    minX <- NA
  }
  
  xHat <- optim(
    1.5,
    logLikelihood,
    lower = minX,
    upper = maxX,
    method = "L-BFGS-B",
    zeroCountAdjustment = zeroCountAdjustment
  )$par
  
  l0 <- integrate(
    vectorLikelihood,
    lower = 1,
    upper = maxRatio,
    zeroCountAdjustment = zeroCountAdjustment
  )$value
  l1 <- integrate(
    vectorLikelihood,
    lower = maxRatio,
    upper = Inf,
    zeroCountAdjustment = zeroCountAdjustment
  )$value
  
  llr <- 2 * (log(l1) - log(l0))
  if (is.nan(llr)) {
    if (xHat > maxRatio) {
      p <- 0
    } else {
      p <- 1
    }
  } else {
    p <- pchisq(llr, 1, lower.tail = FALSE)
  }
  result <- dplyr::tibble(ratio = xHat,
                          p = p,
                          stable = p > alpha)
  return(result)
}


#' Check If Data is Cohort Diagnostics Incidence Rate Data
#'
#' This function checks if the input data is a data frame and contains all the expected columns with the correct data types for cohort diagnostics incidence rate data.
#'
#' @param data A data frame that is expected to contain the columns cohortCount, personYears, gender, ageGroup, calendarYear, incidenceRate, cohortId, and databaseId.
#' @return This function does not return a value. It stops and throws an error if the input data does not meet the expected conditions.
#' @export
checkIfCohortDiagnosticsIncidenceRateData <- function(data) {
  # Check if data is a data frame
  if (!is.data.frame(data)) {
    stop("Input data must be a data frame.")
  }
  
  # Check if data has all the expected columns
  expected_columns <- c(
    "cohortCount",
    "personYears",
    "gender",
    "ageGroup",
    "calendarYear",
    "incidenceRate",
    "cohortId",
    "databaseId"
  )
  if (!all(expected_columns %in% names(data))) {
    stop(paste(
      "Input data must contain the following columns:",
      paste(expected_columns, collapse = ", ")
    ))
  }
  
  # Check data types
  if (!is.numeric(data$cohortCount) |
      !is.numeric(data$personYears) |
      !is.numeric(data$incidenceRate) |
      !is.numeric(data$cohortId) |
      !is.character(data$gender) |
      !is.character(data$ageGroup) |
      !is.character(data$calendarYear) |
      !is.character(data$databaseId)) {
    stop("Data types of some columns are not as expected.")
  }
}


#' Process Cohort Diagnostics Incidence Rate Data
#'
#' This function processes the cohort diagnostics incidence rate data. It filters the data based on certain conditions, calculates the incidence rate, and imputes missing years.
#'
#' @param cohortDiagnosticsIncidenceRateData A data frame that contains the columns cohortId, databaseId, gender, ageGroup, calendarYear, cohortCount, personYears, and incidenceRate.
#' @param zeroCountAdjustment A logical value indicating whether to adjust for zero counts. If TRUE, rows with zero counts or person years are removed.
#' @return A data frame with processed cohort diagnostics incidence rate data.
#' @export
processCohortDiagnosticsIncidenceRateData <- function(cohortDiagnosticsIncidenceRateData,
                                                      zeroCountAdjustment) {
  if (!zeroCountAdjustment) {
    cohortDiagnosticsIncidenceRateData <- cohortDiagnosticsIncidenceRateData |>
      dplyr::filter(cohortCount > 0) |>
      dplyr::filter(personYears > 0)
  }
  
  outputData <- cohortDiagnosticsIncidenceRateData |>
    dplyr::filter(calendarYear != "") |>
    dplyr::filter(ageGroup == '') |> # we are not using age stratified incidence rate for the diagnostic
    dplyr::filter(gender == '') |> # we are not using gender stratified incidence rate for the diagnostic
    dplyr::mutate(
      cohortCount = abs(cohortCount),
      personYears = abs(personYears),
      incidenceRate = abs(cohortCount) / abs(personYears),
      calendarYear = as.integer(calendarYear)
    ) |>
    dplyr::select(
      cohortId,
      databaseId,
      gender,
      ageGroup,
      calendarYear,
      cohortCount,
      personYears,
      incidenceRate
    ) |>
    imputeMissingYears() |>  # we do not need to impute data in general for the Cohort Diagnostics output of incidence rate because it has 0 values
    dplyr::mutate(calendarYear = as.Date(paste0(calendarYear, "-06-01"))) # mid year
  
  return(outputData)
}


#' Check Temporal Stability For Cohort Diagnostics Incidence Rate Data
#'
#' This function checks the temporal stability of the cohort diagnostics incidence rate data. It processes the data, calculates the predicted counts, and compares the likelihoods. The function also adds a flag to indicate whether the data was evaluated and whether it is stable.
#'
#' @param cohortDiagnosticsIncidenceRateData A data frame that contains the columns cohortId, databaseId, gender, ageGroup, calendarYear, cohortCount, personYears, and incidenceRate.
#' @param cohort A data frame that contains the columns cohortId and cohortName.
#' @param includeZeroCount A logical value indicating whether to include zero counts in the analysis. If TRUE, rows with zero counts are included.
#' @param maxNumberOfSplines The maximum number of splines to use in the model. If NULL, the default is 5.
#' @param splineTickInterval The interval at which to place the splines. The default is 3.
#' @return A data frame with the cohortId, databaseId, gender, ageGroup, evaluated flag, stable flag, and isUnstable flag.
#' @export
checkTemporalStabilityForcohortDiagnosticsIncidenceRateData <- function(cohortDiagnosticsIncidenceRateData,
                                                                        cohort,
                                                                        includeZeroCount = TRUE,
                                                                        maxNumberOfSplines = NULL,
                                                                        splineTickInterval = 3) {
  zeroCountAdjustment <- includeZeroCount
  
  checkIfCohortDiagnosticsIncidenceRateData(cohortDiagnosticsIncidenceRateData)
  
  observedCount <- processCohortDiagnosticsIncidenceRateData(cohortDiagnosticsIncidenceRateData = cohortDiagnosticsIncidenceRateData,
                                                             zeroCountAdjustment = includeZeroCount)
  
  combos <- observedCount |>
    dplyr::select(cohortId, databaseId, gender, ageGroup) |>
    dplyr::distinct() |>
    dplyr::mutate(rn = dplyr::row_number())
  
  output <- c()
  
  # Create a progress bar
  pb <- txtProgressBar(min = 0,
                       max = nrow(combos),
                       style = 3)
  
  for (i in (1:nrow(combos))) {
    combo <- combos[i, ]
    
    data <- observedCount |>
      dplyr::inner_join(combo, by = c("cohortId", "databaseId", "gender", "ageGroup"))
    
    output[[i]] <- dplyr::tibble(ratio = 0,
                                 p = 0,
                                 evaluated = FALSE) |>
      tidyr::crossing(combo)
    
    if (nrow(data) > 1) {
      predicted <- getPredictedCount(
        data = data,
        timeSequenceField = "calendarYear",
        countField = "cohortCount",
        personTimeField = "personYears",
        zeroCountAdjustment = zeroCountAdjustment,
        maxNumberOfSplines = NULL,
        splineTickInterval = 3
      )
      
      if (nrow(data) == nrow(predicted)) {
        output[[i]] <- predicted |>
          likelihoodComparison(zeroCountAdjustment = zeroCountAdjustment) |>
          tidyr::crossing(combo) |>
          dplyr::mutate(evaluated = TRUE)
      }
    }
    
    # Update the progress bar
    setTxtProgressBar(pb, i)
  }
  
  # Close the progress bar
  close(pb)
  
  output <- dplyr::bind_rows(output) |>
    dplyr::select(cohortId, databaseId, gender, ageGroup, evaluated, stable) |>
    dplyr::inner_join(cohort |>
                        dplyr::select(cohortId, cohortName), by = "cohortId") |>
    dplyr::relocate(cohortId, cohortName)
  
  return(output)
}



#' Plot Simple Temporal Trend
#'
#' This function generates a line plot of incidence rates over time, grouped by database ID.
#'
#' @param data A data frame containing the data to be plotted.
#' @param xAxisCol A string specifying the column in `data` to use for the x-axis. Default is "calendarYear".
#' @param yAxisCol A string specifying the column in `data` to use for the y-axis. Default is "incidenceRate".
#' @param groupBy A string specifying the column in `data` to use for grouping the data. Default is "databaseId".
#' @param plotTitle A string specifying the title of the plot. Default is "Incidence Rate over time".
#' @param xLabel A string specifying the label for the x-axis. Default is "Calendar year".
#' @param yLabel A string specifying the label for the y-axis. Default is "Incidence Rate".
#' @param useMinimalTheme A logical value indicating whether to use ggplot2's minimal theme. Default is TRUE.
#' @param convertToPlotLy A logical value indicating whether to convert the ggplot object to a plotly object for interactive visualization. Default is FALSE.
#' @return A ggplot or plotly object, depending on the value of `convertToPlotLy`.
#' @export
plotSimpleTemporalTrend <- function(data,
                                    xAxisCol = "calendarYear",
                                    yAxisCol = "incidenceRate",
                                    groupBy = "databaseId",
                                    plotTitle = "Incidence Rate over time",
                                    xLabel = "Calendar year",
                                    yLabel = "Incidence Rate",
                                    useMinimalTheme = TRUE,
                                    convertToPlotLy = FALSE) {
  # Generate color palette based on unique database IDs
  numColors <- length(unique(data[[groupBy]]))
  colors <- OhdsiRPlots::createOhdsiPalette(numColors = numColors)
  colorMapping <- setNames(colors, unique(data[[groupBy]]))
  
  p <-
    ggplot2::ggplot(data, ggplot2::aes(x = .data[[xAxisCol]], y = .data[[yAxisCol]], color = .data[[groupBy]])) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::scale_color_manual(values = colorMapping) +
    ggplot2::labs(
      title = plotTitle,
      x = xLabel,
      y = yLabel,
      color = "Database ID"
    )
  
  if (useMinimalTheme) {
    p <- p + ggplot2::theme_minimal()
  }
  
  if (convertToPlotLy) {
    # Convert ggplot object to plotly for interactive visualization
    p_plotly <-
      plotly::ggplotly(p, tooltip = c("x", "y", groupBy)) # Set tooltip to show data from specific columns
    return(p_plotly)
  } else {
    return(p)
  }
}