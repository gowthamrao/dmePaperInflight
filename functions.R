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
#' @return A data frame with the observed and expected counts.
#' @export
getPredictedCount <- function(data,
                              timeSequenceField,
                              countField,
                              personTimeField,
                              maxNumberOfSplines = NULL,
                              splineTickInterval = 3) {
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
  
  # Set the default value for maxNumberOfSplines if it's NULL
  if (is.null(maxNumberOfSplines)) {
    maxNumberOfSplines <- 5
  }
  
  numberOfXAxisTicks <- length(unique(data$timeId))
  
  # Determine the number of splines based on the number of unique timeId values
  numberOfSplines <-
    min(maxNumberOfSplines,
        ceiling(numberOfXAxisTicks / splineTickInterval)) # we will put a spline for every 3 ticks
  
  # Create a Cyclops data object using the observed counts, timeId, and personTimeField
  # The formula used here is a Poisson regression model with natural splines
  offsetTerm <- paste0("offset(log(", personTimeField, "))")
  
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
  data$expected <- as.double(predictions)
  
  data$numberOfSplinesUsed <- numberOfSplines
  
  # Arrange the data based on timeId and select the 'observed' and 'expected' columns
  data <- data |>
    dplyr::arrange(timeId) |>
    dplyr::select(-timeId) |>
    dplyr::tibble()
  
  return(data)
}




#' Compare Likelihoods in a Data Frame
#'
#' This function compares the likelihoods of observed and expected counts in a data frame. It uses a Poisson regression model with natural splines and calculates the log-likelihood ratio. The function also adds a flag to indicate whether the data was imputed.
#'
#' @param data A data frame that contains the columns 'observed' and 'expected'.
#' @param maxRatio The maximum ratio for the likelihood comparison. If NULL, the default is 1.25.
#' @param alpha The significance level for the likelihood ratio test. If NULL, the default is 0.05.
#' @return A data frame with the ratio, p-value, and stability indicator.
#' @export
likelihoodComparison <- function(data,
                                 maxRatio = 1.25,
                                 alpha = 0.05) {
  observed <- data$observed
  expected <- data$expected
  
  # From https://cdsmithus.medium.com/the-logarithm-of-a-sum-69dd76199790
  smoothMax <- function(x, y) {
    smoothMax <- ifelse(
      test = abs(x - y) > 100,
      yes = pmax(x, y),
      no = x + log(1 + exp(y - x))
    )
    return(smoothMax)
  }
  
  logLikelihood <- function(x) {
    xTerm <- dpois(observed, expected * x, log = TRUE)
    yTerm <- dpois(observed, expected / x, log = TRUE)
    return(-sum(smoothMax(x = xTerm, y = yTerm)))
  }
  
  likelihood <- function(x) {
    return(exp(-logLikelihood(x)))
  }
  
  vectorLikelihood <- function(x) {
    return(sapply(x, likelihood))
  }
  
  x <- seq(1, 10, by = 0.1)
  ll <- sapply(x, logLikelihood)
  
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
    method = "L-BFGS-B"
  )$par
  
  l0 <- integrate(vectorLikelihood, lower = 1, upper = maxRatio)$value
  l1 <- integrate(vectorLikelihood, lower = maxRatio, upper = Inf)$value
  
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
  result <- data |>
    tidyr::crossing(dplyr::tibble(
      ratio = xHat,
      p = p,
      stable = p > alpha
    ))
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
      !is.character(data$calendarYear) |
      !is.character(data$databaseId)) {
    stop("Data types of some columns are not as expected.")
  }
  
  if (data |>
      dplyr::select(databaseId, calendarYear, cohortId) |>
      nrow() != data |>
      dplyr::select(databaseId, calendarYear, cohortId) |>
      dplyr::distinct() |>
      nrow()) {
    stop("data is not unique by databaseId, calendarYear, cohortId")
  }
}


#' Process Cohort Diagnostics Incidence Rate Data
#'
#' This function processes the cohort diagnostics incidence rate data. It filters the data based on certain conditions, calculates the incidence rate, and imputes missing years.
#'
#' @param cohortDiagnosticsIncidenceRateData A data frame that contains the columns cohortId, databaseId, gender, ageGroup, calendarYear, cohortCount, personYears, and incidenceRate.
#' @return A data frame with processed cohort diagnostics incidence rate data.
#' @export
processCohortDiagnosticsIncidenceRateData <- function(cohortDiagnosticsIncidenceRateData) {
  outputData <- cohortDiagnosticsIncidenceRateData |>
    dplyr::filter(personYears > 0) |> #remove any records with 0 or lower (below min cell count for privacy protection) person years
    dplyr::mutate(
      # to adjust for min cell count privacy protection output in CohortDiagnostics
      cohortCount = abs(cohortCount),
      personYears = abs(personYears),
      incidenceRate = abs(incidenceRate),
      calendarYear = as.integer(calendarYear)
    ) |>
    dplyr::select(cohortId,
                  databaseId,
                  calendarYear,
                  cohortCount,
                  personYears,
                  incidenceRate) |>
    dplyr::mutate(calendarYear = as.Date(paste0(calendarYear, "-06-01"))) # mid year
  
  return(outputData)
}


#' Check Temporal Stability For Cohort Diagnostics Incidence Rate Data
#'
#' This function checks the temporal stability of the cohort diagnostics incidence rate data. It processes the data, calculates the predicted counts, and compares the likelihoods.
#'
#' @param cohortDiagnosticsIncidenceRateData A data frame that contains the columns cohortId, databaseId, gender, ageGroup, calendarYear, cohortCount, personYears, and incidenceRate.
#' @param cohort A data frame that contains the columns cohortId and cohortName.
#' @param maxNumberOfSplines The maximum number of splines to use in the model. If NULL, the default is 5.
#' @param splineTickInterval The interval at which to place the splines. The default is 3.
#' @return A data frame with the cohortId, databaseId, gender, ageGroup, stable flag, and isUnstable flag.
#' @param maxRatio The maximum ratio for the likelihood comparison. If NULL, the default is 1.25.
#' @param alpha The significance level for the likelihood ratio test. If NULL, the default is 0.05.
#' @export
checkTemporalStabilityForcohortDiagnosticsIncidenceRateData <- function(cohortDiagnosticsIncidenceRateData,
                                                                        cohort,
                                                                        maxNumberOfSplines = NULL,
                                                                        splineTickInterval = 3,
                                                                        maxRatio = 1.25,
                                                                        alpha = 0.05) {
  checkIfCohortDiagnosticsIncidenceRateData(cohortDiagnosticsIncidenceRateData)
  
  observedCount <- processCohortDiagnosticsIncidenceRateData(cohortDiagnosticsIncidenceRateData = cohortDiagnosticsIncidenceRateData)
  
  combos <- observedCount |>
    dplyr::select(cohortId, databaseId) |>
    dplyr::distinct()
  
  output <- c()
  
  # Create a progress bar
  pb <- txtProgressBar(min = 0,
                       max = nrow(combos),
                       style = 3)
  
  for (i in (1:nrow(combos))) {
    combo <- combos[i, ]
    
    data <- observedCount |>
      dplyr::inner_join(combo, by = c("cohortId", "databaseId"))
    
    if (nrow(data) > 1) {
      predicted <- getPredictedCount(
        data = data,
        timeSequenceField = "calendarYear",
        countField = "cohortCount",
        personTimeField = "personYears",
        maxNumberOfSplines = maxNumberOfSplines,
        splineTickInterval = 3
      )
      
      predicted <- predicted |>
        dplyr::mutate(
          expectedIncidenceRate = round(x = (expected / personYears), digits = 10),
          expected = round(x = expected, digits = 4),
          observedExpectedCountRatio = round(x = observed / expected, digits = 4),
          observedExpectedIncidenceRateRatio = round(x = incidenceRate / expectedIncidenceRate, digits = 4),
          percentageDeviation = (observed - expected) / expected * 100
        )
      
      likelihood <- predicted |>
        likelihoodComparison(maxRatio = maxRatio, alpha = alpha)
      
      output[[i]] <- data |>
        dplyr::left_join(likelihood, by = colnames(data))
    }
    
    # Update the progress bar
    setTxtProgressBar(pb, i)
  }
  
  # Close the progress bar
  close(pb)
  
  output <- dplyr::bind_rows(output) |>
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
                                    xLabel = "Calendar year",
                                    yLabel = "Incidence Rate",
                                    useMinimalTheme = TRUE,
                                    convertToPlotLy = FALSE) {
  # Generate color palette based on unique database IDs
  numColors <- length(unique(data[[groupBy]]))
  colors <- OhdsiRPlots::createOhdsiPalette(numColors = numColors)
  colorMapping <- setNames(colors, unique(data[[groupBy]]))
  
  plotTitle <- paste0("Observed values only\n")
  
  p <- ggplot2::ggplot(data,
                       ggplot2::aes(
                         x = .data[[xAxisCol]],
                         y = .data[[yAxisCol]],
                         color = .data[[groupBy]],
                         group = .data[[groupBy]]
                       )) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::scale_color_manual(values = colorMapping) +
    ggplot2::scale_linetype_manual(values = c("solid")) +
    ggplot2::scale_y_continuous(limits = c(0, NA)) + # Ensure y-axis starts at 0
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





plotTemporalTrendExpectedObserved <- function(data,
                                              xAxisCol = "calendarYear",
                                              yAxisCol = c("expectedIncidenceRate", "incidenceRate"),
                                              groupBy = "databaseId",
                                              xLabel = "Calendar year",
                                              yLabel = "Incidence Rate",
                                              useMinimalTheme = TRUE,
                                              convertToPlotLy = FALSE) {
  # Generate color palette based on unique database IDs
  numColors <- length(unique(data[[groupBy]]))
  colors <- OhdsiRPlots::createOhdsiPalette(numColors = numColors)
  colorMapping <- setNames(colors, unique(data[[groupBy]]))
  
  numberOfSplinesUsed <- 0
  
  if ('numberOfSplinesUsed' %in% colnames(data)) {
    numberOfSplinesUsed <- max(data$numberOfSplinesUsed)
  }
  
  plotTitle <-
    paste0(
      "\nObserved vs Expected. Splines used = ",
      numberOfSplinesUsed,
      ". p-value: ",
      data$p |> min() |> OhdsiHelpers::formatDecimalWithComma(decimalPlaces = 4),
      ". ratio = ",
      data$ratio |> min() |> OhdsiHelpers::formatDecimalWithComma(decimalPlaces = 4)
    )
  
  
  plotTitle <- if (min(data$stable) == 0) {
    plotTitle <- paste0(plotTitle, "\nUNSTABLE (FAILED Diagnostic)")
  } else {
    plotTitle <- paste0(plotTitle, "\nSTABLE (Did not fail diagnostic)")
  }
  
  p <- ggplot2::ggplot(data) +
    ggplot2::geom_line(
      ggplot2::aes(
        x = .data[[xAxisCol]],
        y = .data[[yAxisCol[1]]],
        color = .data[[groupBy]],
        group = .data[[groupBy]],
        linetype = "Expected"
      )
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        x = .data[[xAxisCol]],
        y = .data[[yAxisCol[2]]],
        color = .data[[groupBy]],
        group = .data[[groupBy]],
        linetype = "Observed"
      )
    ) +
    ggplot2::scale_color_manual(values = colorMapping) +
    ggplot2::scale_linetype_manual(values = c("Expected" = "dashed", "Observed" = "solid")) +
    ggplot2::scale_y_continuous(limits = c(0, NA)) + # Ensure y-axis starts at 0
    ggplot2::labs(
      title = plotTitle,
      x = xLabel,
      y = yLabel,
      color = "Database ID",
      linetype = "Incidence Type"
    )
  
  if (useMinimalTheme) {
    p <- p + ggplot2::theme_minimal()
  }
  
  if (convertToPlotLy) {
    # Convert ggplot object to plotly for interactive visualization
    p_plotly <-
      plotly::ggplotly(p, tooltip = c("x", "y", groupBy, "linetype")) # Set tooltip to show data from specific columns
    return(p_plotly)
  } else {
    return(p)
  }
}



createFakeIncidenceRateData <- function(numberOfYears,
                                        cohortId,
                                        databaseId,
                                        pattern = "random",
                                        incidenceRateMin = 0.001,
                                        # Set a positive minimum to avoid division by zero
                                        incidenceRateMax = 0.02,
                                        cohortCountMin = 10,
                                        # Ensure a minimum positive count
                                        cohortCountMax = 10000,
                                        errorMagnitude = 0.01) {
  # Percentage of random error
  library(dplyr)
  library(tibble)
  
  # Generate base years and cohort counts
  calendarYear <- seq(from = 2024,
                      length.out = numberOfYears,
                      by = -1)
  cohortCount <- round(runif(numberOfYears, min = cohortCountMin, max = cohortCountMax))
  
  # Pattern-specific incidence rate generation
  incidenceRate <- numeric(numberOfYears)
  switch(
    pattern,
    stable = {
      rate <- runif(1, incidenceRateMin, incidenceRateMax)
      incidenceRate <- rep(rate, numberOfYears)
    },
    monotonic = {
      incidenceRate <- seq(incidenceRateMin, incidenceRateMax, length.out = numberOfYears)
    },
    random = {
      incidenceRate <- runif(numberOfYears, incidenceRateMin, incidenceRateMax)
    },
    changePoints = {
      changePoint <- sample(2:(numberOfYears - 1), 1)
      incidenceRate[1:changePoint] <- runif(1, incidenceRateMin, incidenceRateMax)
      incidenceRate[(changePoint + 1):numberOfYears] <- runif(1, incidenceRateMin, incidenceRateMax)
    },
    exponential = {
      incidenceRate <- incidenceRateMin * exp(seq(
        log(1),
        log(1 + incidenceRateMax - incidenceRateMin),
        length.out = numberOfYears
      ))
    },
    stepwise = {
      stepChange <- round(numberOfYears / 2)
      incidenceRate[1:stepChange] <- incidenceRateMin
      incidenceRate[(stepChange + 1):numberOfYears] <- incidenceRateMax
    },
    cyclic = {
      period <- 4
      for (i in 1:numberOfYears) {
        incidenceRate[i] <- incidenceRateMin + (incidenceRateMax - incidenceRateMin) * abs(sin(2 * pi * i / period))
      }
    },
    autoregressive = {
      incidenceRate[1] <- runif(1, incidenceRateMin, incidenceRateMax)
      for (i in 2:numberOfYears) {
        incidenceRate[i] <- max(incidenceRateMin,
                                0.5 * incidenceRate[i - 1] + runif(1, -0.01, 0.01))
      }
    },
    randomWalk = {
      incidenceRate[1] <- runif(1, incidenceRateMin, incidenceRateMax)
      for (i in 2:numberOfYears) {
        incidenceRate[i] <- incidenceRate[i - 1] + runif(1, -0.01, 0.01)
        incidenceRate[i] <- max(incidenceRateMin,
                                min(incidenceRate[i], incidenceRateMax))
      }
    }
  )
  
  # Add random error to each pattern
  error <- runif(numberOfYears, -errorMagnitude, errorMagnitude) * incidenceRate
  incidenceRate <- pmax(incidenceRate + error, 0.0001)  # Apply error and ensure positivity
  
  # Calculate person years and correct if zero to prevent division by zero
  personYears <- round(cohortCount / incidenceRate)
  personYears <- ifelse(personYears == 0, 1, personYears)  # Ensure personYears is not zero
  
  # Generate the tibble
  output <- tibble(
    calendarYear = as.character(calendarYear),
    cohortId = cohortId,
    databaseId = databaseId,
    incidenceRate = incidenceRate,
    cohortCount = cohortCount,
    personYears = personYears
  )
  
  return(output)
}
