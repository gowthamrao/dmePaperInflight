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
  
  numberOfXAxisTicks <- length(unique(data$timeId))
  
  # Determine the number of splines based on the number of unique timeId values
  numberOfSplines <-
    min(maxNumberOfSplines,
        ceiling(numberOfXAxisTicks / splineTickInterval)) # we will put a spline for every 3 ticks
  
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
  data$expected <- as.double(predictions)
  
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
#' @return A data frame with processed cohort diagnostics incidence rate data.
#' @export
processCohortDiagnosticsIncidenceRateData <- function(cohortDiagnosticsIncidenceRateData) {
  outputData <- cohortDiagnosticsIncidenceRateData |>
    dplyr::filter(personYears > 0) |> #remove any records with 0 or lower (below min cell count for privacy protection) person years
    dplyr::filter(calendarYear != "") |>
    dplyr::filter(ageGroup == '') |> # we are not using age stratified incidence rate for the diagnostic
    dplyr::filter(gender == '') |> # we are not using gender stratified incidence rate for the diagnostic
    dplyr::mutate(
      # to adjust for min cell count privacy protection output in CohortDiagnostics
      cohortCount = abs(cohortCount),
      personYears = abs(personYears),
      incidenceRate = abs(incidenceRate),
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
    dplyr::mutate(calendarYear = as.Date(paste0(calendarYear, "-06-01"))) # mid year
  
  return(outputData)
}


#' Check Temporal Stability For Cohort Diagnostics Incidence Rate Data
#'
#' This function checks the temporal stability of the cohort diagnostics incidence rate data. It processes the data, calculates the predicted counts, and compares the likelihoods. The function also adds a flag to indicate whether the data was evaluated and whether it is stable.
#'
#' @param cohortDiagnosticsIncidenceRateData A data frame that contains the columns cohortId, databaseId, gender, ageGroup, calendarYear, cohortCount, personYears, and incidenceRate.
#' @param cohort A data frame that contains the columns cohortId and cohortName.
#' @param filterZeroAndLowCountRecords A logical value indicating whether to remove records from the observed data that have zero persons in the cohort during the calendar period or records that fall below the privacy protecting minimum threshold specified during the cohort diagnostics run. If set to TRUE, the function will exclude these records from the analysis. Default is FALSE.
#' @param maxNumberOfSplines The maximum number of splines to use in the model. If NULL, the default is 5.
#' @param splineTickInterval The interval at which to place the splines. The default is 3.
#' @return A data frame with the cohortId, databaseId, gender, ageGroup, evaluated flag, stable flag, and isUnstable flag.
#' @export
checkTemporalStabilityForcohortDiagnosticsIncidenceRateData <- function(cohortDiagnosticsIncidenceRateData,
                                                                        cohort,
                                                                        filterZeroAndLowCountRecords = FALSE,
                                                                        maxNumberOfSplines = NULL,
                                                                        splineTickInterval = 3) {
  zeroCountAdjustment <- TRUE
  checkIfCohortDiagnosticsIncidenceRateData(cohortDiagnosticsIncidenceRateData)
  
  # If filterZeroAndLowCountRecords is TRUE, remove records with zero persons in the cohort or records below the minimum threshold
  if (filterZeroAndLowCountRecords) {
    cohortDiagnosticsIncidenceRateData <- cohortDiagnosticsIncidenceRateData |>
      dplyr::filter(cohortCount > 0) # this also removes records that have negative values i.e. below threshold
    writeLines(
      "Cohort periods with zero counts or those below the privacy protecting minimum threshold value are removed from analysis."
    )
  }
  
  observedCount <- processCohortDiagnosticsIncidenceRateData(cohortDiagnosticsIncidenceRateData = cohortDiagnosticsIncidenceRateData)
  
  combos <- observedCount |>
    dplyr::select(cohortId, databaseId, gender, ageGroup) |>
    dplyr::distinct()
  
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
      
      predicted <- predicted |>
        dplyr::mutate(
          expectedIncidenceRate = round(x = (expected / personYears) * 1000, digits = 4),
          expected = round(x = expected, digits = 4),
          observedExpectedCountRatio = round(x = observed / expected, digits = 4),
          observedExpectedIncidenceRateRatio = round(x = incidenceRate / expectedIncidenceRate, digits = 4),
          percentageDeviation = (observed - expected) / expected * 100,
          deviationDirection = dplyr::if_else(observed < expected, "Under", "Equal or Over")
        )
      
      if (nrow(data) == nrow(predicted)) {
        output[[i]] <- predicted |>
          likelihoodComparison(zeroCountAdjustment = zeroCountAdjustment) |>
          dplyr::mutate(evaluated = TRUE)
      }
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


createTemporalStabilityPlot <- function(data) {
  # Format the title with unique data values
  title <- paste0(
    paste("Database Id:", unique(data$databaseId)),
    "\n",
    paste("Cohort: ", unique(data$cohortName), " (", data$cohortId |> unique(), ")"),
    "\n",
    paste(
      "Overall ratio: ",
      unique(data$ratio) |> OhdsiHelpers::formatDecimalWithComma(decimalPlaces = 2),
      ", p-value",
      unique(data$p) |> OhdsiHelpers::formatDecimalWithComma(decimalPlaces = 3)
    )
    # ,
    # "\n",
    # paste("Stable: ", as.logical(unique(data$stable)))
  )
  
  bgColor = "white"
  
  color <- ifelse(max(data$stable) == 1, 'white', 'red')
  
  # Define a minimal theme with the custom background color and no grid lines
  custom_theme <- function(bgColor) {
    theme_minimal() +
      theme(
        plot.background = element_rect(fill = bgColor, color = bgColor),
        panel.background = element_rect(fill = bgColor, color = bgColor),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = 10, r = 10, b = 0, l = 10)
      )
  }
  
  # Create the upper plot
  upper_plot <- ggplot(data, aes(x = calendarYear)) +
    geom_line(aes(y = incidenceRate), color = "blue") +
    geom_line(aes(y = expectedIncidenceRate), color = "gray", linetype = "dotted") +
    labs(y = "Incidence Rate", x = NULL) +
    custom_theme(color)
  
  # Create the lower plot
  lower_plot <- ggplot(data, aes(x = calendarYear, y = abs(percentageDeviation))) +
    geom_bar(stat = "identity", aes(fill = deviationDirection)) +
    scale_fill_manual(values = c("Under" = "orange", "Over" = "gray")) +
    scale_y_continuous(limits = c(0, 100), labels = scales::percent_format(scale = 1)) +
    labs(y = "Abs Percentage Deviation", x = "Calendar Year") +
    custom_theme(color) +
    theme(legend.position = "none")
  
  # Combine the plots
  final_plot <- invisible(
    gridExtra::grid.arrange(
      upper_plot,
      lower_plot,
      ncol = 1,
      heights = c(4, 1),
      top = grid::textGrob(
        x = 0,
        y = 0,
        title,
        gp = grid::gpar(fontsize = 14, fontface = "bold"),
        just = "left",
        vjust = 1
      )
    )
  )
  
  # Return the final plot object
  return(final_plot)
}




# Wrapper function that processes the data and generates the trellis of plots
createTrellisOfPlots <- function(data) {
  library(ggplot2)
  library(gridExtra)
  library(dplyr)
  # Unique combinations of cohortName and databaseName
  combinations <- unique(data[c("cohortName", "databaseName")])
  
  # List to store plots
  plot_list <- vector("list", length(combinations))
  
  # Generate a plot for each combination
  for (i in seq_along(combinations$cohortName)) {
    subset_data <- data %>% 
      filter(cohortName == combinations$cohortName[i], databaseName == combinations$databaseName[i])
    
    if (nrow(subset_data) > 0) {
      plot_list[[i]] <- createTemporalStabilityPlot(subset_data)
    } else {
      plot_list[[i]] <- NULL # Leave empty if no data
    }
  }
  
  # Grid dimensions
  n_cols <- 3
  n_rows <- ceiling(length(plot_list) / n_cols)
  
  # Create a trellis layout of the plots
  plot_grid <- gridExtra::grid.arrange(
    grobs = plot_list,
    ncol = n_cols,
    nrow = n_rows,
    top = "Trellis of Temporal Stability Plots"
  )
  
  return(plot_grid)
}

