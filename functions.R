
source("performGAMAnalysis.R")

# Reference population proportions (should come from an external source)
referencePopulation <- data.frame(
  ageGroup = c(
    '0-9',
    '10-19',
    '20-29',
    '30-39',
    '40-49',
    '50-59',
    '60-69',
    '70-79',
    '80-89',
    '90-99',
    '100-109'
  ),
  gender = rep(c('Female', 'Male'), each = 11),
  proportion = c(
    0.030,
    0.032,
    0.035,
    0.035,
    0.031,
    0.035,
    0.033,
    0.025,
    0.011,
    0.004,
    0.002,
    0.030,
    0.032,
    0.035,
    0.034,
    0.031,
    0.035,
    0.033,
    0.024,
    0.010,
    0.002,
    0.001
  )
) |> dplyr::tibble()

drawPopulationPyramid <- function(data) {
  # Ensure ggplot2 is installed
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install the ggplot2 package to use this function.")
  }
  
  # Load ggplot2
  library(ggplot2)
  
  # Extract the start of the age range for sorting
  data$ageStart <- as.numeric(sub("-.*", "", data$ageGroup))
  
  # Order the data by the numeric start of the age group
  data <- data[order(data$ageStart), ]
  
  # Adjust proportions for plotting (make male proportions negative)
  data$proportion <- ifelse(data$gender == 'Male', -data$proportion, data$proportion)
  
  # Create the plot
  p <- ggplot(data, aes(x = reorder(ageGroup, ageStart), y = proportion, fill = gender)) +
    geom_bar(stat = "identity") +
    coord_flip() +  # Flip coordinates for pyramid style
    scale_y_continuous(labels = abs) +  # Show positive y-axis labels
    labs(title = "Population Pyramid", x = "Age Group", y = "Proportion of Population") +
    scale_fill_manual(values = c("Male" = "blue", "Female" = "pink")) +
    theme_minimal()
  
  # Display the plot
  return(p)
}

# Example usage with the provided reference_population data frame
drawPopulationPyramid(referencePopulation)





#' Get Predicted Count in a Data Frame
#'
#' This function predicts the count of events in a data frame using a Poisson regression model with natural splines. It groups the data by timeSequenceField and fills in missing counts with predicted values.
#' The function also adds a flag to indicate whether the data was imputed.
#'
#' @param data A data frame that contains the columns specified in the function parameters.
#' @param timeSequenceField The field in the data frame that represents the time sequence.
#' @param countField The field in the data frame that represents the count of events.
#' @param personTimeField The field in the data frame that represents the person time.
#' @param maxNumberOfSplines The maximum number of splines to use in the model. If NULL, the default is 3.
#' @param splineTickInterval The interval at which to place the splines. The default is 3.
#' @param maxRatio The maximum ratio for the likelihood comparison. If NULL, the default is 1.25.
#' @param alpha The significance level for the likelihood ratio test. If NULL, the default is 0.05.
#' @return A data frame with the observed and expected counts.
#' @export
getPredictedCount <- function(data,
                              timeSequenceField,
                              countField,
                              personTimeField,
                              maxNumberOfSplines = NULL,
                              splineTickInterval = 3,
                              alpha = 0.05,
                              maxRatio = 1.25) {
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
    maxNumberOfSplines <- 3
  }
  
  numberOfXAxisTicks <- length(unique(data$timeId))
  
  # Determine the number of splines based on the number of unique timeId values
  numberOfSplines <-
    min(maxNumberOfSplines,
        ceiling(numberOfXAxisTicks / splineTickInterval)) # we will put a spline for every 3 ticks
  
  data$numberOfSplinesUsed <- numberOfSplines
  
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
  data$cyclopsExpected <- as.double(predictions)
  
  data <- data |>
    likelihoodComparison(maxRatio = maxRatio, alpha = alpha) |>
    dplyr::rename(
      cyclopsRatio = ratio,
      cyclopsPValue = p,
      cyclopsStable = stable
    )
  
  # We cannot get confidence interval because Cyclops does not provide std errors.
  # to get confidence intervals we will have to use another package other than cyclops like Glm
  
  
  # ########
  # # Fit the Poisson regression using glm() with the correct reference to the personTimeField column
  # # Try fitting Poisson regression using glm() and handle errors
  tryCatch({
    # assign default values
    data$glmExpected <- as.double(NA)
    data$glmExpectedLowerBound <- as.double(NA)
    data$glmExpectedUpperBound <- as.double(NA)
    data$glmDevianceValue <- as.double(NA)
    data$glmDegreesOfFreedom <- as.double(NA)
    data$glmPValueDeviance <- as.double(NA)
    data$glmPearsonChiSquare <- as.double(NA)
    data$glmPPValuePearson <- as.double(NA)
    
    modelGlm <- glm(
      observed ~ splines::ns(timeId, df = numberOfSplines) + offset(log(data[[personTimeField]])),
      data = data,
      family = poisson
    )
    glmPredictions <- stats::predict(modelGlm,
                                     newdata = data,
                                     type = "link",
                                     se.fit = TRUE)
    glmPredictedLog <- glmPredictions$fit
    seLog <- glmPredictions$se.fit
    
    zValue <- qnorm(1 - alpha / 2)
    
    glmLowerBoundLog <- glmPredictedLog - zValue * seLog
    glmUpperBoundLog <- glmPredictedLog + zValue * seLog
    
    # Exponentiate to get the predicted counts and CIs on the original scale
    data$glmExpected <- as.double(exp(glmPredictedLog))
    data$glmExpectedLowerBound <- as.double(exp(glmLowerBoundLog))
    data$glmExpectedUpperBound <- as.double(exp(glmUpperBoundLog))
    
    #using Martijn's likelikhood custom function to calculate p-value
    glmLikelihood <- data |>
      dplyr::select(glmExpected, observed) |>
      dplyr::rename(cyclopsExpected = glmExpected) |>
      likelihoodComparison(maxRatio = maxRatio, alpha = alpha) |>
      dplyr::rename(glmRatio = ratio,
                    glmPValue = p,
                    glmStable = stable) |>
      dplyr::select(glmRatio, glmPValue, glmStable) |>
      dplyr::distinct()
    
    data <- data |>
      tidyr::crossing(glmLikelihood)
    
    # These tests, although valid, give different results from Martijn's likelihood function.
    
    # Deviance Test (G-test)
    # In Poisson regression, the deviance measures the difference between the observed and expected counts under the model.
    # Get the deviance from the fitted model
    data$glmDevianceValue <- modelGlm$deviance
    # Get the degrees of freedom (difference between the number of observations and the number of parameters)
    data$glmDegreesOfFreedom <- modelGlm$df.residual
    # Compute the p-value from the chi-square distribution
    data$glmPValueDeviance <- 1 - pchisq(modelGlm$deviance, modelGlm$df.residual)
    
    # Pearson Chi-Squared Test - This test sums the squared differences between the observed and expected counts, scaled by the expected counts.
    # Compute Pearson's chi-squared test statistic
    data$glmPearsonChiSquare <- sum((data$observed - data$glmExpected) ^
                                      2 / data$glmExpected)
    # Compute p-value for the Pearson chi-squared test
    data$glmPPValuePearson <- 1 - pchisq(data$glmPearsonChiSquare, modelGlm$df.residual)
    
    data$glmPearsonStable <- min(data$glmPPValuePearson) > alpha
    data$glmDevianceStable <- min(data$glmPValueDeviance) > alpha
    
  }, error = function(e) {
    # If there's an error, skip the glm part
    # message("\nError in glm fitting: ", e$message)
  })
  
  # browser()
  # report <- data |> dplyr::select(
  #   observed,
  #   cyclopsExpected,
  #   glmExpected,
  #   cyclopsPValue,
  #   glmPValue,
  #   glmExpectedLowerBound,
  #   glmExpectedUpperBound,
  #   glmPPValuePearson,
  #   glmPValueDeviance,
  #   cyclopsStable,
  #   glmStable
  # )
  
  # Arrange the data based on timeId and select the 'observed' and 'expected' columns
  data <- data |>
    dplyr::arrange(timeId) |>
    dplyr::select(-timeId) |>
    dplyr::tibble()
  
  return(data)
}




#' Compare Likelihoods in a Data Frame
#'
#' This function compares the likelihoods of observed and expected counts in a data frame and calculates the log-likelihood ratio.
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
  expected <- data$cyclopsExpected
  
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
#' @param removeOutlierZeroCounts (Default TRUE) Do you want to remove continuous zero counts in the beginning and end in the time series?
#' @param removeOutlierUnstableDataSource (Default TRUE) Do you want to remove the outlier time windows that have dramatically lower denominator counts?
#' @return A data frame with processed cohort diagnostics incidence rate data.
#' @export
processCohortDiagnosticsIncidenceRateData <- function(cohortDiagnosticsIncidenceRateData,
                                                      removeOutlierZeroCounts = TRUE,
                                                      removeOutlierUnstableDataSource = TRUE) {
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
    dplyr::mutate(calendarYear = as.Date(paste0(calendarYear, "-06-01"))) |>  # mid year
    dplyr::arrange(cohortId, databaseId, calendarYear) |>
    dplyr::mutate(zeroRecordLeadRemoved = 0,
                  zeroRecordTrailRemoved = 0)
  
  if (removeOutlierZeroCounts) {
    outputData <- outputData |>
      dplyr::group_by(cohortId, databaseId) |>
      dplyr::arrange(calendarYear) |>
      dplyr::mutate(
        # Identify leading and trailing zeros
        isZero = incidenceRate == 0,
        leadZeros = cumsum(!isZero) == 0,
        # Cumulative sum for leading zeros
        trailZeros = rev(cumsum(rev(!isZero))) == 0  # Reverse cumulative sum for trailing zeros
      ) |>
      dplyr::mutate(
        zeroRecordLeadRemoved = sum(leadZeros),
        # Count lead zeros removed
        zeroRecordTrailRemoved = sum(trailZeros) # Count trail zeros removed
      ) |>
      dplyr::filter(!(leadZeros |
                        trailZeros)) |>   # Remove rows that are leading or trailing zeros
      dplyr::select(-isZero, -leadZeros, -trailZeros) |>  # Clean up extra columns
      dplyr::ungroup()
  }
  
  combos <- outputData |>
    dplyr::select(cohortId, databaseId) |>
    dplyr::distinct()
  
  outputDataByCombo <- c()
  for (i in (1:nrow(combos))) {
    outputDataByCombo[[i]] <- outputData |>
      dplyr::inner_join(combos[i, ], by = c("cohortId", "databaseId"))
    
    outputDataByCombo[[i]] <- checkFirstLastYearPersonYearStability(data = outputDataByCombo[[i]])
  }
  
  outputData <- dplyr::bind_rows(outputDataByCombo) |>
    dplyr::arrange(cohortId, databaseId, calendarYear)
  
  return(outputData)
}



# Define a helper function to check stability for each group
checkFirstLastYearPersonYearStability <- function(data) {
  # Exclude first and last calendarYear for model training
  trainData <- data |>
    dplyr::filter(calendarYear != max(calendarYear),
                  calendarYear != min(calendarYear))
  
  # Fit a linear model using trainData
  if (nrow(trainData) > 1) {
    # Ensure there are enough points to fit a model
    model <- lm(personYears ~ calendarYear, data = trainData)
    
    # Calculate the residual standard error and significance threshold
    residualStandardError <- summary(model)$sigma
    
    #keep a very high significance threshold of 2 SE
    significanceThreshold <- 2 * residualStandardError
    
    # Predict the expected personYears for the last calendar year
    lastYear <- max(data$calendarYear)
    predictedPersonYearsLast <- predict(model, newdata = data.frame(calendarYear = lastYear))
    actualPersonYearsLast <- data |>
      dplyr::filter(calendarYear == lastYear) |>
      dplyr::pull(personYears)
    differenceActualToLast <- actualPersonYearsLast - predictedPersonYearsLast
    stableLast <- abs(differenceActualToLast) > significanceThreshold
    
    # Predict the expected personYears for the first calendar year
    firstYear <- min(data$calendarYear)
    predictedPersonYearsFirst <- predict(model, newdata = data.frame(calendarYear = firstYear))
    actualPersonYearsFirst <- data |>
      dplyr::filter(calendarYear == firstYear) |>
      dplyr::pull(personYears)
    differenceActualToFirst <- actualPersonYearsFirst - predictedPersonYearsFirst
    stableFirst <- abs(differenceActualToFirst) > significanceThreshold
    
    outputTestData <- dplyr::tibble(
      calendarYear = c(lastYear, firstYear),
      useYearData = c(stableLast |> as.integer(), stableFirst |> as.integer())
    )
    
    data <- data |>
      dplyr::left_join(outputTestData, by = c("calendarYear")) |>
      tidyr::replace_na(list(useYearData = 1))
    
    # Return the stability result and difference
    return(data)
  } else {
    return(data |>
             dplyr::mutate(useYearData = as.integer(NA)))
  }
}





#' Check Temporal Stability For Cohort Diagnostics Incidence Rate Data
#'
#' This function checks the temporal stability of the cohort diagnostics incidence rate data. It processes the data, calculates the predicted counts, and compares the likelihoods.
#'
#' @param cohortDiagnosticsIncidenceRateData A data frame that contains the columns cohortId, databaseId, gender, ageGroup, calendarYear, cohortCount, personYears, and incidenceRate.
#' @param cohort A data frame that contains the columns cohortId and cohortName.
#' @param maxNumberOfSplines The maximum number of splines to use in the model. If NULL, the default is 5.
#' @param splineTickInterval The interval at which to place the splines. The default is 3.
#' @param removeOutlierZeroCounts (Default TRUE) Do you want to remove continous zero counts in the beginning and end in the time series?
#' @param maxRatio The maximum ratio for the likelihood comparison. If NULL, the default is 1.25.
#' @param alpha The significance level for the likelihood ratio test. If NULL, the default is 0.05.
#' @return A data frame with the cohortId, databaseId, gender, ageGroup, stable flag, and isUnstable flag.
#' @export
checkTemporalStabilityForcohortDiagnosticsIncidenceRateData <- function(cohortDiagnosticsIncidenceRateData,
                                                                        cohort,
                                                                        maxNumberOfSplines = NULL,
                                                                        splineTickInterval = 3,
                                                                        maxRatio = 1.25,
                                                                        alpha = 0.05,
                                                                        removeOutlierZeroCounts = TRUE) {
  checkIfCohortDiagnosticsIncidenceRateData(cohortDiagnosticsIncidenceRateData)
  
  observedCount <- processCohortDiagnosticsIncidenceRateData(
    cohortDiagnosticsIncidenceRateData = cohortDiagnosticsIncidenceRateData,
    removeOutlierZeroCounts = removeOutlierZeroCounts
  )
  
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
      dplyr::inner_join(combo, by = c("cohortId", "databaseId")) |>
      dplyr::arrange(calendarYear) |>
      dplyr::mutate(
        firstYearToCensor = as.Date(NA),
        firstYearCensoredPersonYears = as.numeric(NA),
        lastYearToCensor = as.Date(NA),
        lastYearCensoredPersonYears = as.numeric(NA)
      )
    
    lastYear <- data |>
      dplyr::filter(useYearData == 0) |>
      dplyr::filter(calendarYear > min(data$calendarYear))
    
    if (nrow(lastYear) > 0) {
      data$lastYearToCensor <- lastYear$calendarYear
      data$lastYearCensoredPersonYears <- lastYear$personYears
    }
    
    firstYear <- data |>
      dplyr::filter(useYearData == 0) |>
      dplyr::filter(calendarYear == min(data$calendarYear))
    
    if (nrow(firstYear) > 0) {
      data$firstYearToCensor <- firstYear$calendarYear
      data$firstYearCensoredPersonYears <- firstYear$personYears
    }
    
    data <- data |>
      dplyr::filter(useYearData == 1)
    
    if (nrow(data) > 1) {
      predicted <- getPredictedCount(
        data = data,
        timeSequenceField = "calendarYear",
        countField = "cohortCount",
        personTimeField = "personYears",
        maxNumberOfSplines = maxNumberOfSplines,
        splineTickInterval = splineTickInterval,
        maxRatio = maxRatio,
        alpha = alpha
      )
      
      output[[i]] <- predicted |>
        dplyr::mutate(
          cyclopsExpectedIncidenceRate = round(
            x = (cyclopsExpected / personYears) * 1000,
            digits = 10
          ),
          cyclopsObservedExpectedCountRatio = round(x = observed / cyclopsExpected, digits = 4),
          cyclopsObservedExpectedIncidenceRateRatio = round(
            x = incidenceRate / cyclopsExpectedIncidenceRate,
            digits = 4
          ),
          cyclopsPercentageDeviation = (observed - cyclopsExpected) / cyclopsExpected * 100,
          cylopsExpected = round(x = cyclopsExpected, digits = 4),
          
          ####
          glmExpectedIncidenceRate = round(
            x = (glmExpected / personYears) * 1000,
            digits = 10
          ),
          glmExpectedIncidenceRateUpperBound = round(
            x = (glmExpectedUpperBound / personYears) * 1000,
            digits = 10
          ),
          glmExpectedIncidenceRateLowerBound = round(
            x = (glmExpectedLowerBound / personYears) * 1000,
            digits = 10
          ),
          glmObservedExpectedCountRatio = round(x = observed / glmExpected, digits = 4),
          glmObservedExpectedIncidenceRateRatio = round(
            x = incidenceRate / glmExpectedIncidenceRate,
            digits = 4
          ),
          glmPercentageDeviation = (observed - glmExpected) / glmExpected * 100,
          glmExpected = round(x = glmExpected, digits = 4)
        )
    }
    
    # Update the progress bar
    setTxtProgressBar(pb, i)
  }
  
  # Close the progress bar
  close(pb)
  
  output <- dplyr::bind_rows(output)
  
  commonFields <- intersect(colnames(observedCount), colnames(output))
  
  output <- observedCount |>
    dplyr::left_join(output, by = commonFields) |>
    dplyr::inner_join(cohort |>
                        dplyr::select(cohortId, cohortName), by = "cohortId") |>
    dplyr::relocate(cohortId, cohortName)
  
  return(output)
}



#' Plot Simple Temporal Trend
#'
#' This function generates a line plot of incidence rates over time.
#'
#' @param data A data frame containing the data to be plotted.
#' @param yAxisCol A string specifying the column in `data` to use for the y-axis. Default is "incidenceRate".
#' @param plotTitle A string specifying the title of the plot. Default is "Incidence Rate over time".
#' @param convertToPlotLy A logical value indicating whether to convert the ggplot object to a plotly object for interactive visualization. Default is FALSE.
#' @return A ggplot or plotly object, depending on the value of `convertToPlotLy`.
#' @export
plotSimpleTemporalTrend <- function(data,
                                    yAxisCol = "incidenceRate",
                                    xLabel = "Calendar year",
                                    convertToPlotLy = FALSE) {
  plotTitle <- "Observed values only\n"
  
  p <- ggplot2::ggplot(data, ggplot2::aes(x = .data[["calendarYear"]], y = .data[[yAxisCol]])) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::scale_y_continuous(limits = c(0, NA)) + # Ensure y-axis starts at 0
    ggplot2::labs(title = plotTitle, x = xLabel, y = "Incidence Rate") +
    ggplot2::theme_minimal()
  
  if (convertToPlotLy) {
    # Convert ggplot object to plotly for interactive visualization
    p_plotly <- plotly::ggplotly(p, tooltip = c("x", "y"))
    return(p_plotly)
  } else {
    return(p)
  }
}




#' Plot Temporal Trend with Expected and Observed Incidence Rates
#'
#' This function generates a line plot comparing expected and observed incidence rates over time.
#'
#' @param data A data frame containing the data to be plotted.
#' @param convertToPlotLy A logical value indicating whether to convert the ggplot object to a plotly object for interactive visualization. Default is FALSE.
#' @param plotSource Options 'cyclops', 'glm' do you want to plot from cyclops or glm?
#' @return A ggplot or plotly object, depending on the value of `convertToPlotLy`.
#' @export
plotTemporalTrendExpectedObserved <- function(data,
                                              convertToPlotLy = FALSE,
                                              plotSource = 'default') {
  if (nrow(data) == 0) {
    return(NULL)
  }
  
  if (plotSource == 'cyclops') {
    data$expected <- data$cyclopsExpected
    data$expectedIncidenceRate <- data$cyclopsExpectedIncidenceRate
    data$p <- data$cyclopsPValue |> max(na.rm = TRUE)
    data$ratio <- data$cyclopsRatio |> max(na.rm = TRUE)
    data$stable <- data$cyclopsStable |> max(na.rm = TRUE)
  } else if (plotSource == 'glm') {
    data$expected <- data$glmExpected
    data$expectedIncidenceRate <- data$glmExpectedIncidenceRate
    data$p <- data$glmPValue |> max(na.rm = TRUE)
    data$ratio <- data$glmRatio |> max(na.rm = TRUE)
    data$stable <- data$glmStable |> max(na.rm = TRUE)
    data$expectedIncidenceRateUpperBound <- data$glmExpectedIncidenceRateUpperBound
    data$expectedIncidenceRateLowerBound <- data$glmExpectedIncidenceRateLowerBound
  } else {
    data$expected <- dplyr::coalesce(data$glmExpected, data$cyclopsExpected)
    data$expectedIncidenceRate <- dplyr::coalesce(data$glmExpectedIncidenceRate,
                                                  data$cyclopsExpectedIncidenceRate)
    data$p <- dplyr::coalesce(data$glmPValue, data$cyclopsPValue) |> max(na.rm = TRUE)
    data$ratio <- dplyr::coalesce(data$glmRatio, data$cyclopsRatio) |> max(na.rm = TRUE)
    data$stable <- dplyr::coalesce(data$glmStable, data$cyclopsStable) |> max(na.rm = TRUE)
    data$expectedIncidenceRateUpperBound <- data$glmExpectedIncidenceRateUpperBound
    data$expectedIncidenceRateLowerBound <- data$glmExpectedIncidenceRateLowerBound
  }
  
  numberOfSplinesUsed <- if ('numberOfSplinesUsed' %in% colnames(data)) {
    max(data$numberOfSplinesUsed |> max(na.rm = TRUE))
  }
  else {
    0
  }
  
  plotTitle <- paste0(
    min(data$cohortName),
    " (id: ",
    min(data$cohortId),
    ")\n",
    min(data$databaseId),
    "\nObserved vs Expected. Splines used = ",
    numberOfSplinesUsed,
    ". p-value: ",
    OhdsiHelpers::formatDecimalWithComma(min(data$p), decimalPlaces = 4),
    ". ratio = ",
    OhdsiHelpers::formatDecimalWithComma(min(data$ratio), decimalPlaces = 4),
    "\n",
    "\n - ",
    min(data$zeroRecordLeadRemoved),
    " records with lead (left) zero counts found and has been removed",
    "\n - ",
    min(data$zeroRecordTrailRemoved),
    " records with trail (right) zero counts found and has been removed"
  )
  
  if (data |> dplyr::filter(useYearData == 0) |> nrow() > 0) {
    reportData <- data |>
      dplyr::filter(useYearData == 1)
    
    expectedMedian <- median(reportData$personYears[-1])
    
    if (any(!is.na(reportData$firstYearToCensor))) {
      plotTitle <- paste0(
        plotTitle ,
        "\n - ",
        " removed data from first calendar year ",
        format(max(
          reportData$firstYearToCensor, na.rm = TRUE
        ), "%Y"),
        " with ",
        max(reportData$firstYearCensoredPersonYears) |> OhdsiHelpers::formatIntegerWithComma(),
        " person years of data. Expected around: ",
        expectedMedian |> OhdsiHelpers::formatIntegerWithComma()
      )
    }
    
    if (any(!is.na(reportData$lastYearToCensor))) {
      plotTitle <- paste0(
        plotTitle ,
        "\n - ",
        " removed data from last calendar year ",
        format(max(
          reportData$lastYearToCensor, na.rm = TRUE
        ), "%Y"),
        " with ",
        max(reportData$lastYearCensoredPersonYears) |> OhdsiHelpers::formatIntegerWithComma(),
        " person years of data. Expected around: ",
        expectedMedian |> OhdsiHelpers::formatIntegerWithComma()
      )
    }
  } else {
    plotTitle <- paste0(plotTitle, "\n - using all calendar years")
  }
  
  if (any(data$stable == FALSE, na.rm = TRUE)) {
    stableTitle <- "\n\nUNSTABLE (FAILED LogLikelihood)"
  } else {
    stableTitle <- "\n\nSTABLE (Did not fail LogLikelihood)"
  }
  
  if (any(!is.na(data$expectedIncidenceRateLowerBound))) {
    stableTitlePearson <- NULL
    if (any(data$glmPearsonStable == FALSE, na.rm = TRUE)) {
      stableTitlePearson <- " [Pearson FAILED]"
    } else {
      stableTitlePearson <- " [Pearson did not fail]"
    }
    
    stableTitleDeviance <- NULL
    if (any(data$glmDevianceStable == FALSE, na.rm = TRUE)) {
      stableTitleDeviance <- " [G-test deviance FAILED]"
    } else {
      stableTitlePearson <- " [G-test deviance did not fail]"
    }
    
    stableTitle <- paste0(stableTitle, stableTitlePearson, stableTitleDeviance)
  }
  
  plotTitle <- if (any(data$stable == FALSE, na.rm = TRUE)) {
    paste0(plotTitle, "\n\n", stableTitle)
  } else {
    paste0(plotTitle, "\n\n", stableTitle)
  }
  
  # Dynamically generate the column names for the confidence interval
  lowerBoundCol <- "expectedIncidenceRateLowerBound"
  upperBoundCol <- "expectedIncidenceRateUpperBound"
  
  # Check if the dynamically named UpperBound and LowerBound columns exist in the data
  if (lowerBoundCol %in% colnames(data) &&
      upperBoundCol %in% colnames(data)) {
    p <- ggplot2::ggplot(data) +
      ggplot2::geom_ribbon(
        ggplot2::aes(
          x = .data[["calendarYear"]],
          ymin = .data[[lowerBoundCol]],
          ymax = .data[[upperBoundCol]]
        ),
        fill = "lightgray",
        alpha = 0.5
      ) + # Confidence interval shading
      ggplot2::geom_line(ggplot2::aes(
        x = .data[["calendarYear"]],
        y = .data[["expectedIncidenceRate"]],
        linetype = "Expected"
      )) +
      ggplot2::geom_line(ggplot2::aes(
        x = .data[["calendarYear"]],
        y = .data[["incidenceRate"]],
        linetype = "Observed"
      )) +
      ggplot2::scale_linetype_manual(values = c("Expected" = "dashed", "Observed" = "solid")) +
      ggplot2::scale_y_continuous(limits = c(0, NA)) + # Ensure y-axis starts at 0
      ggplot2::labs(
        title = plotTitle,
        x = "Calendar Year",
        y = "Incidence Rate",
        linetype = "Incidence Type"
      )
  } else {
    p <- ggplot2::ggplot(data) +
      ggplot2::geom_line(ggplot2::aes(
        x = .data[["calendarYear"]],
        y = .data[["expectedIncidenceRate"]],
        linetype = "Expected"
      )) +
      ggplot2::geom_line(ggplot2::aes(
        x = .data[["calendarYear"]],
        y = .data[["incidenceRate"]],
        linetype = "Observed"
      )) +
      ggplot2::scale_linetype_manual(values = c("Expected" = "dashed", "Observed" = "solid")) +
      ggplot2::scale_y_continuous(limits = c(0, NA)) + # Ensure y-axis starts at 0
      ggplot2::labs(
        title = plotTitle,
        x = "Calendar Year",
        y = "Incidence Rate",
        linetype = "Incidence Type"
      )
  }
  
  p <- p + ggplot2::theme_minimal()
  
  if (convertToPlotLy) {
    # Convert ggplot object to plotly for interactive visualization
    p_plotly <- plotly::ggplotly(p, tooltip = c("x", "y", "linetype"))
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

createPlotsByDatabaseId <- function(data, cohortId) {
  # Extract unique databaseIds
  databaseIds <- data$databaseId |> unique() |> sort()
  
  # Create a list of plots for each databaseId
  plots <- lapply(databaseIds, function(dbId) {
    # Filter data for the given cohortId and databaseId
    subsetData <- data |>
      dplyr::filter(cohortId == !!cohortId & databaseId == !!dbId)
    
    # Create the plot using your function plotTemporalTrendExpectedObserved
    plotTemporalTrendExpectedObserved(subsetData)
  })
  
  return(plots)
}


summarizeTemporalStability <- function(temporalStabilityOutput,
                                       cylopsStatisticallyUnstable) {
  # Calculate the counts for tested cohorts, databases, and combinations
  cohortsTested <- temporalStabilityOutput$cohortId |> unique() |> length()
  databasesTested <- temporalStabilityOutput$databaseId |> unique() |> length()
  combosTested <- temporalStabilityOutput |> dplyr::select(cohortId, databaseId) |> dplyr::distinct() |> nrow()
  
  # Calculate the counts for failed cohorts, databases, and combinations
  cohortsFailed <- length(cylopsStatisticallyUnstable$cohortId |> unique())
  databasesFailed <- length(cylopsStatisticallyUnstable$databaseId |> unique())
  combosFailed <- cylopsStatisticallyUnstable |> dplyr::select(cohortId, databaseId) |> dplyr::distinct() |> nrow()
  
  # Print the summary
  cat(
    "\n",
    "-----",
    "\n",
    "Number of cohorts tested:",
    cohortsTested,
    "\n",
    "Number of databases tested:",
    databasesTested,
    "\n",
    "Number of combos tested:",
    combosTested,
    "\n",
    "-----",
    "\n",
    "Number of cohorts that failed in atleast one datasource:",
    OhdsiHelpers::formatCountPercent(count = cohortsFailed, percent = cohortsFailed /
                                       cohortsTested),
    "\n",
    "Number of databases with atleast one cohort failed:",
    OhdsiHelpers::formatCountPercent(count = databasesFailed, percent = databasesFailed /
                                       databasesTested),
    "\n",
    "Number of combos failed:",
    OhdsiHelpers::formatCountPercent(count = combosFailed, percent = combosFailed /
                                       combosTested),
    "\n",
    "-----",
    "\n"
  )
}
