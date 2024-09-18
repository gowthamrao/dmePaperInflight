performGAMAnalysis <- function(data) {
  library(mgcv)
  library(ggplot2)
  
  # Ensure that calendarYear is in date format if not already
  if (!inherits(data$calendarYear, "Date")) {
    stop("calendarYear must be in Date format.")
  }
  
  # Check that required columns exist
  required_cols <- c("databaseId", "calendarYear", "expectedObservedRatio")
  if (!all(required_cols %in% names(data))) {
    stop("The data must contain columns: databaseId, calendarYear, expectedObservedRatio.")
  }
  
  # Convert the 'calendarYear' date to numeric for the GAM model
  data$calendarYearNumeric <- as.numeric(data$calendarYear)
  
  # Fit a GAM with factor-smooth interaction for nonlinear trends across calendarYear
  gam_model <- gam(expectedObservedRatio ~ s(
    calendarYearNumeric,
    bs = "fs",
    k = 5,
    m = 1
  ),
  data = data)
  
  # Predict the fitted values and standard errors
  predictions <- predict(gam_model, newdata = data, se.fit = TRUE)
  
  # Add fitted values, standard errors, and confidence intervals to the data
  data$fitted_values <- predictions$fit
  data$upper_bound_95 <- data$fitted_values + 1.96 * predictions$se.fit  # 95% upper bound
  data$lower_bound_95 <- data$fitted_values - 1.96 * predictions$se.fit  # 95% lower bound
  data$upper_bound_99 <- data$fitted_values + 2.576 * predictions$se.fit  # 99% upper bound
  data$lower_bound_99 <- data$fitted_values - 2.576 * predictions$se.fit  # 99% lower bound
  
  # Calculate residuals to assess deviations from the model
  data$residuals <- residuals(gam_model)
  
  # Calculate residual z-scores
  data$residual_zscore <- (data$residuals - mean(data$residuals)) / sd(data$residuals)
  
  # Threshold for outliers (|z| > 2 is a common cutoff)
  outliers <- data[abs(data$residual_zscore) > 2, ]
  
  # Identify unique databaseIds that are outliers
  outlier_databaseIds <- unique(outliers$databaseId)
  outlier_message <- if (length(outlier_databaseIds) > 0) {
    paste(paste(outlier_databaseIds, collapse = ",\n"))
  } else {
    "No outliers detected"
  }
  
  # Plot 1: Spaghetti plot with confidence intervals
  plot1 <- ggplot(
    data,
    aes(
      x = calendarYear,
      y = expectedObservedRatio,
      group = databaseId,
      color = databaseId
    )
  ) +
    geom_line() +
    geom_line(aes(y = fitted_values), linetype = "dashed") +
    geom_ribbon(aes(
      ymin = lower_bound_99,
      ymax = upper_bound_99,
      fill = databaseId
    ),
    alpha = 0.05) +
    geom_ribbon(aes(
      ymin = lower_bound_95,
      ymax = upper_bound_95,
      fill = databaseId
    ),
    alpha = 0.15) +
    labs(
      title = paste("Failed: ", outlier_message),
      y = "Expected/Observed Ratio",
      x = "Calendar Year",
      color = "Database ID",
      fill = "Database ID"
    ) +
    scale_x_date(date_labels = "%Y", date_breaks = "1 year") +  # Properly format the date axis
    theme_minimal()
  
  # Plot 2: Spaghetti plot with outliers highlighted
  plot2 <- ggplot(
    data,
    aes(
      x = calendarYear,
      y = expectedObservedRatio,
      group = databaseId,
      color = databaseId
    )
  ) +
    geom_line() +
    geom_line(aes(y = fitted_values), linetype = "dashed") +
    geom_ribbon(aes(
      ymin = lower_bound_99,
      ymax = upper_bound_99,
      fill = databaseId
    ),
    alpha = 0.05) +
    geom_ribbon(aes(
      ymin = lower_bound_95,
      ymax = upper_bound_95,
      fill = databaseId
    ),
    alpha = 0.15) +
    geom_point(
      data = outliers,
      aes(x = calendarYear, y = expectedObservedRatio),
      color = "red",
      size = 3
    ) +
    labs(
      title = paste("Failed: ", outlier_message),
      y = "Expected/Observed Ratio",
      x = "Calendar Year",
      color = "Database ID",
      fill = "Database ID"
    ) +
    scale_x_date(date_labels = "%Y", date_breaks = "1 year") +  # Properly format the date axis
    theme_minimal()
  
  # Return the model and the plots in a list
  return(list(
    model = gam_model,
    plot1 = plot1,
    plot2 = plot2
  ))
}
