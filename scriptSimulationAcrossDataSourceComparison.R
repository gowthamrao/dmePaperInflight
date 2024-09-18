# Install the mgcv package if not already installed
# install.packages("mgcv")

library(mgcv)
library(ggplot2)

# Enhanced Sample data with 10 years of data per database using mid-year dates
set.seed(123)  # For reproducibility

years <- seq(as.Date("2010-06-01"), as.Date("2019-06-01"), by = "year")  # Mid-year dates
databases <- paste0('DB', 1:10)

# Generate random expected and observed counts for each database and year
data <- expand.grid(
  databaseId = databases,
  calendarYear = years
)

# Create synthetic counts based on Poisson distribution
data$expectedCount <- round(runif(nrow(data), 100, 200))  # Random expected counts between 100 and 200
data$observedCount <- rpois(nrow(data), lambda = data$expectedCount)  # Observed counts from Poisson

# Calculate the rate ratio
data$rateRatio <- data$observedCount / data$expectedCount

# Convert the 'calendarYear' date to numeric for the GAM model (necessary for the spline function)
data$calendarYearNumeric <- as.numeric(data$calendarYear)

# Fit a GAM with factor-smooth interaction for nonlinear trends across calendarYear
gam_model <- gam(rateRatio ~ s(calendarYearNumeric, bs = "fs", k = 5, m = 1), data = data)

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

# Visualize: Spaghetti plot showing the ratio and fitted values over time for each database
ggplot(data, aes(x = calendarYear, y = rateRatio, group = databaseId, color = databaseId)) +
  geom_line() +
  geom_line(aes(y = fitted_values), linetype = "dashed") +
  
  # Add 99% confidence interval ribbon with high transparency
  geom_ribbon(aes(ymin = lower_bound_99, ymax = upper_bound_99, fill = databaseId), alpha = 0.01) +
  
  # Add 95% confidence interval ribbon with moderate transparency
  geom_ribbon(aes(ymin = lower_bound_95, ymax = upper_bound_95, fill = databaseId), alpha = 0.05) +
  
  labs(title = "Observed vs. Fitted Rate Ratios with 95% and 99% Confidence Intervals",
       y = "Rate Ratio (Observed / Expected)",
       x = "Calendar Year",
       color = "Database ID", fill = "Database ID") +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year") +  # Properly format the date axis
  theme_minimal()

# Identify outliers based on residuals (e.g., using z-scores for residuals)
data$residual_zscore <- (data$residuals - mean(data$residuals)) / sd(data$residuals)

# Threshold for outliers (|z| > 2 is a common cutoff)
outliers <- data[abs(data$residual_zscore) > 2, ]

# Highlight outliers
ggplot(data, aes(x = calendarYear, y = rateRatio, group = databaseId, color = databaseId)) +
  geom_line() +
  geom_line(aes(y = fitted_values), linetype = "dashed") +
  
  # Add 99% confidence interval ribbon with high transparency
  geom_ribbon(aes(ymin = lower_bound_99, ymax = upper_bound_99, fill = databaseId), alpha = 0.01) +
  
  # Add 95% confidence interval ribbon with moderate transparency
  geom_ribbon(aes(ymin = lower_bound_95, ymax = upper_bound_95, fill = databaseId), alpha = 0.05) +
  
  # Highlight outliers in red
  geom_point(data = outliers, aes(x = calendarYear, y = rateRatio), color = "red", size = 3) +
  
  labs(title = "Rate Ratios with Outliers Highlighted and 95% and 99% Confidence Intervals",
       y = "Rate Ratio (Observed / Expected)",
       x = "Calendar Year",
       color = "Database ID", fill = "Database ID") +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year") +  # Properly format the date axis
  theme_minimal()
