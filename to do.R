





# plotData <- function(data,
#                      baseLineCol = NULL,
#                      baseLineAdjustment = NULL,
#                      xAxisCol,
#                      yAxisCol,
#                      groupBy = "databaseId",
#                      logTransform = FALSE,
#                      rescaleAbsolute = NULL,
#                      plotTitle = "Incidence Rate over time",
#                      xLabel = "Calendar year",
#                      yLabel = "Incidence Rate",
#                      colorMapping,  # Added colorMapping as a parameter
#                      useMinimalTheme = TRUE) {
#
#   # Group data by specified column
#   data <- data |>
#     dplyr::group_by(!!rlang::sym(groupBy))
#
#   # Check if colorMapping has enough colors for the unique IDs
#   uniqueIds <- unique(data[[groupBy]])
#   if (length(colorMapping) < length(uniqueIds)) {
#     stop("Insufficient colors provided for the number of unique IDs.")
#   }
#
#   # Check rescaleAbsolute option values
#   if (!is.null(rescaleAbsolute)) {
#     checkmate::assert_choice(tolower(rescaleAbsolute), choices = c("minmax", "zscore"))
#   }
#
#   # Apply rescaling transformations
#   if (!is.null(rescaleAbsolute)) {
#     data <- data |>
#       dplyr::mutate(!!rlang::sym(yAxisCol) := {
#         if (tolower(rescaleAbsolute) == "minmax") {
#           # Min-Max scaling to range -1 to 1
#           -1 + 2 * ((.data[[yAxisCol]] - min(.data[[yAxisCol]], na.rm = TRUE)) /
#                       (max(.data[[yAxisCol]], na.rm = TRUE) - min(.data[[yAxisCol]], na.rm = TRUE)))
#         } else if (tolower(rescaleAbsolute) == "zscore") {
#           # Z-score normalization
#           (.data[[yAxisCol]] - mean(.data[[yAxisCol]], na.rm = TRUE)) / sd(.data[[yAxisCol]], na.rm = TRUE)
#         } else {
#           .data[[yAxisCol]]
#         }
#       })
#   }
#
#   # Baseline adjustments
#   if (!is.null(baseLineCol) && !is.null(baseLineAdjustment)) {
#     checkmate::assert_choice(tolower(baseLineAdjustment), choices = c("ancova", "change", "percent", "mixed"))
#     data <- data |>
#       dplyr::mutate(!!rlang::sym(yAxisCol) := {
#         baseline <- .data[[baseLineCol]]
#         current <- .data[[yAxisCol]]
#         if (tolower(baseLineAdjustment) == "change") {
#           current - baseline
#         } else if (tolower(baseLineAdjustment) == "percent") {
#           100 * (current - baseline) / baseline
#         } else {
#           current  # Placeholder for 'ancova' and 'mixed', which are more complex and typically handled differently
#         }
#       })
#   }
#
#   # Adjust data for log transformation if necessary
#   if (logTransform) {
#     minValue <- min(data[[yAxisCol]][data[[yAxisCol]] > 0], na.rm = TRUE)
#     shiftValue <- ifelse(minValue <= 0, 0.0001, minValue)
#     data <- data |>
#       dplyr::mutate(!!rlang::sym(yAxisCol) := log10(ifelse(.data[[yAxisCol]] <= 0, shiftValue, .data[[yAxisCol]])))
#   }
#
#   # Start plotting
#   p <- ggplot2::ggplot(data, ggplot2::aes_string(x = xAxisCol, y = yAxisCol, color = groupBy)) +
#     ggplot2::geom_line() +
#     ggplot2::geom_point() +
#     ggplot2::scale_color_manual(values = colorMapping)
#
#   # Conditionally apply log transformation on the y-axis
#   if (logTransform) {
#     p <- p + ggplot2::scale_y_log10(
#       breaks = scales::trans_breaks("log10", function(x) 10 ^ x),
#       labels = scales::trans_format("log10", scales::math_format(10^.x))
#     )
#   }
#
#   # Add titles and labels
#   p <- p + ggplot2::labs(
#     title = plotTitle,
#     x = xLabel,
#     y = yLabel,
#     color = "Database ID"
#   )
#
#   # Optionally use a minimal theme
#   if (useMinimalTheme) {
#     p <- p + ggplot2::theme_minimal()
#   }
#
#   return(p)
# }



#
#
#
#
#
#
#
#
#
#
# #P1: no rescale----
# data |>
#   processData() |>
#   plotData(plotTitle = "1 IR raw")
#
# #P2: abs change plot----
# a <- data |>
#   processData(rescale = "absoluteChange") |>
#   plotData(plotTitle = "2 IR rescale absolute change to baseline")
#
# #P3: percent change plot----
# data |>
#   processData(rescale = "percentChange") |>
#   plotData(plotTitle = "3 IR rescale percent change to baseline")
#
# #P4: percent change plot raw data, no adjustment----
# data |>
#   processData(rescale = "minMax") |>
#   plotData(plotTitle = "4 IR rescale min max")
#
# #P5: log change plot----
# data |>
#   processData(rescale = "log") |>
#   plotData(plotTitle = "5 IR rescale min max")
#
# #P6: zscore plot----
# data |>
#   processData(rescale = "zScore") |>
#   plotData(plotTitle = "6 IR rescale zscore")
#
# #P7: zscore on absoluteChange----
# data |>
#   processData(rescale = "minMax") |>
#   processData(rescale = "zScore") |>
#   plotData(plotTitle = "7 IR rescale zscore")
#
# #P8: minMax and percentChange----
# data |>
#   processData(rescale = "minMax") |>
#   processData(rescale = "percentChange") |>
#   plotData(plotTitle = "8 IR min max percent change")
#
#
#
#
# #P11: no rescale----
# data |>
#   processData() |>
#   plotData(plotTitle = "11 IR raw - deviation from expected", yAxisCol = "deviationFromExpected")
#
# #P12: abs change plot----
# data |>
#   processData(rescale = "absoluteChange") |>
#   plotData(plotTitle = "12 IR rescale absolute change to baseline - deviation from expected",
#            yAxisCol = "deviationFromExpected")
#
# #P13: percent change plot----
# data |>
#   processData(rescale = "percentChange") |>
#   plotData(plotTitle = "13 IR rescale percent change to baseline - deviation from expected",
#            yAxisCol = "deviationFromExpected")
#
# #P14: percent change plot raw data, no adjustment----
# data |>
#   processData(rescale = "minMax") |>
#   plotData(plotTitle = "14 IR rescale min max - deviation from expected",
#            yAxisCol = "deviationFromExpected")
#
# #P15: log change plot----
# data |>
#   processData(rescale = "log") |>
#   plotData(plotTitle = "15 IR rescale min max - deviation from expected",
#            yAxisCol = "deviationFromExpected")
#
# #P16: zscore plot----
# data |>
#   processData(rescale = "zScore") |>
#   plotData(plotTitle = "16 IR rescale zscore - deviation from expected",
#            yAxisCol = "deviationFromExpected")
#
# #P17: zscore on absoluteChange----
# data |>
#   processData(rescale = "minMax") |>
#   processData(rescale = "zScore") |>
#   plotData(plotTitle = "17 IR rescale zscore - deviation from expected",
#            yAxisCol = "deviationFromExpected")
#
# #P18: minMax and percentChange----
# data |>
#   processData(rescale = "minMax") |>
#   processData(rescale = "percentChange") |>
#   plotData(plotTitle = "18 IR min max percent change - deviation from expected",
#            yAxisCol = "deviationFromExpected")


