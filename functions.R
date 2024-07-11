processData <- function(data,
                        xAxisCol = "calendarYear",
                        baseLineCol = "baselineIncidenceRate",
                        yAxisCol = "incidenceRate",
                        groupBy = "databaseId",
                        rescaleMethod = "none",
                        numberOfSplines = NULL) {
  data <- data |>
    dplyr::select(
      !!rlang::sym(groupBy),!!rlang::sym(xAxisCol),!!rlang::sym(yAxisCol),!!rlang::sym(baseLineCol)
    )
  
  # Updated the list of acceptable methods to include 'zScore'
  acceptableMethods <-
    c("none",
      "minMax",
      "absoluteChange",
      "percentChange",
      "log",
      "zScore")
  
  # Validate rescaleMethod is acceptable
  if (!rescaleMethod %in% acceptableMethods) {
    stop(paste(
      "rescaleMethod must be one of the following:",
      toString(acceptableMethods)
    ))
  }
  
  # Only apply transformations if rescaleMethod is not 'none'
  if (rescaleMethod != "none") {
    data <- data |>
      dplyr::group_by(!!rlang::sym(groupBy)) |>
      dplyr::mutate(!!yAxisCol := {
        dplyr::case_when(
          rescaleMethod == "minMax" ~ {
            -1 + 2 * ((.data[[yAxisCol]] - min(.data[[yAxisCol]], na.rm = TRUE)) /
                        (max(.data[[yAxisCol]], na.rm = TRUE) - min(.data[[yAxisCol]], na.rm = TRUE)))
          },
          rescaleMethod == "absoluteChange" &&
            !is.null(baseLineCol) ~ {
              .data[[yAxisCol]] - .data[[baseLineCol]]
            },
          rescaleMethod == "percentChange" &&
            !is.null(baseLineCol) ~ {
              100 * (.data[[yAxisCol]] - .data[[baseLineCol]]) / .data[[baseLineCol]]
            },
          rescaleMethod == "log" ~ {
            minValue <-
              min(.data[[yAxisCol]][.data[[yAxisCol]] > 0], na.rm = TRUE)
            shiftValue <-
              ifelse(minValue <= 0, 0.000000001, minValue)
            log10(ifelse(.data[[yAxisCol]] <= 0, shiftValue, .data[[yAxisCol]]))
          },
          rescaleMethod == "zScore" ~ {
            # Z-score normalization
            (.data[[yAxisCol]] - mean(.data[[yAxisCol]], na.rm = TRUE)) / sd(.data[[yAxisCol]], na.rm = TRUE)
          },
          TRUE ~ .data[[yAxisCol]]  # Retain original data if no other condition matches
        )
      })
  }
  
  # creates a new column called timeId based on xAxisCol
  data <- data |>
    dplyr::inner_join(
      data |>
        dplyr::select(!!rlang::sym(xAxisCol)) |>
        dplyr::distinct() |>
        dplyr::arrange(!!rlang::sym(xAxisCol)) |>
        dplyr::mutate(timeId = dplyr::row_number())
    )
  
  # modify such that for each group the output of stats::predict creates a new column called "predicted_"yAxisCol
  data_grouped <- split(data, data[[groupBy]])
  predictions <- lapply(data_grouped, function(df) {
    # Dynamic assignment of numberOfSplines
    numberOfSplinesForGroup <- numberOfSplines
    if (is.null(numberOfSplinesForGroup)) {
      numberOfSplinesForGroup <- 5
    }
    numberOfSplinesForGroup <-
      max(numberOfSplinesForGroup, ceiling(length(unique(df$timeId)) / 2))
    
    model <- Cyclops::createCyclopsData(formula = as.formula(
      paste(
        yAxisCol,
        "~ splines::ns(x = timeId, df =",
        numberOfSplinesForGroup,
        ")"
      )
    ),
    data = df,
    modelType = "pr") |>
      Cyclops::fitCyclopsModel()
    # can prediction have confidence intervals?
    stats::predict(model, newdata = df)
  })
  
  # Merge predictions back into the original data frame
  predictedColName <- paste("predicted_", yAxisCol, sep = "")
  for (i in names(predictions)) {
    data_grouped[[i]][[predictedColName]] <- predictions[[i]]
  }
  
  data <- do.call(rbind, data_grouped) |>
    dplyr::tibble() |>
    dplyr::mutate(deviationFromExpected = !!rlang::sym(predictedColName) -!!rlang::sym(yAxisCol)) |>
    dplyr::ungroup() |>
    dplyr::select(-timeId)
  
  return(data)
}


plotData <- function(data,
                     xAxisCol = "calendarYear",
                     yAxisCol = "incidenceRate",
                     groupBy = "databaseId",
                     plotTitle = "Incidence Rate over time",
                     xLabel = "Calendar year",
                     yLabel = "Incidence Rate",
                     useMinimalTheme = TRUE) {
  library(ggplot2)
  library(plotly)
  # Generate color palette based on unique database IDs
  numColors <- length(unique(data[[groupBy]]))
  colors <- OhdsiRPlots::createOhdsiPalette(numColors = numColors)
  colorMapping <- setNames(colors, unique(data[[groupBy]]))
  
  p <-
    ggplot(data, aes(x = .data[[xAxisCol]], y = .data[[yAxisCol]], color = .data[[groupBy]])) +
    geom_line() +
    geom_point() +
    scale_color_manual(values = colorMapping) +
    labs(
      title = plotTitle,
      x = xLabel,
      y = yLabel,
      color = "Database ID"
    )
  
  if (useMinimalTheme) {
    p <- p + theme_minimal()
  }
  
  # Convert ggplot object to plotly for interactive visualization
  p_plotly <-
    ggplotly(p, tooltip = c("x", "y", groupBy)) # Set tooltip to show data from specific columns
  return(p_plotly)
}




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