#set up----
setwd(rstudioapi::getActiveDocumentContext()$path |> dirname())
source("functions.R")

# connectionDetails <- DatabaseConnector::createConnectionDetails(
#   dbms = "sqlite",
#   server = "E:\\studyResults\\ohdsiPlDme\\CohortDiagnostics\\MergedCohortDiagnosticsData.sqlite"
# )
# 
# connection <-
#   DatabaseConnector::connect(connectionDetails = connectionDetails)
# 
# databases <- DatabaseConnector::renderTranslateQuerySql(connection = connection,
#                                                         sql = "SELECT * FROM database;",
#                                                         snakeCaseToCamelCase = TRUE) |> 
#   dplyr::tibble()

databaseIdsOfInterest <- c(
  'cdm_ims_australia_lpd_v3006',
  'cdm_ims_france_v2914',
  'cdm_ims_germany_v2915',
  'cdm_iqvia_pharmetrics_plus_v3403',
  'cdm_jmdc_v3044',
  'cdm_optum_ehr_v3037',
  'cdm_optum_extended_dod_v3039',
  'cdm_truven_ccae_v3406',
  'cdm_truven_mdcd_v3038',
  'cdm_truven_mdcr_v3045'
)

yearRange <- c(2001:2024) |> as.character()

#extract data----
# incidenceRate <- DatabaseConnector::renderTranslateQuerySql(
#   connection = connection,
#   sql = "SELECT * FROM main.incidence_rate;",
#   snakeCaseToCamelCase = TRUE
# ) |> dplyr::tibble()
# saveRDS(object = incidenceRate, file = "incidenceRate.rds")
incidenceRate <- readRDS(file = "incidenceRate.rds")|>
  dplyr::filter(databaseId %in% databaseIdsOfInterest,
                calendarYear %in% yearRange)|>
  dplyr::filter(calendarYear != "") |>
  tidyr::replace_na(list(ageGroup = '', gender = '')) |> 
  dplyr::filter(ageGroup == '') |> # we are not using age stratified incidence rate for the diagnostic
  dplyr::filter(gender == '') # we are not using gender stratified incidence rate for the diagnostic

# cohort <- DatabaseConnector::renderTranslateQuerySql(
#   connection = connection,
#   sql = "SELECT * FROM @results_database_schema.cohort;",
#   results_database_schema = resultsSchema,
#   snakeCaseToCamelCase = TRUE
# ) |> dplyr::tibble()
# saveRDS(object = cohort, file = "cohort.rds")
cohort <- readRDS(file = "cohort.rds")


# cohortCount <- DatabaseConnector::renderTranslateQuerySql(
#   connection = connection,
#   sql = "SELECT * FROM main.cohort_count;",
#   snakeCaseToCamelCase = TRUE
# ) |> dplyr::tibble()
# saveRDS(object = cohortCount, file = "cohortCount.rds")
cohortCount <- readRDS(file = "cohortCount.rds")

# cohort <- DatabaseConnector::renderTranslateQuerySql(
#   connection = connection,
#   sql = "SELECT * FROM main.cohort;",
#   snakeCaseToCamelCase = TRUE
# ) |> dplyr::tibble()
# saveRDS(object = cohort, file = "cohort.rds")
cohort <- readRDS(file = "cohort.rds")


#R1: cohort counts----
cohortCount |>
  tidyr::pivot_wider(
    id_cols = c(cohortId),
    values_from = cohortSubjects,
    names_from = databaseId
  ) |>
  dplyr::inner_join(cohort |>
                      dplyr::select(cohortId, cohortName)) |>
  dplyr::relocate(cohortId, cohortName) |>
  dplyr::arrange(dplyr::desc(cdm_optum_extended_dod_v3039)) |>
  View()


#R2: stability----
temporalStabilityOutput <- checkTemporalStabilityForcohortDiagnosticsIncidenceRateData(
  cohortDiagnosticsIncidenceRateData = incidenceRate,
  cohort = cohort,
  maxNumberOfSplines = 5,
  splineTickInterval = 3
)

report <- temporalStabilityOutput |>
  dplyr::select(cohortId, cohortName, databaseId, stable) |>
  dplyr::distinct()  |>
  dplyr::mutate(stable = dplyr::case_when(
    is.na(stable) ~ "Not tested",
    stable == FALSE ~ "Unstable",
    TRUE ~ ""
  )) |> 
  tidyr::pivot_wider(
    id_cols = c(cohortId, cohortName),
    names_from = databaseId,
    values_from = stable
  )

if (all(na.omit(temporalStabilityOutput$stable))) {
  writeLines("All incidence rate data are stable over time")
} else {
  writeLines("Some incidence rate data are not stable over time")
  temporalStabilityOutput |>
    dplyr::filter(stable == FALSE) |>
    View()
}


temporalStabilityOutput |>
  dplyr::select(cohortId, cohortName, databaseId, ratio, p, stable) |>
  dplyr::distinct() |>
  # dplyr::filter(stable == FALSE) |>
  dplyr::arrange(stable, cohortId, cohortName, databaseId) |>
  dplyr::mutate(
    ratio = OhdsiHelpers::formatDecimalWithComma(ratio),
    p = OhdsiHelpers::formatDecimalWithComma(p, decimalPlaces = 5)
  ) |> View()


#R3 plot observed and expected

# R3: do some simple plots - all database

dataToTryPlot <- temporalStabilityOutput |> 
  dplyr::mutate(databaseName = databaseId) |> 
  dplyr::filter(cohortId == 207)


undebug(createTrellisOfPlots)
createTrellisOfPlots(data = dataToTryPlot)





cohortIds <- temporalStabilityOutput$cohortId |> unique() |> sort()
databaseIds <- temporalStabilityOutput$databaseId |> unique() |> sort()

plotsIncidenceRateCompareObservedExpected <- c()
for (i in (1:length(cohortIds))) {
  
  data <- temporalStabilityOutput |>
    dplyr::filter(cohortId == cohortIdToReport,
                  databaseId == databaseIdToReport) |>
    dplyr::filter(evaluated == TRUE)
  
  
  
  plotsIncidenceRateCompareObservedExpected[[cohortIdToReport]][[databaseIdToReport]] <- createTemporalStabilityPlot(data, bgColor = color)
  
}



# R3: do some simple plots - all database
cohortIds <- incidenceRate$cohortId |> unique() |> sort()
plotsIncidenceRateSimple <- c()
for (i in (1:length(cohortIds))) {
  cohortIdToReport <- cohortIds[i]
  
  cohortIdToReport <- 222
  
  data <- incidenceRate |>
    dplyr::filter(cohortId == cohortIdToReport)
  
  checkIfCohortDiagnosticsIncidenceRateData(data)
  
  incidenceRateData <- processCohortDiagnosticsIncidenceRateData(cohortDiagnosticsIncidenceRateData = data)
  
  plotTitle <- paste0(
    paste("Cohort ID:", cohortIdToReport),
    "\n",
    paste0(
      cohort |>
        dplyr::filter(cohortId == cohortIdToReport) |>
        dplyr::pull(cohortName)
    )
  )
  
  plotsIncidenceRateSimple[[i]] <- incidenceRateData |>
    plotSimpleTemporalTrend(plotTitle = plotTitle, yAxisCol = "incidenceRate")
}

# Create a list of plots
plot_list <- plotsIncidenceRateSimple

# Use marrangeGrob to arrange the plots across multiple pages
combined_plot <- gridExtra::marrangeGrob(plot_list, ncol = 1, nrow = 1)

# Now you can save the combined plot as PDF
ggplot2::ggsave("simpleTrendPlots.pdf",
                combined_plot,
                width = 11,
                height = 8.5)  # Adjust the width and height as needed


# R4: do some simple plots - remove premier
cohortIds <- incidenceRate$cohortId |> unique() |> sort()
plotsIncidenceRateSimple <- c()
for (i in (1:length(cohortIds))) {
  data <- incidenceRate |>
    dplyr::filter(cohortId == cohortIds[i]) |>
    dplyr::filter(databaseId != 'cdm_premier_v2543')
  
  checkIfCohortDiagnosticsIncidenceRateData(data)
  
  incidenceRateData <- processCohortDiagnosticsIncidenceRateData(cohortDiagnosticsIncidenceRateData = data)
  
  plotTitle <- paste0(
    paste("Cohort ID:", cohortIds[i]),
    "\n",
    paste0(
      cohort |>
        dplyr::filter(cohortId == cohortIds[i]) |>
        dplyr::pull(cohortName)
    )
  )
  
  plotsIncidenceRateSimple[[i]] <- incidenceRateData |>
    plotSimpleTemporalTrend(plotTitle = plotTitle, yAxisCol = "incidenceRate")
}

# Create a list of plots
plot_list <- plotsIncidenceRateSimple

# Use marrangeGrob to arrange the plots across multiple pages
combined_plot <- gridExtra::marrangeGrob(plot_list, ncol = 1, nrow = 1)

# Now you can save the combined plot as PDF
ggplot2::ggsave(
  "simpleTrendPlotsNoPremier.pdf",
  combined_plot,
  width = 11,
  height = 8.5
)  # Adjust the width and height as needed


# todo
# plot expectedObserved to illustrate reasons for deviation.
