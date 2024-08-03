#set up----
setwd(rstudioapi::getActiveDocumentContext()$path |> dirname())
source("functions.R")


# paste0("'", incidenceRate$databaseId |> unique(), "'") |> paste0(collapse = ", ")
databaseIdsOfInterest <- c(
  'cdm_ims_australia_lpd_v2353',
  'cdm_ims_france_v2354',
  'cdm_ims_germany_v2352',
  'cdm_iqvia_pharmetrics_plus_v2525',
  'cdm_jmdc_v2554',
  'cdm_optum_ehr_v2577',
  'cdm_optum_extended_dod_v2553',
  'cdm_premier_v2543',
  'cdm_truven_ccae_v2542',
  'cdm_truven_mdcd_v2565',
  'cdm_truven_mdcr_v2540'
)

yearRange <- c(2001:2024) |> as.character()


# what is the name of the schema you want to upload to?
resultsSchema <- 'Icpe2023Dme' # change to your schema

# Postgres server: connection details to OHDSI Phenotype library. Please change to your postgres connection details
# connectionDetails <- DatabaseConnector::createConnectionDetails(
#   dbms = Sys.getenv("shinydbDbms", unset = "postgresql"),
#   server = paste(
#     Sys.getenv("phenotypeLibraryServer"),
#     Sys.getenv("phenotypeLibrarydb"),
#     sep = "/"
#   ),
#   port = Sys.getenv("phenotypeLibraryDbPort"),
#   user = Sys.getenv("phenotypeLibrarydbUser"),
#   password = Sys.getenv("phenotypeLibrarydbPw")
# )
#
# connection <-
#   DatabaseConnector::connect(connectionDetails = connectionDetails)


#extract data----
# incidenceRate <- DatabaseConnector::renderTranslateQuerySql(
#   connection = connection,
#   sql = "SELECT * FROM @results_database_schema.incidence_rate;",
#   results_database_schema = resultsSchema,
#   snakeCaseToCamelCase = TRUE
# ) |> dplyr::tibble()
# saveRDS(object = incidenceRate, file = "incidenceRate.rds")
incidenceRate <- readRDS(file = "incidenceRate.rds")|>
  dplyr::filter(databaseId %in% databaseIdsOfInterest,
                calendarYear %in% yearRange)|>
  dplyr::filter(calendarYear != "") |>
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
#   sql = "SELECT * FROM @results_database_schema.cohort_count;",
#   results_database_schema = resultsSchema,
#   snakeCaseToCamelCase = TRUE
# ) |> dplyr::tibble()
# saveRDS(object = cohortCount, file = "cohortCount.rds")
cohortCount <- readRDS(file = "cohortCount.rds")


#R1: cohort counts----
# cohortCount |>
#   tidyr::pivot_wider(
#     id_cols = c(cohortId),
#     values_from = cohortSubjects,
#     names_from = databaseId
#   ) |>
#   dplyr::inner_join(cohort |>
#                       dplyr::select(cohortId, cohortName)) |>
#   dplyr::relocate(cohortId, cohortName) |>
#   dplyr::arrange(dplyr::desc(cdm_optum_extended_dod_v2553)) |>
#   View()



#R2: stability----
temporalStabilityOutputWithZeroCount <- checkTemporalStabilityForcohortDiagnosticsIncidenceRateData(
  cohortDiagnosticsIncidenceRateData = incidenceRate,
  cohort = cohort,
  maxNumberOfSplines = NULL,
  splineTickInterval = 3
)

if (all(na.omit(temporalStabilityOutputWithoutZeroCount$stable))) {
  writeLines("All incidence rate data are stable over time")
} else {
  writeLines("Some incidence rate data are not stable over time")
  temporalStabilityOutputWithZeroCount |>
    dplyr::filter(stable == FALSE) |>
    View()
}

temporalStabilityOutputWithoutZeroCount <- checkTemporalStabilityForcohortDiagnosticsIncidenceRateData(
  cohortDiagnosticsIncidenceRateData = incidenceRate,
  cohort = cohort,
  filterZeroAndLowCountRecords = TRUE,
  maxNumberOfSplines = NULL,
  splineTickInterval = 3
)
if (all(na.omit(temporalStabilityOutputWithoutZeroCount$stable))) {
  writeLines("All incidence rate data are stable over time")
} else {
  writeLines("Some incidence rate data are not stable over time")
  temporalStabilityOutputWithoutZeroCount |>
    dplyr::filter(stable == FALSE) |>
    View()
}

temporalStabilityOutput <- dplyr::bind_rows(
  temporalStabilityOutputWithZeroCount |> dplyr::mutate(useZeroCount = 1),
  temporalStabilityOutputWithoutZeroCount |> dplyr::mutate(useZeroCount = 0)
)

# R3: do some simple plots
cohortIds <- incidenceRate$cohortId |> unique() |> sort()
plotsIncidenceRateSimple <- c()
for (i in (1:length(cohortIds))) {
  data <- incidenceRate |>
    dplyr::filter(cohortId == cohortIds[i])
  
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
ggplot2::ggsave("simpleTrendPlots.pdf",
                combined_plot,
                width = 11,
                height = 8.5)  # Adjust the width and height as needed

