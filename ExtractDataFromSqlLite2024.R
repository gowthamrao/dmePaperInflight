
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
#   sql = "SELECT * FROM cohort;",
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