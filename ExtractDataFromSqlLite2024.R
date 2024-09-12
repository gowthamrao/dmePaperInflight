

connectionDetails <- DatabaseConnector::createConnectionDetails(dbms = "sqlite", server = "E:\\studyResults\\ohdsiPlDme\\CohortDiagnostics\\MergedCohortDiagnosticsData.sqlite")

connection <-
  DatabaseConnector::connect(connectionDetails = connectionDetails)

#extract data----
incidenceRate2024 <- DatabaseConnector::renderTranslateQuerySql(connection = connection,
                                                                sql = "SELECT * FROM main.incidence_rate;",
                                                                snakeCaseToCamelCase = TRUE) |> dplyr::tibble()
saveRDS(object = incidenceRate2024, file = "incidenceRate2024.RDS")

cohort2024 <- DatabaseConnector::renderTranslateQuerySql(connection = connection,
                                                         sql = "SELECT * FROM cohort;",
                                                         snakeCaseToCamelCase = TRUE) |> dplyr::tibble()
saveRDS(object = cohort2024, file = "cohort2024.RDS")


cohortCount2024 <- DatabaseConnector::renderTranslateQuerySql(connection = connection,
                                                              sql = "SELECT * FROM main.cohort_count;",
                                                              snakeCaseToCamelCase = TRUE) |> dplyr::tibble()
saveRDS(object = cohortCount2024, file = "cohortCount2024.RDS")

DatabaseConnector::disconnect(connection)
