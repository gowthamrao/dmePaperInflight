#set up----
setwd(rstudioapi::getActiveDocumentContext()$path |> dirname())
source("functions.R")

# what is the name of the schema you want to upload to?
resultsSchema <- 'Icpe2023Dme' # change to your schema

# Postgres server: connection details to OHDSI Phenotype library. Please change to your postgres connection details
connectionDetails <- DatabaseConnector::createConnectionDetails(
  dbms = Sys.getenv("shinydbDbms", unset = "postgresql"),
  server = paste(
    Sys.getenv("phenotypeLibraryServer"),
    Sys.getenv("phenotypeLibrarydb"),
    sep = "/"
  ),
  port = Sys.getenv("phenotypeLibraryDbPort"),
  user = Sys.getenv("phenotypeLibrarydbUser"),
  password = Sys.getenv("phenotypeLibrarydbPw")
)

connection <-
  DatabaseConnector::connect(connectionDetails = connectionDetails)

#extract data----
incidenceRateIcpe2023 <- DatabaseConnector::renderTranslateQuerySql(
  connection = connection,
  sql = "SELECT * FROM @results_database_schema.incidence_rate;",
  results_database_schema = resultsSchema,
  snakeCaseToCamelCase = TRUE
) |> dplyr::tibble()
saveRDS(object = incidenceRateIcpe2023, file = "incidenceRateIcpe2023.RDS")


cohortIcpe2023 <- DatabaseConnector::renderTranslateQuerySql(
  connection = connection,
  sql = "SELECT * FROM @results_database_schema.cohort;",
  results_database_schema = resultsSchema,
  snakeCaseToCamelCase = TRUE
) |> dplyr::tibble()
saveRDS(object = cohortIcpe2023, file = "cohortIcpe2023.RDS")

cohortCountIcpe2023 <- DatabaseConnector::renderTranslateQuerySql(
  connection = connection,
  sql = "SELECT * FROM @results_database_schema.cohort_count;",
  results_database_schema = resultsSchema,
  snakeCaseToCamelCase = TRUE
) |> dplyr::tibble()
saveRDS(object = cohortCountIcpe2023, file = "cohortCountIcpe2023.RDS")
