#set up----
setwd(rstudioapi::getActiveDocumentContext()$path |> dirname())
source("functions.R")

#create cohort table
cohort <- dplyr::tibble(cohortId = 1, cohortName = 'Fake cohort')

# scenario 1 - simple stable straight line
fakeDatabaseId1 <- createFakeIncidenceRateData(
  pattern = "stable",
  numberOfYears = 5,
  cohortId = 1,
  databaseId = "fake1",
  incidenceRateMin = 0.0001,
  incidenceRateMax = 0.02,
  cohortCountMin = 5,
  cohortCountMax = 10000
)

predicedIncidenceRateData <- checkTemporalStabilityForcohortDiagnosticsIncidenceRateData(
  cohortDiagnosticsIncidenceRateData = fakeDatabaseId1,
  cohort = cohort,
  maxNumberOfSplines = 5,
  splineTickInterval = 3,
  maxRatio = 1.25,
  alpha = 0.05
) |> dplyr::mutate(cyclopsExpectedIncidenceRate = cyclopsExpectedIncidenceRate / 1000)
plotTemporalTrendExpectedObserved(data = predicedIncidenceRateData)


# scenario 2 - simple stable monotonic line
fakeDatabaseId2 <- createFakeIncidenceRateData(
  pattern = "monotonic",
  numberOfYears = 5,
  cohortId = 1,
  databaseId = "fake1",
  incidenceRateMin = 0.01,
  incidenceRateMax = 0.8,
  cohortCountMin = 5,
  cohortCountMax = 10000
)
predicedIncidenceRateData <- checkTemporalStabilityForcohortDiagnosticsIncidenceRateData(
  cohortDiagnosticsIncidenceRateData = fakeDatabaseId2,
  cohort = cohort,
  maxNumberOfSplines = 5,
  splineTickInterval = 3,
  maxRatio = 1.25,
  alpha = 0.05
) |> dplyr::mutate(cyclopsExpectedIncidenceRate = cyclopsExpectedIncidenceRate / 1000)
plotTemporalTrendExpectedObserved(data = predicedIncidenceRateData)



# scenario 3 - change point
fakeDatabaseId3 <- createFakeIncidenceRateData(
  pattern = "changePoints",
  numberOfYears = 5,
  cohortId = 1,
  databaseId = "fake1",
  incidenceRateMin = 0.0001,
  incidenceRateMax = 0.02,
  cohortCountMin = 5,
  cohortCountMax = 10000
)
predicedIncidenceRateData <- checkTemporalStabilityForcohortDiagnosticsIncidenceRateData(
  cohortDiagnosticsIncidenceRateData = fakeDatabaseId3,
  cohort = cohort,
  maxNumberOfSplines = 5,
  splineTickInterval = 3,
  maxRatio = 1.25,
  alpha = 0.05
) |> dplyr::mutate(cyclopsExpectedIncidenceRate = cyclopsExpectedIncidenceRate / 1000)
plotTemporalTrendExpectedObserved(data = predicedIncidenceRateData)




# scenario 4 - stepwise
fakeDatabaseId4 <- createFakeIncidenceRateData(
  pattern = "stepwise",
  numberOfYears = 5,
  cohortId = 1,
  databaseId = "fake1",
  incidenceRateMin = 0.0001,
  incidenceRateMax = 0.02,
  cohortCountMin = 5,
  cohortCountMax = 10000
)
predicedIncidenceRateData <- checkTemporalStabilityForcohortDiagnosticsIncidenceRateData(
  cohortDiagnosticsIncidenceRateData = fakeDatabaseId4,
  cohort = cohort,
  maxNumberOfSplines = 5,
  splineTickInterval = 1,
  maxRatio = 1.25,
  alpha = 0.05
) |> dplyr::mutate(cyclopsExpectedIncidenceRate = cyclopsExpectedIncidenceRate / 1000)
plotTemporalTrendExpectedObserved(data = predicedIncidenceRateData)

