#create cohort table
cohort <- dplyr::tibble(cohortId = 1, cohortName = 'Fake cohort')


# scenario 1 - simple stable straight line
fakeDatabaseId <- createFakeIncidenceRateData(
  pattern = "stable",
  numberOfYears = 5,
  cohortId = 1,
  databaseId = "fake1",
  incidenceRateMin = 0.0001,
  incidenceRateMax = 0.02,
  cohortCountMin = 5,
  cohortCountMax = 10000
)
plotSimpleTemporalTrend(data = fakeDatabaseId)
predicedIncidenceRateData <- checkTemporalStabilityForcohortDiagnosticsIncidenceRateData(
  cohortDiagnosticsIncidenceRateData = fakeDatabaseId,
  cohort = cohort,
  maxNumberOfSplines = 5,
  splineTickInterval = 3,
  maxRatio = 1.25,
  alpha = 0.05
)
plotTemporalTrendExpectedObserved(data = predicedIncidenceRateData)



# scenario 2 - simple stable monotonic line
fakeDatabaseId <- createFakeIncidenceRateData(
  pattern = "monotonic",
  numberOfYears = 5,
  cohortId = 1,
  databaseId = "fake1",
  incidenceRateMin = 0.01,
  incidenceRateMax = 0.8,
  cohortCountMin = 5,
  cohortCountMax = 10000
)
plotSimpleTemporalTrend(data = fakeDatabaseId)
predicedIncidenceRateData <- checkTemporalStabilityForcohortDiagnosticsIncidenceRateData(
  cohortDiagnosticsIncidenceRateData = fakeDatabaseId,
  cohort = cohort,
  maxNumberOfSplines = 5,
  splineTickInterval = 3,
  maxRatio = 1.25,
  alpha = 0.05
)
plotTemporalTrendExpectedObserved(data = predicedIncidenceRateData)



# scenario 3 - change point
fakeDatabaseId <- createFakeIncidenceRateData(
  pattern = "changePoints",
  numberOfYears = 5,
  cohortId = 1,
  databaseId = "fake1",
  incidenceRateMin = 0.0001,
  incidenceRateMax = 0.02,
  cohortCountMin = 5,
  cohortCountMax = 10000
)
plotSimpleTemporalTrend(data = fakeDatabaseId)
predicedIncidenceRateData <- checkTemporalStabilityForcohortDiagnosticsIncidenceRateData(
  cohortDiagnosticsIncidenceRateData = fakeDatabaseId,
  cohort = cohort,
  maxNumberOfSplines = 5,
  splineTickInterval = 3,
  maxRatio = 1.25,
  alpha = 0.05
)
plotTemporalTrendExpectedObserved(data = predicedIncidenceRateData)




# scenario 4 - stepwise
fakeDatabaseId <- createFakeIncidenceRateData(
  pattern = "stepwise",
  numberOfYears = 5,
  cohortId = 1,
  databaseId = "fake1",
  incidenceRateMin = 0.0001,
  incidenceRateMax = 0.02,
  cohortCountMin = 5,
  cohortCountMax = 10000
)
plotSimpleTemporalTrend(data = fakeDatabaseId)
predicedIncidenceRateData <- checkTemporalStabilityForcohortDiagnosticsIncidenceRateData(
  cohortDiagnosticsIncidenceRateData = fakeDatabaseId,
  cohort = cohort,
  maxNumberOfSplines = 5,
  splineTickInterval = 1,
  maxRatio = 1.25,
  alpha = 0.05
)
plotTemporalTrendExpectedObserved(data = predicedIncidenceRateData)


