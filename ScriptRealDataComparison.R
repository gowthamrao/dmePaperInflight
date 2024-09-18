#set up----
setwd(rstudioapi::getActiveDocumentContext()$path |> dirname())
source("functions.R")

#get data ----
##icpe2023Data----
cohort <- readRDS("cohortIcpe2023.RDS")
incidenceRate <- readRDS("incidenceRateIcpe2023.RDS") |>
  dplyr::filter(gender != '', ageGroup != '', calendarYear != '')
cohortCount <- readRDS("cohortCountIcpe2023.RDS")

##2024Data-----
cohort <- readRDS("cohort2024.RDS")
incidenceRate <- readRDS("incidenceRate2024.RDS")  |>
  dplyr::filter(!is.na(gender), !is.na(ageGroup), !is.na(calendarYear)) |>
  dplyr::filter(cohortCount > 0) |>
  dplyr::mutate(
    cohortCount = abs(cohortCount),
    incidenceRate = abs(cohortCount) / abs(personYears)
  )

incidenceRate <- incidenceRate |>
  dplyr::left_join(referencePopulation, by = c("gender", "ageGroup")) |>
  dplyr::left_join(
    incidenceRate |>
      dplyr::group_by(cohortId, databaseId, calendarYear) |>
      dplyr::summarise(
        totalPersonYears = sum(personYears),
        .groups = "keep"
      ) |>
      dplyr::ungroup(),
    by = c("cohortId", "databaseId", "calendarYear")
  ) |>
  dplyr::mutate(
    standardizedPersonYears = totalPersonYears * proportion,
    standardizedCohortCount = incidenceRate * standardizedPersonYears
  ) |>
  dplyr::group_by(cohortId, databaseId, calendarYear) |>
  dplyr::summarise(
    personYears = sum(personYears),
    cohortCount = sum(cohortCount),
    standardizedPersonYears = sum(standardizedPersonYears),
    standardizedCohortCount = sum(standardizedCohortCount),
    .groups = "keep"
  ) |>
  dplyr::ungroup() |>
  dplyr::arrange(cohortId, databaseId, calendarYear) |>
  dplyr::mutate(
    cohortCount = round(cohortCount),
    personYears = round(personYears),
    standardizedCohortCount = round(standardizedCohortCount),
    standardizedPersonYears = round(standardizedPersonYears)
  ) |>
  dplyr::mutate(
    incidenceRate = (cohortCount / personYears) * 1000,
    standardizedIncidenceRate = (standardizedCohortCount / standardizedPersonYears) * 1000
  ) |>
  dplyr::mutate(ratioOfStandardizedToReportedIncidenceRate = incidenceRate /
                  standardizedIncidenceRate)

#R2: stability----
temporalStabilityOutputRawData <- checkTemporalStabilityForcohortDiagnosticsIncidenceRateData(
  cohortDiagnosticsIncidenceRateData = incidenceRate,
  cohort = cohort,
  removeOutlierZeroCounts = TRUE,
  maxNumberOfSplines = 3,
  splineTickInterval = 3,
  maxRatio = 1.25,
  alpha = 0.05
)

standardizedIncidenceRate <-  incidenceRate |>
  dplyr::select(-incidenceRate, -cohortCount, -personYears) |>
  dplyr::rename(
    cohortCount = standardizedCohortCount,
    personYears = standardizedPersonYears,
    incidenceRate = standardizedIncidenceRate
  )

temporalStabilityOutputStandardizedData <- checkTemporalStabilityForcohortDiagnosticsIncidenceRateData(
  cohortDiagnosticsIncidenceRateData = standardizedIncidenceRate,
  cohort = cohort,
  removeOutlierZeroCounts = TRUE,
  maxNumberOfSplines = 3,
  splineTickInterval = 3,
  maxRatio = 1.25,
  alpha = 0.05
)

#statistically unstable----
cylopsStatisticallyUnstableRawData <- temporalStabilityOutputRawData |>
  dplyr::filter(cyclopsStable == FALSE) |>
  dplyr::select(cohortId,
                cohortName,
                databaseId,
                cyclopsStable,
                cyclopsRatio,
                cyclopsPValue) |>
  dplyr::distinct() |>
  dplyr::mutate(
    cyclopsRatio = OhdsiHelpers::formatDecimalWithComma(cyclopsRatio, decimalPlaces = 5),
    cyclopsPValue = OhdsiHelpers::formatDecimalWithComma(cyclopsPValue, decimalPlaces = 5)
  ) |>
  dplyr::arrange(cohortId, databaseId)
cylopsStatisticallyUnstableRawData


cylopsStatisticallyUnstableStandardizedData <- temporalStabilityOutputStandardizedData |>
  dplyr::filter(cyclopsStable == FALSE) |>
  dplyr::select(cohortId,
                cohortName,
                databaseId,
                cyclopsStable,
                cyclopsRatio,
                cyclopsPValue) |>
  dplyr::distinct() |>
  dplyr::mutate(
    cyclopsRatio = OhdsiHelpers::formatDecimalWithComma(cyclopsRatio, decimalPlaces = 5),
    cyclopsPValue = OhdsiHelpers::formatDecimalWithComma(cyclopsPValue, decimalPlaces = 5)
  ) |>
  dplyr::arrange(cohortId, databaseId)
cylopsStatisticallyUnstableStandardizedData

#R3: pass or fail rate-----
summarizeTemporalStability(temporalStabilityOutputRawData,
                           cylopsStatisticallyUnstableRawData)
summarizeTemporalStability(
  temporalStabilityOutputStandardizedData,
  cylopsStatisticallyUnstableStandardizedData
)


incidenceRatesComparison <- temporalStabilityOutputRawData |> dplyr::select(cohortId, databaseId, cohortName, calendarYear, incidenceRate) |>
  dplyr::rename(incidenceRateRaw = incidenceRate) |>
  dplyr::left_join(
    temporalStabilityOutputStandardizedData |> dplyr::select(cohortId, databaseId, calendarYear, incidenceRate) |>
      dplyr::rename(incidenceRateStandardized = incidenceRate),
    by = c("cohortId", "databaseId", "calendarYear")
  ) |>
  dplyr::left_join(
    readRDS("incidenceRate2024.RDS")  |>
      dplyr::filter(is.na(gender), is.na(ageGroup), !is.na(calendarYear)) |>
      dplyr::select(cohortId, databaseId, calendarYear, incidenceRate) |>
      dplyr::rename(incidenceRateOriginal = incidenceRate) |>
      dplyr::mutate(calendarYear = as.Date(paste0(
        calendarYear, "-06-01"
      ))),
    by = c("cohortId", "databaseId", "calendarYear")
  ) |>
  dplyr::mutate(
    originalToStandardized = incidenceRateStandardized / incidenceRateOriginal,
    originalToSumFromStrata = incidenceRateRaw / incidenceRateOriginal
  )

#clear rstudio plots
if (!is.null(dev.list())) {
  dev.off()
}


#statistic fails
plotsPureRedCellAplasiaRawData <- createPlotsByDatabaseId(data = temporalStabilityOutputRawData, cohortId = 207)
plotsPureRedCellAplasiaStandardizedData <- createPlotsByDatabaseId(data = temporalStabilityOutputStandardizedData, cohortId = 207)
plotsPureRedCellAplasiaRawData
plotsPureRedCellAplasiaStandardizedData

plotsPolyMorphicVentricularTachycardiaRawData <- createPlotsByDatabaseId(data = temporalStabilityOutputRawData, cohortId = 275)
plotsPolyMorphicVentricularTachycardiaStandardizedData <- createPlotsByDatabaseId(data = temporalStabilityOutputStandardizedData, cohortId = 275)
plotsPolyMorphicVentricularTachycardiaRawData
plotsPolyMorphicVentricularTachycardiaStandardizedData

plotsAcuteHepaticFailureRawData <- createPlotsByDatabaseId(data = temporalStabilityOutputRawData, cohortId = 207)
plotsAcuteHepaticFailureStandardizedData <- createPlotsByDatabaseId(data = temporalStabilityOutputStandardizedData, cohortId = 207)
plotsAcuteHepaticFailureRawData
plotsAcuteHepaticFailureStandardizedData

plotsAutoimmuneHepatitisRawData <- createPlotsByDatabaseId(data = temporalStabilityOutputRawData, cohortId = 729)
plotsAutoimmuneHepatitisStandardizedData <- createPlotsByDatabaseId(data = temporalStabilityOutputStandardizedData, cohortId = 729)
plotsAutoimmuneHepatitisRawData
plotsAutoimmuneHepatitisStandardizedData

plotsAutoimmuneHemolyticAnemiaRawData <- createPlotsByDatabaseId(data = temporalStabilityOutputRawData, cohortId = 728)
plotsAutoimmuneHemolyticAnemiaStandardizedData <- createPlotsByDatabaseId(data = temporalStabilityOutputStandardizedData, cohortId = 728)
plotsAutoimmuneHemolyticAnemiaRawData
plotsAutoimmuneHemolyticAnemiaStandardizedData

plotsDressRawData <- createPlotsByDatabaseId(data = temporalStabilityOutputRawData, cohortId = 733)
plotsDressStandardizedData <- createPlotsByDatabaseId(data = temporalStabilityOutputStandardizedData, cohortId = 733)
plotsDressRawData
plotsDressStandardizedData

plotsAliRawData <- createPlotsByDatabaseId(data = temporalStabilityOutputRawData, cohortId = 736)
plotsAliStandardizedData <- createPlotsByDatabaseId(data = temporalStabilityOutputStandardizedData, cohortId = 736)
plotsAliRawData
plotsAliStandardizedData






## comparative
cohort729StandardizedData <- temporalStabilityOutputStandardizedData |>
  dplyr::filter(cohortId == 729) |>
  dplyr::mutate(expected = dplyr::coalesce(glmExpected, cyclopsExpected),
                observed = cohortCount) |>
  dplyr::mutate(expectedObservedRatio = expected / observed) |>
  dplyr::select(databaseId, calendarYear, expectedObservedRatio) |>
  dplyr::filter(!is.na(expectedObservedRatio) &
                  is.finite(expectedObservedRatio))

performGAMAnalysis(cohort729StandardizedData)
