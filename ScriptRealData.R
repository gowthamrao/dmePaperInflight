#set up----
setwd(rstudioapi::getActiveDocumentContext()$path |> dirname())
source("functions.R")

#get data ----
##icpe2023Data----
cohort <- readRDS("cohortIcpe2023.RDS")
incidenceRate <- readRDS("incidenceRateIcpe2023.RDS") |>
  dplyr::filter(gender == '', ageGroup == '', calendarYear != '')
cohortCount <- readRDS("cohortCountIcpe2023.RDS")

##2024Data-----
cohort <- readRDS("cohort2024.RDS")
incidenceRate <- readRDS("incidenceRate2024.RDS") |>
  dplyr::filter(is.na(gender), is.na(ageGroup), !is.na(calendarYear))
cohortCount <- readRDS("cohortCount2024.RDS")

#R1: cohort counts----
reportCohortCount <- cohortCount |>
  dplyr::arrange(cohortId, dplyr::desc(cohortSubjects)) |>
  tidyr::pivot_wider(
    id_cols = c(cohortId),
    values_from = cohortSubjects,
    names_from = databaseId
  ) |>
  dplyr::inner_join(cohort |>
                      dplyr::select(cohortId, cohortName), by = "cohortId") |>
  dplyr::relocate(cohortId, cohortName)

#R2: stability----
temporalStabilityOutput <- checkTemporalStabilityForcohortDiagnosticsIncidenceRateData(
  cohortDiagnosticsIncidenceRateData = incidenceRate,
  cohort = cohort,
  removeOutlierZeroCounts = TRUE,
  maxNumberOfSplines = 3,
  splineTickInterval = 3,
  maxRatio = 1.25,
  alpha = 0.05
)

#statistically unstable----
cylopsStatisticallyUnstable <- temporalStabilityOutput |>
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
cylopsStatisticallyUnstable

# #statistically unstable----
glmStatisticallyUnstable <- temporalStabilityOutput |>
  dplyr::filter(glmStable == FALSE) |>
  dplyr::select(cohortId,
                cohortName,
                databaseId,
                glmStable,
                glmRatio,
                glmPValue,
                glmPPValuePearson) |>
  dplyr::distinct() |>
  dplyr::mutate(glmPPValuePearson = OhdsiHelpers::formatDecimalWithComma(glmPPValuePearson, decimalPlaces = 5)) |>
  dplyr::arrange(cohortId, databaseId)
glmStatisticallyUnstable


#R3: pass or fail rate-----
summarizeTemporalStability(temporalStabilityOutput, cylopsStatisticallyUnstable)
summarizeTemporalStability(temporalStabilityOutput, glmStatisticallyUnstable)

#clear rstudio plots
if (!is.null(dev.list())) {
  dev.off()
}

debug(plotTemporalTrendExpectedObserved)
#statistic fails
plotsPureRedCellAplasia <- createPlotsByDatabaseId(data = temporalStabilityOutput, cohortId = 207)
plotsPureRedCellAplasia

plotsPolyMorphicVentricularTachycardia <- createPlotsByDatabaseId(data = temporalStabilityOutput, cohortId = 275)
plotsPolyMorphicVentricularTachycardia

plotFirstAcuteHepaticFailure <- createPlotsByDatabaseId(data = temporalStabilityOutput, cohortId = 724)
plotFirstAcuteHepaticFailure

plotAutoimmuneHepatitis <- createPlotsByDatabaseId(data = temporalStabilityOutput, cohortId = 729)
plotAutoimmuneHepatitis

plotAutoimmuneHemolyticAnemia <- createPlotsByDatabaseId(data = temporalStabilityOutput, cohortId = 728)
plotAutoimmuneHemolyticAnemia

plotsDress <- createPlotsByDatabaseId(data = temporalStabilityOutput, cohortId = 733)
plotsDress

plotsAli <- createPlotsByDatabaseId(data = temporalStabilityOutput, cohortId = 736)
plotsAli
