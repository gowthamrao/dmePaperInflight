#set up----
setwd(rstudioapi::getActiveDocumentContext()$path |> dirname())
source("functions.R")

#icpe2023Data
cohort <- readRDS("cohortIcpe2023.rds")
incidenceRate <- readRDS("incidenceRateIcpe2023.rds") |>
  dplyr::filter(gender == '', ageGroup == '', calendarYear != '')
cohortCount <- readRDS("cohortCountIcpe2023.rds")

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
  maxNumberOfSplines = 5,
  splineTickInterval = 3,
  maxRatio = 1.25,
  alpha = 0.05
)

#statistically unstable----
statisticallyUnstable <- temporalStabilityOutput |>
  dplyr::filter(stable == FALSE) |>
  dplyr::select(cohortId, cohortName, databaseId, stable, ratio, p) |>
  dplyr::distinct() |>
  dplyr::mutate(
    ratio = OhdsiHelpers::formatDecimalWithComma(ratio, decimalPlaces = 5),
    p = OhdsiHelpers::formatDecimalWithComma(p, decimalPlaces = 5)
  )
statisticallyUnstable

plots <- createPlotsByDatabaseId(data = temporalStabilityOutput, cohortId = 723)
plots
