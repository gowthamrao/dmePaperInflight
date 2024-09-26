#set up----
setwd(rstudioapi::getActiveDocumentContext()$path |> dirname())
source("functions.R")

#get data ----
##icpe2023Data----
cohort <- readRDS("cohortIcpe2023.RDS")
incidenceRate <- readRDS("incidenceRateIcpe2023.RDS")
cohortCount <- readRDS("cohortCountIcpe2023.RDS")

##2024Data-----
cohort <- readRDS("cohort2024.RDS")
incidenceRate <- readRDS("incidenceRate2024.RDS")
cohortCount <- readRDS("cohortCount2024.RDS")


## step 1 - identify outlier years
incidenceRateDataFiltered <- incidenceRate |>
  dplyr::filter(is.na(gender), is.na(ageGroup), !is.na(calendarYear))


#step 1 - find the calendar years to use
checkIfCohortDiagnosticsIncidenceRateDataStratifiedOnlyByCalendarYear(incidenceRateDataFiltered)
calendarYearsToUse <- processCohortDiagnosticsIncidenceRateDataStratifiedOnlyByCalendarYear(
  cohortDiagnosticsIncidenceRateData = incidenceRateDataFiltered,
  removeOutlierZeroCounts = TRUE
) |> 
  dplyr::filter(useYearData == 1) |> 
  dplyr::select(cohortId,
                databaseId,
                calendarYear) |> 
  dplyr::distinct()

cumulativeIncidenceRate <- incidenceRate |>
  dplyr::filter(!is.na(calendarYear)) |>
  dplyr::filter(!is.na(gender), !is.na(ageGroup)) |>
  dplyr::mutate(calendarYear = as.Date(paste0(calendarYear, "-06-01"))) |>
  dplyr::inner_join(calendarYearsToUse,
                    by = c("cohortId", "databaseId", "calendarYear")) |>
  dplyr::group_by(cohortId, databaseId, gender, ageGroup) |>
  dplyr::summarise(
    personYears = sum(abs(ifelse(
      personYears < 0, 1, personYears
    ))),
    cohortCount = sum(abs(ifelse(
      cohortCount < 0, 1, cohortCount
    ))),
    .groups = "keep"
  ) |>
  dplyr::ungroup() |>
  dplyr::rename(outcomes = cohortCount)

cumulativeIncidenceRate <- cumulativeIncidenceRate |>
  dplyr::mutate(sourceKey = databaseId) |>
  dplyr::mutate(
    databaseId = stringr::str_replace_all(
      string = sourceKey,
      pattern = "^cdm_|_v[0-9]+$",
      replacement = ""
    )
  ) |>
  dplyr::mutate(databaseId = SqlRender::camelCaseToTitleCase(SqlRender::snakeCaseToCamelCase(databaseId))) |> 
  dplyr::mutate(incidenceRate = (outcomes/personYears)*1000) |> 
  dplyr::mutate(incidenceRateSe = sqrt(incidenceRate / personYears)) # Standard error for the incidence rate

library(meta)
library(dplyr)

# Define the function to extract meta-analysis results
rand_te <- function(data, .group) {
  random.meta <- meta::metarate(
    data = data,
    event = data$outcomes,
    time = data$personYears,
    studlab = data$databaseId,
    sm = "IRLN",
    comb.random = TRUE,
    method.tau = "DL"
  )
  
  # Extract relevant statistics from the meta-analysis result
  random.te <- random.meta[["TE.random"]]
  random.te.lower <- random.meta[["lower.random"]]
  random.te.upper <- random.meta[["upper.random"]]
  seTE.random <- random.meta[["seTE.random"]]
  tau2 <- random.meta[["tau2"]]
  se.tau2 <- random.meta[["se.tau2"]]
  tau <- random.meta[["tau"]]
  lower.predict <- random.meta[["lower.predict"]]
  upper.predict <- random.meta[["upper.predict"]]
  seTE.predict <- random.meta[["seTE.predict"]]
  
  # Create a data frame with the meta-analysis results
  meta.out <-
    data.frame(
      random.te,
      random.te.lower,
      random.te.upper,
      seTE.random,
      tau2,
      se.tau2,
      tau,
      lower.predict,
      upper.predict,
      seTE.predict
    ) |>
    dplyr::mutate(
      ir.rand = exp(random.te) * 1000,
      ir.rand.l = exp(random.te.lower) * 1000,
      ir.rand.u = exp(random.te.upper) * 1000,
      ir.predict.lower = exp(lower.predict) * 1000,
      ir.predict.upper = exp(upper.predict) * 1000,
      ir.predict.seTE = exp(seTE.predict) * 1000 
    ) |>
    dplyr::tibble()
  
  return(meta.out)
}


metaAnalyticEstimate <- c()
databaseIds <- cumulativeIncidenceRate$databaseId |> unique() |> sort()


for (i in (1:length(databaseIds))) {
  targetDatabaseId <- databaseIds[[i]]
  allOtherDatabaseId <- databaseIds[-i] |> sort()
  
  dataToPredict <- cumulativeIncidenceRate |>
    dplyr::filter(databaseId %in% allOtherDatabaseId)
  
  metaAnalyticEstimate[[targetDatabaseId]] <- dataToPredict |>
    dplyr::inner_join(
      cumulativeIncidenceRate |>
        dplyr::filter(databaseId == targetDatabaseId) |>
        dplyr::select(ageGroup, gender, cohortId) |>
        dplyr::distinct()
    ) |>
    dplyr::group_by(cohortId, gender, ageGroup) |>
    dplyr::group_modify( ~ rand_te(.x, .y)) |>
    dplyr::mutate(
      databaseIdTarget = targetDatabaseId,
      databaseIdComparator = paste0(allOtherDatabaseId, collapse = ", ")
    )
}


dplyr::bind_rows(metaAnalyticEstimate) |> saveRDS("metaAnalyticEstimate.RDS")
metaAnalyticEstimate <- readRDS("metaAnalyticEstimate.RDS")




output <- metaAnalyticEstimate |>
  dplyr::left_join(cohort |>
                     dplyr::select(cohortId, cohortName)) |>
  dplyr::left_join(
    cumulativeIncidenceRate |>
      dplyr::rename(databaseIdTarget = databaseId),
    by = c("cohortId", "databaseIdTarget", "gender", "ageGroup")
  ) |> dplyr::mutate(
    isWithinPredictionInterval = ifelse(
      incidenceRate >= ir.predict.lower &
        incidenceRate <= ir.predict.upper,
      TRUE,
      FALSE
    )
  ) |>
  dplyr::mutate(pValueRaw = 2 * (1 - pnorm(abs((incidenceRate - (ir.predict.lower + ir.predict.upper) / 2) / incidenceRateSe
  )))) |>
  dplyr::mutate(pValue = ifelse(isWithinPredictionInterval, 1, 2 * (1 - pnorm(
    abs((
      incidenceRate - (ir.predict.lower + ir.predict.upper) / 2
    ) / incidenceRateSe)
  )))) |>
  dplyr::select(
    cohortId,
    cohortName,
    databaseIdTarget,
    databaseIdComparator,
    outcomes,
    personYears,
    ir.predict.lower,
    ir.predict.upper,
    incidenceRate,
    pValue,
    pValueRaw
  ) |>
  dplyr::rename(cohortCount = outcomes) |>
  dplyr::mutate(
    cohortCount = OhdsiHelpers::formatIntegerWithComma(cohortCount),
    personYears = OhdsiHelpers::formatIntegerWithComma(personYears),
    ir.predict.lower = OhdsiHelpers::formatDecimalWithComma(ir.predict.lower, decimalPlaces = 6),
    ir.predict.upper = OhdsiHelpers::formatDecimalWithComma(ir.predict.upper, decimalPlaces = 6),
    incidenceRate = OhdsiHelpers::formatDecimalWithComma(incidenceRate, decimalPlaces = 6),
    pValue = OhdsiHelpers::formatDecimalWithComma(pValue, decimalPlaces = 3),
    pValueRaw = OhdsiHelpers::formatDecimalWithComma(pValueRaw, decimalPlaces = 3)
  ) |>
  dplyr::mutate(reject = ifelse(
    pValue < 0.05,
    1,
    0
  ))

output |> 
  dplyr::group_by(cohortId,cohortName, databaseIdTarget, reject)
