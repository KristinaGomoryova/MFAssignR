load_expected <- function(mode) {
  mono <- readRDS(file.path("test-data", paste0(mode, "_recal_mono.rds")))
  iso <- readRDS(file.path("test-data", paste0(mode, "_recal_iso.rds")))
  recallist <- readRDS(file.path("test-data", paste0(mode, "_recal_recallist.rds")))
  return(list(Mono = mono, Iso = iso, RecalList = recallist))
}

update_expected <- function(actual, mode) {
  saveRDS(actual$Mono, file.path("test-data", paste0(mode, "_recal_mono.rds")))
  saveRDS(actual$Iso, file.path("test-data", paste0(mode, "_recal_iso.rds")))
  saveRDS(actual$RecalList, file.path("test-data", paste0(mode, "_recal_recallist.rds")))
}

patrick::with_parameters_test_that("Recal works",
  {
    peaks <- readRDS(file.path("test-data", paste0(mode, "_iso.rds")))
    unambig <- readRDS(file.path("test-data", paste0(mode, "_cho_unambig.rds")))
    recallist <- readRDS(file.path("test-data", paste0(mode, "_recallist.rds"))) |> dplyr::arrange_at("Series.Score")

    actual <- MFAssignR::Recal(
      unambig,
      peaks = peaks$Mono,
      isopeaks = peaks$Iso,
      mzRange = 30,
      mode = mode,
      SN = 0,
      series1 = recallist$Series[1],
      series2 = recallist$Series[2],
      series3 = recallist$Series[3],
      series4 = recallist$Series[4],
      series5 = recallist$Series[5]
    )

    expected <- load_expected(mode)

    expect_equal(actual$Mono, expected$Mono)
    expect_equal(actual$Iso, expected$Iso)
    expect_equal(actual$RecalList, expected$RecalList)
  },
  mode = c("neg", "pos"),
)

patrick::with_parameters_test_that("RecalList works",
  {
    expected_path <- file.path("test-data", paste0(mode, "_recallist.rds"))
    unambig <- readRDS(file.path("test-data", paste0(mode, "_cho_unambig.rds")))

    actual <- MFAssignR::RecalList(unambig)
    expected <- readRDS(expected_path)

    actual <- actual %>% dplyr::select(-"Series.Index")
    expected <- expected %>% dplyr::select(-"Series.Index")
    actual_sorted <- dplyr::arrange_at(actual, "Series")
    expected_sorted <- dplyr::arrange_at(expected, "Series")
    
    expect_equal(actual_sorted, expected_sorted)
  },
  mode = c("neg", "pos"),
)

test_that("Replicate the Recal bug", {
  unambig <- read.delim("test-data/bug_unambig.tabular")
  mono <- read.delim("test-data/bug_mono.tabular")
  iso <- read.delim("test-data/bug_iso.tabular")
  recallist <- read.delim("test-data/bug_recalseries.tabular")

  actual <- MFAssignR::Recal(
      df = unambig,
      peaks = mono,
      isopeaks = iso,
      mode = "neg",
      SN = 6*346.0706,
      series1 = recallist$Series[1],
      series2 = recallist$Series[2],
      series3 = recallist$Series[3],
      series4 = recallist$Series[4],
      series5 = recallist$Series[5],
      step_O = 3,
      step_H2 = 5,
      mzRange = 50,
      CalPeak = 150
    )
    expected <- readRDS("test-data/bug_recal_expected.rds")
    expect_equal(actual, expected)
})

patrick::with_parameters_test_that("Replicate Recal isopeaks error", {
  mode <- "neg"
  peaks <- readRDS(file.path("test-data", paste0(mode, "_iso.rds")))
  unambig <- readRDS(file.path("test-data", paste0(mode, "_cho_unambig.rds")))
  recallist <- readRDS(file.path("test-data", paste0(mode, "_recallist.rds"))) |> dplyr::arrange_at("Series.Score")

  actual <- MFAssignR::Recal(
      df = unambig,
      peaks = peaks$Mono,
      isopeaks = iso,
      mode = mode,
      series1 = recallist$Series[1],
      series2 = recallist$Series[2],
      series3 = recallist$Series[3],
      series4 = recallist$Series[4],
      series5 = recallist$Series[5],
      mzRange = 80
    )
    expected <- readRDS("test-data/recal_isopeaks.rds")

    expect_equal(actual, expected)
},
iso = c(NA, 
        data.frame(exp_mass = NA, abundance = NA, tag = NA), 
        data.frame(exp_mass = c(NA, NA), abundance = c(NA, NA), tag = c(NA, NA)))
)

test_that("Column names in the Recal function", {
    mode = "pos"
    unambig <- readRDS(file.path("test-data", paste0(mode, "_cho_unambig.rds")))   
    expected <- c("Series", "Number.Observed", "Series.Index", "Mass.Range", "Tall.Peak", "Abundance.Score", "Peak.Score", "Peak.Distance", "Series.Score")

    actual <- colnames(MFAssignR::RecalList(unambig))

    expect_equal(expected, actual)
})