skip_map <- function() {
  if (!exists(".map_aggregated_columns", mode = "function"))
    skip("dge.R not sourced")
}

meta_fixture <- function(conditions, samples) {
  data.frame(
    cond = rep(conditions, each = length(samples)),
    samp = rep(samples, times = length(conditions)),
    stringsAsFactors = FALSE
  )
}

test_that("underscored condition names map correctly", {
  skip_map()
  meta <- meta_fixture(c("Severe_DISEASE", "Healthy"), c("P001", "P002"))
  cols <- c("Severe_DISEASE_P001", "Severe_DISEASE_P002", "Healthy_P001", "Healthy_P002")
  out  <- .map_aggregated_columns(cols, meta, c("cond", "samp"))

  expect_equal(nrow(out), 4L)

  expect_equal(out[["cond"]], c("Severe_DISEASE", "Severe_DISEASE", "Healthy", "Healthy"))
  expect_equal(out[["samp"]], c("P001", "P002", "P001", "P002"))
})

test_that("underscores in the SECOND group are handled too", {
  skip_map()

  meta <- meta_fixture(c("Healthy"), c("P001_v2", "P002_v2"))
  cols <- c("Healthy_P001_v2", "Healthy_P002_v2")
  out  <- .map_aggregated_columns(cols, meta, c("cond", "samp"))
  expect_equal(out[["cond"]], c("Healthy", "Healthy"))
  expect_equal(out[["samp"]], c("P001_v2", "P002_v2"))
})

test_that("Seurat's g-prefix for numeric levels is recovered", {
  skip_map()
  meta <- meta_fixture(c("3", "12"), c("40552"))
  out  <- .map_aggregated_columns(c("g3_g40552", "g12_g40552"), meta, c("cond", "samp"))
  expect_equal(out[["cond"]], c("3", "12"))
  expect_equal(out[["samp"]], c("40552", "40552"))
})

test_that("an unmappable column is a hard error, not a silent guess", {
  skip_map()
  meta <- meta_fixture(c("Healthy"), c("P001"))
  expect_error(
    .map_aggregated_columns(c("Healthy_P001", "Totally_Unknown_Column"),
                            meta, c("cond", "samp")),
    "Could not map"
  )
})

test_that("three grouping columns still resolve", {
  skip_map()
  meta <- data.frame(
    cond = c("Severe_DISEASE", "Healthy"),
    tp   = c("Day_0", "Day_0"),
    samp = c("P001", "P002"),
    stringsAsFactors = FALSE
  )
  cols <- c("Severe_DISEASE_Day_0_P001", "Healthy_Day_0_P002")
  out  <- .map_aggregated_columns(cols, meta, c("cond", "tp", "samp"))
  expect_equal(out[["cond"]], c("Severe_DISEASE", "Healthy"))
  expect_equal(out[["tp"]],   c("Day_0", "Day_0"))
  expect_equal(out[["samp"]], c("P001", "P002"))
})

test_that("COMBAT Source levels survive aggregation round-trip", {
  skip_map()
  source_levels <- c("HV", "DISEASE_MILD", "DISEASE_SEV", "DISEASE_CRIT",
                     "DISEASE_HCW_MILD", "DISEASE_LDN", "Sepsis", "Flu")
  ids <- c("G05061", "G05064", "G05073")

  meta <- expand.grid(cond = source_levels, samp = ids,
                      stringsAsFactors = FALSE)
  cols <- paste(meta$cond, meta$samp, sep = "_")

  out <- .map_aggregated_columns(cols, meta, c("cond", "samp"))

  expect_equal(out[["cond"]], meta$cond)
  expect_equal(out[["samp"]], meta$samp)

  expect_equal(sort(unique(out[["cond"]])), sort(source_levels))
  expect_equal(length(unique(out[["cond"]])), 8L)
})

test_that("Yoshida batch ids with underscores round-trip", {
  skip_map()
  meta <- data.frame(
    cond = c("DISEASE+", "Healthy", "Post-DISEASE"),
    samp = c("CV001_KM8853544", "CV001_KM8853545", "CV001_KM8853546"),
    stringsAsFactors = FALSE
  )
  cols <- paste(meta$cond, meta$samp, sep = "_")
  out  <- .map_aggregated_columns(cols, meta, c("cond", "samp"))
  expect_equal(out[["cond"]], meta$cond)
  expect_equal(out[["samp"]], meta$samp)
})

test_that("digit-leading levels resolve under both g- and X-prefixing", {
  skip_map()
  meta <- data.frame(
    Terminal_Group   = c("24mO", "24mR", "3mR", "12mR", "C24mO"),
    Participant_ID   = c("40035", "40093", "40186", "40552", "40865"),
    Timepoint_months = c("24", "12", "3", "1", "24"),
    stringsAsFactors = FALSE
  )
  gp <- c("Terminal_Group", "Participant_ID", "Timepoint_months")

  raw <- apply(meta[, gp], 1, paste, collapse = "_")
  gpx <- apply(meta[, gp], 1, function(r)
    paste(ifelse(grepl("^[0-9]", r), paste0("g", r), r), collapse = "_"))
  xpx <- apply(meta[, gp], 1, function(r) paste(make.names(r), collapse = "_"))

  for (cols in list(raw, gpx, xpx)) {
    out <- .map_aggregated_columns(cols, meta, gp)
    expect_equal(out[["Terminal_Group"]],   meta$Terminal_Group)
    expect_equal(out[["Participant_ID"]],   meta$Participant_ID)
    expect_equal(out[["Timepoint_months"]], meta$Timepoint_months)
  }
})
