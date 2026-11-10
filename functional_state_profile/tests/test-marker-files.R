.find_shipped <- function(fn) {
  cands <- c(fn, file.path("..", fn),
             file.path(dirname(dirname(FIXTURES)), fn),
             file.path(dirname(FIXTURES), "..", fn))
  for (p in cands) if (file.exists(p)) return(normalizePath(p))
  stop("cannot locate shipped file: ", fn)
}
RNA_MARKERS  <- .find_shipped("cell_state_markers.json")
ADT_MARKERS  <- .find_shipped("cell_state_markers_adt.json")

read_states <- function(path) {
  d <- jsonlite::fromJSON(path, simplifyVector = TRUE)
  d[!startsWith(names(d), "_")]
}

test_that("both marker files parse and every state has a non-empty core", {
  for (path in c(RNA_MARKERS, ADT_MARKERS)) {
    states <- read_states(path)
    expect_gt(length(states), 0)
    for (st in names(states)) {
      expect_true(length(states[[st]]$core) > 0,
                  info = sprintf("%s: state '%s' has an empty core", basename(path), st))
    }
  }
})

test_that("no state carries duplicate symbols across core + additional", {
  for (path in c(RNA_MARKERS, ADT_MARKERS)) {
    states <- read_states(path)
    for (st in names(states)) {
      genes <- c(states[[st]]$core, states[[st]]$additional)
      dups  <- unique(genes[duplicated(genes)])
      expect_length(dups, 0)
      if (length(dups) > 0) {
        fail(sprintf("%s: state '%s' duplicates: %s",
                     basename(path), st, paste(dups, collapse = ", ")))
      }
    }
  }
})

test_that("RNA marker symbols are plausibly HGNC", {

  not_hgnc <- c("CD25", "CD127", "CD62L", "CD95", "CD122", "CD137", "CD57",
                "CD16", "CD8", "CD103", "CD49A", "CD49D", "CD45RA", "CD45RO",
                "HLA-DR", "PD-1", "PD1", "TIM3", "TIM-3", "CTLA-4", "OX40",
                "GITR", "4-1BB", "LAMP-1", "CD107A", "FR4", "CD11B", "CD94",
                "NKG2D", "CD31", "CD71", "HELIOS", "FOXP-3")

  states <- read_states(RNA_MARKERS)
  for (st in names(states)) {
    genes <- c(states[[st]]$core, states[[st]]$additional)
    bad_alias <- intersect(toupper(genes), not_hgnc)
    expect_length(bad_alias, 0)
    if (length(bad_alias) > 0) {
      fail(sprintf("RNA state '%s' contains protein alias names: %s",
                   st, paste(bad_alias, collapse = ", ")))
    }

    ok <- grepl("^[A-Z][A-Z0-9]*(-[A-Z0-9]+)*$", genes)
    expect_true(all(ok),
                info = sprintf("RNA state '%s' has malformed symbols: %s",
                               st, paste(genes[!ok], collapse = ", ")))
  }
})

test_that("RNA file is free of housekeeping-family genes", {

  fam_regex <- c("^RP[LS]", "^MT-", "^YWHA", "^TUB[AB]", "^HIST")
  fam_exact <- c("UBB", "UBC", "UBA52", "RPS27A", "ACTB", "ACTG1", "GAPDH",
                 "B2M", "PTPRC", "HMGB1", "HSPA8", "EEF1A1", "EEF2",
                 "TMSB4X", "TMSB10", "MALAT1", "FTL", "FTH1", "PFN1", "CFL1")
  states <- read_states(RNA_MARKERS)
  for (st in names(states)) {
    genes <- c(states[[st]]$core, states[[st]]$additional)
    hits  <- genes[genes %in% fam_exact |
                   Reduce(`|`, lapply(fam_regex, function(p) grepl(p, genes)))]
    expect_length(hits, 0)
    if (length(hits) > 0) {
      fail(sprintf("RNA state '%s' contains housekeeping genes: %s",
                   st, paste(hits, collapse = ", ")))
    }
  }
})

test_that("ADT states are a subset of RNA states (concordance requirement)", {
  rna <- names(read_states(RNA_MARKERS))
  adt <- names(read_states(ADT_MARKERS))
  extra <- setdiff(adt, rna)
  expect_length(extra, 0)
  if (length(extra) > 0) {
    fail(sprintf("ADT-only states break RNA/ADT concordance: %s",
                 paste(extra, collapse = ", ")))
  }
})

test_that("ADT _aliases maps to non-empty character vectors with no self-collisions", {
  d <- jsonlite::fromJSON(ADT_MARKERS, simplifyVector = TRUE)
  aliases <- d[["_aliases"]]
  expect_true(is.list(aliases) && length(aliases) > 0)
  expect_true(all(vapply(aliases, function(x) is.character(x) && length(x) > 0,
                         logical(1))))

  collisions <- intersect(unlist(aliases, use.names = FALSE), names(aliases))
  expect_length(collisions, 0)
})
