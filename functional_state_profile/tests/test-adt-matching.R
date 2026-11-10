skip_adt <- function() {
  if (!exists(".adt_match_features", mode = "function"))
    skip("classify_functional_state.R not sourced")
}

.find_shipped_am <- function(fn) {
  cands <- c(fn, file.path("..", fn),
             file.path(dirname(dirname(FIXTURES)), fn),
             file.path(dirname(FIXTURES), "..", fn))
  for (p in cands) if (file.exists(p)) return(normalizePath(p))
  stop("cannot locate shipped file: ", fn)
}
ADT_MARKERS <- .find_shipped_am("cell_state_markers_adt.json")

expand_with <- function(aliases) {
  function(mk) unique(c(mk, unlist(aliases[mk], use.names = FALSE)))
}

decorate <- function(x, i) {
  style <- i %% 5L
  if (style == 0L) paste0("Ab.", gsub("[^A-Za-z0-9]+", ".", x))
  else if (style == 1L) paste0(x, "-TotalSeqC")
  else if (style == 2L) paste0("Hu-", x)
  else if (style == 3L) paste0(gsub("[^A-Za-z0-9]+", "_", x), "_AbSeq")
  else paste0("anti-human ", x)
}

test_that("every ADT core marker resolves on a decorated canonical-name panel", {
  skip_adt()
  d       <- jsonlite::fromJSON(ADT_MARKERS, simplifyVector = TRUE)
  aliases <- d[["_aliases"]]
  states  <- d[!startsWith(names(d), "_")]
  expand  <- expand_with(aliases)

  all_markers <- unique(unlist(lapply(states, function(s) c(s$core, s$additional)),
                               use.names = FALSE))
  panel <- vapply(seq_along(all_markers),
                  function(i) decorate(all_markers[i], i), character(1))

  for (st in names(states)) {
    for (mk in states[[st]]$core) {
      hits <- .adt_match_features(panel, expand(mk))
      expect_true(length(hits) >= 1,
                  info = sprintf("state '%s': core marker '%s' resolved nothing on the decorated panel",
                                 st, mk))
    }
  }
})

test_that("ADT core markers also resolve when the panel uses ALIAS spellings", {
  skip_adt()
  d       <- jsonlite::fromJSON(ADT_MARKERS, simplifyVector = TRUE)
  aliases <- d[["_aliases"]]
  expand  <- expand_with(aliases)

  for (canonical in c("PD-1", "TIM3", "CTLA-4", "OX40", "4-1BB", "CD127",
                      "CD62L", "CD95", "CD57", "CD16")) {
    alias_names <- aliases[[canonical]]
    expect_true(length(alias_names) >= 1,
                info = sprintf("no alias entry for %s", canonical))
    panel <- paste0(alias_names, "-TotalSeqC")
    hits  <- .adt_match_features(panel, expand(canonical))
    expect_true(length(hits) >= 1,
                info = sprintf("'%s' did not resolve alias spellings: %s",
                               canonical, paste(alias_names, collapse = ", ")))
  }
})

test_that("whole-token matching does not cross-match related names", {
  skip_adt()

  panel <- c("CD34-TotalSeqC", "Ab.CD38", "ABCD3", "PLCD3", "Hu-CD33")
  expect_length(.adt_match_features(panel, "CD3"), 0)

  expect_length(.adt_match_features(c("CD141_AbSeq"), "CD14"), 0)

  expect_identical(.adt_match_features(c("Ab.CD3", panel), "CD3"), "Ab.CD3")
})

test_that("punctuated markers match dotted and glued panel spellings", {
  skip_adt()
  expect_identical(.adt_match_features(c("Ab.HLA.DR"), "HLA-DR"), "Ab.HLA.DR")
  expect_identical(.adt_match_features(c("HLADR-TotalSeqC"), "HLA-DR"),
                   "HLADR-TotalSeqC")
  expect_identical(.adt_match_features(c("PD-1.TotalSeqC"), "PD-1"),
                   "PD-1.TotalSeqC")
})
