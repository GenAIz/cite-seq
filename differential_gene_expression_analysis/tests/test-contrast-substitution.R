skip_sub <- function() {
  if (!exists(".substitute_group_names", mode = "function"))
    skip("dge.R not sourced")
}

groups  <- c("CD4+", "CD4+ Tcm", "Severe_DISEASE", "Healthy",
             "Tcm (CD27+)", "Non-classical")
valid   <- make.names(groups)
revmap  <- setNames(valid, groups)

test_that("a trailing + does not leak into the rewritten formula", {
  skip_sub()
  out <- .substitute_group_names("CD4+ - Healthy", groups, revmap)
  expect_equal(out, paste(revmap[["CD4+"]], "-", revmap[["Healthy"]]))

  expect_false(grepl("\\.\\+", out))
})

test_that("labels containing parentheses are substituted, not skipped", {
  skip_sub()
  out <- .substitute_group_names("Tcm (CD27+) - Non-classical", groups, revmap)
  expect_equal(out, paste(revmap[["Tcm (CD27+)"]], "-", revmap[["Non-classical"]]))

  expect_false(grepl("(", out, fixed = TRUE))
  expect_false(grepl(")", out, fixed = TRUE))
  expect_false(grepl("+", out, fixed = TRUE))

  operands <- trimws(strsplit(out, "-", fixed = TRUE)[[1]])
  expect_equal(operands, make.names(operands))
})

test_that("underscored labels survive intact", {
  skip_sub()
  out <- .substitute_group_names("Severe_DISEASE - Healthy", groups, revmap)
  expect_equal(out, paste(revmap[["Severe_DISEASE"]], "-", revmap[["Healthy"]]))
})

test_that("the longest matching label wins", {
  skip_sub()

  out <- .substitute_group_names("`CD4+ Tcm` - `CD4+`", groups, revmap)
  expect_equal(out, paste(revmap[["CD4+ Tcm"]], "-", revmap[["CD4+"]]))
})

test_that("arithmetic in the formula is preserved", {
  skip_sub()
  out <- .substitute_group_names("(CD4+ + Healthy)/2 - Severe_DISEASE", groups, revmap)
  expect_equal(out, sprintf("(%s + %s)/2 - %s",
                            revmap[["CD4+"]], revmap[["Healthy"]], revmap[["Severe_DISEASE"]]))
})

test_that("substitution is single-pass and idempotent in effect", {
  skip_sub()

  g <- c("A", "A.B")
  rv <- setNames(c("A.B", "A.B.C"), g)
  out <- .substitute_group_names("A.B - A", g, rv)
  expect_equal(out, "A.B.C - A.B")
})

test_that("an unknown token is left alone rather than mangled", {
  skip_sub()
  out <- .substitute_group_names("Healthy - Unmapped", groups, revmap)
  expect_equal(out, paste(revmap[["Healthy"]], "- Unmapped"))
})

test_that("Yoshida DISEASE_status levels do not corrupt each other", {
  skip_sub()
  lv  <- c("DISEASE+", "Healthy", "Post-DISEASE")
  rv  <- setNames(make.names(lv), lv)
  coefs <- unname(rv)

  for (fm in c("DISEASE+ - Healthy", "Post-DISEASE - Healthy", "DISEASE+ - Post-DISEASE",
               "(DISEASE+ + Post-DISEASE)/2 - Healthy")) {
    out <- .substitute_group_names(fm, lv, rv)

    toks <- regmatches(out, gregexpr("[A-Za-z0-9._]+", out))[[1]]
    toks <- toks[!grepl("^[0-9.]+$", toks)]
    expect_true(all(toks %in% coefs),
                info = sprintf("formula '%s' produced unknown coefficient(s): %s",
                               fm, paste(setdiff(toks, coefs), collapse = ", ")))
  }
})

test_that("C24mO is not corrupted by the 24mO level", {
  skip_sub()
  lv <- c("12mR", "24mO", "24mR", "3mR", "C24mO")
  rv <- setNames(make.names(lv), lv)
  out <- .substitute_group_names("C24mO - 24mO", lv, rv)
  expect_equal(out, paste(rv[["C24mO"]], "-", rv[["24mO"]]))

  expect_true(startsWith(out, "C24mO"))
})
