.sf_error_context <- new.env(parent = emptyenv())
.sf_error_context$smart_function <- "unknown"
.sf_error_context$step           <- "startup"

#' Register which smart function is running (goes into error.json verbatim)
sf_error_init <- function(smart_function) {
  .sf_error_context$smart_function <- smart_function
  invisible(NULL)
}

#' Record the pipeline stage currently running, for error attribution
sf_step <- function(step) {
  .sf_error_context$step <- step
  invisible(NULL)
}

#' Resolve the smart function's output directory
#'
#' @return The value of `SF_OUTPUT_PATH`, or `"run/out"` when unset (local
#'   runs outside the platform).
.sf_output_dir <- function() {
  Sys.getenv("SF_OUTPUT_PATH", unset = "run/out")
}

#' Write $SF_OUTPUT_PATH/error.json; never throws, and returns the path invisibly.
sf_write_error_json <- function(code, message, details = character(0)) {
  out_dir <- .sf_output_dir()
  tryCatch({
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }
  }, error = function(e) NULL)

  payload <- list(
    status         = "error",
    smart_function = .sf_error_context$smart_function,
    error_code     = as.character(code),
    message        = as.character(message),
    details        = as.character(details),
    step           = .sf_error_context$step,
    timestamp      = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )
  path <- file.path(out_dir, "error.json")
  tryCatch({
    writeLines(jsonlite::toJSON(payload, auto_unbox = TRUE, pretty = TRUE),
               path)
  }, error = function(e) {
    cat("WARNING: could not write error.json: ", conditionMessage(e), "\n",
        sep = "", file = stderr())
  })

  tryCatch({
    port_dir <- file.path(out_dir, "output")
    if (!dir.exists(port_dir)) {
      dir.create(port_dir, recursive = TRUE, showWarnings = FALSE)
    }
    file.copy(path, file.path(port_dir, "error.json"), overwrite = TRUE)
  }, error = function(e) NULL)
  invisible(path)
}

#' Report a structured error, write error.json and the FAILURE status file, and exit.
sf_fail <- function(code, message, details = character(0), status = 1) {
  cat(sprintf("\n[%s] %s\n", code, message))
  if (length(details) > 0) {
    cat(paste0("  - ", details, collapse = "\n"), "\n", sep = "")
  }
  sf_write_error_json(code, message, details)
  status_file <- Sys.getenv("SF_STATUS_FILE", unset = "")
  if (nzchar(status_file)) {
    tryCatch(writeLines("FAILURE", status_file), error = function(e) NULL)
  }
  quit(status = status, save = "no")
}

#' Backstop handler that writes an error.json (PIPELINE_FAILED) for any uncaught error.
sf_uncaught_error_handler <- function() {
  msg <- geterrmessage()
  sf_fail("PIPELINE_FAILED",
          "The analysis stopped due to an unexpected error.",
          details = c(sub("\n$", "", msg),
                      sprintf("step: %s", .sf_error_context$step)))
}
