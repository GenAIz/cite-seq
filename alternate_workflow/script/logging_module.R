`%||%` <- function(a, b) if (!is.null(a)) a else b

.logger <- list(
  log_file = NULL,
  start_time = NULL,
  initialized = FALSE,
  con = NULL,
  tee = FALSE
)

initialize_logger <- function(config_file) {
  tryCatch({
    if (is.null(config_file) || !file.exists(config_file)) {
      .logger$log_file <<- "tcell_workflow_default.log"
    } else {
      config <- jsonlite::fromJSON(config_file)

      env_logs_dir <- Sys.getenv("SF_LOG_DIR", unset = "")
      logs_dir <- if (nzchar(env_logs_dir)) env_logs_dir else (config$paths$logs_dir %||% "logs")

      if (!dir.exists(logs_dir)) {
        dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
      }

      proj <- config$metadata$project_name %||% config$load_data$project_name %||% "workflow"
      proj <- gsub("[^A-Za-z0-9._-]+", "_", as.character(proj)[1])
      proj <- gsub("^[._-]+|[._-]+$", "", proj)
      if (!nzchar(proj)) proj <- "workflow"

      .logger$log_file <<- file.path(logs_dir,
        sprintf("tcell_workflow_%s_%s.log", proj, format(Sys.time(), "%Y%m%d_%H%M%S")))
    }
  }, error = function(e) {
    .logger$log_file <<- "tcell_workflow_default.log"
  })
  
  .logger$start_time <<- Sys.time()
  .logger$initialized <<- TRUE
  
  header <- sprintf("===================================================\n")
  header <- paste0(header, "T cell Workflow Log\n")
  header <- paste0(header, sprintf("Started: %s\n", format(.logger$start_time, "%Y-%m-%d %H:%M:%S")))
  header <- paste0(header, "===================================================\n\n")
  
  write(header, file = .logger$log_file, append = FALSE)
}

start_log_tee <- function() {
  if (!isTRUE(.logger$initialized) || isTRUE(.logger$tee)) return(invisible(FALSE))
  ok <- tryCatch({
    con <- file(.logger$log_file, open = "at")
    sink(con, split = TRUE)
    sink(con, type = "message")
    .logger$con <<- con
    .logger$tee <<- TRUE
    TRUE
  }, error = function(e) {
    .logger$con <<- NULL; .logger$tee <<- FALSE
    cat(sprintf("Log tee setup failed (%s); continuing with milestone-only logging.\n",
                conditionMessage(e)))
    FALSE
  })
  invisible(ok)
}

stop_log_tee <- function() {
  if (!isTRUE(.logger$tee)) return(invisible(FALSE))
  try(sink(type = "message"), silent = TRUE)
  try(sink(), silent = TRUE)
  try(if (!is.null(.logger$con)) close(.logger$con), silent = TRUE)
  .logger$tee <<- FALSE
  .logger$con <<- NULL
  invisible(TRUE)
}

log_milestone <- function(milestone, details = "") {
  if (!.logger$initialized) return(invisible(NULL))
  
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  elapsed <- difftime(Sys.time(), .logger$start_time, units = "secs")
  
  message <- sprintf("[%s] [OK] %s", timestamp, milestone)
  if (details != "") {
    message <- paste0(message, " | ", details)
  }
  message <- paste0(message, sprintf(" (+%.1fs)\n", as.numeric(elapsed)))
  
  cat(message)
  if (!isTRUE(.logger$tee)) write(message, file = .logger$log_file, append = TRUE)
}

log_warning <- function(category, message) {
  if (!.logger$initialized) return(invisible(NULL))
  
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_entry <- sprintf("[%s] [WARN] %s - %s\n", timestamp, category, message)
  
  cat(log_entry)
  if (!isTRUE(.logger$tee)) write(log_entry, file = .logger$log_file, append = TRUE)
}

log_error <- function(category, message) {
  if (!.logger$initialized) return(invisible(NULL))
  
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_entry <- sprintf("[%s] [ERROR] ERROR: %s - %s\n", timestamp, category, message)
  
  cat(log_entry)
  if (!isTRUE(.logger$tee)) write(log_entry, file = .logger$log_file, append = TRUE)
}

log_data_summary <- function(obj, description = "") {
  if (!.logger$initialized) return(invisible(NULL))
  
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  summary_text <- sprintf("\n[%s] %s\n", timestamp, description)
  
  if (inherits(obj, "Seurat")) {
    summary_text <- paste0(summary_text, 
      sprintf(" - Cells: %d\n", ncol(obj)),
      sprintf(" - Features: %d\n", nrow(obj)),
      sprintf(" - Metadata columns: %d\n", ncol(obj@meta.data)))
  }
  
  summary_text <- paste0(summary_text, "\n")
  
  cat(summary_text)
  if (!isTRUE(.logger$tee)) write(summary_text, file = .logger$log_file, append = TRUE)
}

log_config_snapshot <- function(config_file) {
  if (!.logger$initialized) return(invisible(NULL))

  section <- paste0(
    "\n---------------------------------------------------\n",
    "CONFIG SNAPSHOT\n",
    sprintf("File: %s\n", config_file),
    sprintf("Captured: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "---------------------------------------------------\n"
  )

  json_text <- tryCatch(
    paste(readLines(config_file, warn = FALSE), collapse = "\n"),
    error = function(e) sprintf("(could not read config file: %s)", e$message)
  )

  section <- paste0(section, json_text, "\n",
    "---------------------------------------------------\n\n")

  if (isTRUE(.logger$tee) && !is.null(.logger$con)) {
    cat(section, file = .logger$con); flush(.logger$con)
  } else {
    write(section, file = .logger$log_file, append = TRUE)
  }
}

finalize_log <- function() {
  if (!.logger$initialized) return(invisible(NULL))
  
  elapsed <- difftime(Sys.time(), .logger$start_time, units = "secs")
  
  footer <- sprintf("\n===================================================\n")
  footer <- paste0(footer, sprintf("Completed: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  footer <- paste0(footer, sprintf("Total time: %.1f seconds (%.1f minutes)\n", 
                                   as.numeric(elapsed), as.numeric(elapsed)/60))
  footer <- paste0(footer, "===================================================\n")
  
  cat(footer)
  if (!isTRUE(.logger$tee)) write(footer, file = .logger$log_file, append = TRUE)
  stop_log_tee()
}
