# DATA holds the input tables described in data/README.md.

DATA <- "data"

step_output <- function(step) {
  d <- file.path(getwd(), "analysis", step, "output")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

seq_labels <- function(df) {
  nm <- colnames(df)
  if (!any(grepl("^R3[PW]", nm))) {
    lib <- grepl("^R[12][PW][0-9]$", nm)
    nm[lib] <- sub("^R2", "R3", nm[lib])
    nm[lib] <- sub("^R1", "R2", nm[lib])
    nm[nm == "S1P1"] <- "S1P4"
    colnames(df) <- nm
  }
  df
}
