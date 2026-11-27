# Shared paths for the R scripts. Run every script from the repository root.
#
# DATA holds the input tables described in data/README.md.
# step_output("06_FEAST") returns analysis/06_FEAST/output as an absolute path and
# creates it. The path is absolute because FEAST changes the working directory
# while it runs, which breaks relative output paths.

DATA <- "data"

step_output <- function(step) {
  d <- file.path(getwd(), "analysis", step, "output")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

# The scripts use the sequencing labels (R2 = manuscript R1, R3 = manuscript R2,
# S1P4 = manuscript S1P1). A table with manuscript labels is converted on read.
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
