# Shared paths for the R figure scripts.
#
# DATA holds the input tables described in the top-level README.md.
# OUT is where figures and derived files are written.

DATA <- "data"
OUT  <- "output"

dir.create(OUT, showWarnings = FALSE)

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
