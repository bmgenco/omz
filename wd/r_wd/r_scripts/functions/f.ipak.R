# f.ipak: load packages, installing any that are missing
# usage: f.ipak(c("sf", "oce", "ggplot2"))
f.ipak <- function(pkg) {
  new <- setdiff(pkg, rownames(installed.packages()))
  if (length(new)) install.packages(new, dependencies = TRUE)
  invisible(lapply(pkg, library, character.only = TRUE,
                   quietly = TRUE, warn.conflicts = FALSE))
}
