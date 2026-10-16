# f.load_argo: source OneArgo-R and initialize with project directories
# requires paths.R (argo_code, argo_dir, argo_prof, argo_index) and f.ipak
# usage: f.load_argo()            # reuse existing profile files
#        f.load_argo(update = 1)  # re-download updated profiles
f.load_argo <- function(update = 0) {
  f.ipak(c("gsw", "R.utils", "Matrix"))

  # source OneArgo functions into the global env (they use <<-)
  main <- list.files(argo_code, pattern = "\\.R$", full.names = TRUE)
  main <- main[!basename(main) %in% c("Tutorial.R", "oneargo_r_license.R")]
  aux  <- list.files(file.path(argo_code, "auxil"), pattern = "\\.R$", full.names = TRUE)
  for (f in c(main, aux)) source(f, local = .GlobalEnv)

  dir.create(argo_prof,  recursive = TRUE, showWarnings = FALSE)
  dir.create(argo_index, recursive = TRUE, showWarnings = FALSE)

  # initialize_argo() writes Index/ relative to the working directory
  old <- setwd(argo_dir); on.exit(setwd(old), add = TRUE)
  initialize_argo()

  # must follow init, which rebuilds Setting from scratch
  Setting$prof_dir  <<- paste0(argo_prof, "/")
  Setting$hosts     <<- "https://data-argo.ifremer.fr/"   # older submodule commits list a dead GODAE host first
  Setting$use_pause <<- 0
  Setting$update    <<- update                          # 0 = reuse existing profile files
  invisible(TRUE)
}
