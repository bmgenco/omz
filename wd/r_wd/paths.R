# wd/r_wd/paths.R : single source of truth for directories
# Source at the top of every script: source(here::here("wd", "r_wd", "paths.R"))
# All paths are absolute, resolved from the project root (omz/).
library(here)

# --- code ---
r_wd        <- here("wd", "r_wd")
scripts_d   <- here("wd", "r_wd", "r_scripts")
fun_d       <- here("wd", "r_wd", "r_scripts", "functions")

# --- project data & outputs ---
robj        <- here("wd", "r_wd", "r_objects")
fig         <- here("figures")
data_d      <- here("data")
gis_data    <- here("data", "gis_data")
ocean_color <- here("data", "ocean_color_bud")
bathy_d     <- here("data", "bathy")

# --- OneArgo-R ---
argo_code   <- here("wd", "r_wd", "r_scripts", "submodules", "OneArgo-R")
argo_dir    <- here("wd", "r_wd")               # initialize_argo() runs here -> Index/ lands in wd/r_wd/Index
argo_prof   <- here("wd", "r_wd", "Profiles")   # existing downloads (gitignored)
argo_index  <- here("wd", "r_wd", "Index")

# --- machine-specific paths (e.g. Dropbox) in an untracked file ---
local_paths <- here("wd", "r_wd", "local_paths.R")
if (file.exists(local_paths)) source(local_paths)

# --- core utility functions ---
source(file.path(fun_d, "f.ipak.R"))
source(file.path(fun_d, "f.load_argo.R"))

# --- create output dirs if missing ---
for (d in c(robj, fig)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
