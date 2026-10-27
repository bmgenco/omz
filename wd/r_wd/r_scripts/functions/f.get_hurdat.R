
# f.get_hurdat: download NHC HURDAT2, parse, delete the download
# Replaces HURDAT::get_hurdat() (package archived 2021; its URLs are dead).
#
# basin: "EP" = NE & North-Central Pacific (includes CP storms), "AL" = Atlantic
# url:   optional; overrides the default release for that basin
#        (newer releases: https://www.nhc.noaa.gov/data/hurdat/)
#
# Returns all 21 HURDAT2 fields plus Key and Name:
#   Key, Name, DateTime, Record, Status, Lat, Lon, Wind (kt), Pressure (mb),
#   NE34..NW34, NE50..NW50, NE64..NW64 (wind radii, n mi), RMW (radius of max wind, n mi)
# Missing values (-999) are NA. Wind radii begin in 2004; RMW in 2021.
# Lon: west negative, east positive (points past the dateline stay positive)
#
# usage: h <- f.get_hurdat()                 # EP, default release
#        h <- f.get_hurdat(basin = "AL")
#        h <- f.get_hurdat(basin = "EP", url = "<newer release URL>")
f.get_hurdat <- function(basin = c("EP", "AL"), url = NULL) {
  
  basin <- match.arg(basin)
  
  # default releases (NHC, posted 2026-09-30)
  urls <- c(
    EP = "https://www.nhc.noaa.gov/data/hurdat/hurdat2-nepac-1949-2025-092926.txt",
    AL = "https://www.nhc.noaa.gov/data/hurdat/hurdat2-1851-2025-092326.txt"
  )
  if (is.null(url)) url <- urls[[basin]]
  
  # download to a temp file; delete it however the function exits
  tmp <- tempfile(fileext = ".txt")
  on.exit(unlink(tmp), add = TRUE)
  download.file(url, tmp, quiet = TRUE, mode = "wb")
  
  x   <- readr::read_lines(tmp)
  hdr <- grepl("^[A-Z]{2}[0-9]{6},", x)            # storm header lines, e.g. "EP032018, BUD, 30,"
  
  hp   <- stringr::str_split_fixed(x[hdr], ",", 4)
  keys <- trimws(hp[, 1])
  nms  <- trimws(hp[, 2])
  id   <- cumsum(hdr)[!hdr]                         # storm index for each data line
  
  d <- stringr::str_split_fixed(x[!hdr], ",", 21)   # all 21 fields
  
  lat <- trimws(d[, 5]); lon <- trimws(d[, 6])
  num <- function(v) { v <- suppressWarnings(as.numeric(trimws(v))); v[v == -999] <- NA; v }
  
  h <- data.frame(
    Key      = keys[id],
    Name     = nms[id],
    DateTime = lubridate::ymd_hm(paste(trimws(d[, 1]), trimws(d[, 2])), tz = "UTC"),
    Record   = trimws(d[, 3]),
    Status   = trimws(d[, 4]),
    Lat      = as.numeric(sub("[NS]$", "", lat)) * ifelse(grepl("S$", lat), -1, 1),
    Lon      = as.numeric(sub("[EW]$", "", lon)) * ifelse(grepl("W$", lon), -1, 1),
    Wind     = num(d[, 7]),
    Pressure = num(d[, 8]),
    stringsAsFactors = FALSE
  )
  
  # wind radii: 34, 50, 64 kt by quadrant (fields 9-20), then RMW (field 21)
  radii <- paste0(rep(c("NE", "SE", "SW", "NW"), 3), rep(c(34, 50, 64), each = 4))
  for (i in seq_along(radii)) h[[radii[i]]] <- num(d[, 8 + i])
  h$RMW <- num(d[, 21])
  
  attr(h, "basin")  <- basin
  attr(h, "source") <- url     # record which release this came from
  h
}