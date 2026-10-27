# f.raster_kml: classify a lon/lat raster into bands and write styled polygons to KML
# r:       SpatRaster in EPSG:4326
# breaks:  band edges (last can be Inf); bands are [a, b)
# labels:  one label per band
# colours: one hex colour ("#RRGGBB") per band
# layer:   KML layer name
# file:    output .kml path
# alpha:   fill opacity as 2-digit hex ("B3" ~ 70%)
# usage: f.raster_kml(p, c(0.05, 0.1, Inf), c("5–10%", "≥10%"), c("#c7e9c0", "#238b45"), "p_annual", "p.kml")
f.raster_kml <- function(r, breaks, labels, colours, layer, file, alpha = "B3", description = "") {
  n   <- length(breaks) - 1
  rcl <- cbind(breaks[-(n + 1)], breaks[-1], seq_len(n))
  cls <- terra::classify(r, rcl, right = FALSE, include.lowest = TRUE, others = NA)

  pol <- sf::st_as_sf(terra::as.polygons(cls, dissolve = TRUE))
  names(pol)[1] <- "band"
  col <- substr(colours[pol$band], 1, 7)

  pol$Name        <- labels[pol$band]
  pol$Description <- paste0(layer, ": ", labels[pol$band], ifelse(description == "", "", paste0(". ", description)))
  pol$OGR_STYLE   <- paste0("BRUSH(fc:", col, alpha, ");PEN(c:", col, "FF,w:0.5px)")
  pol <- pol[, c("Name", "Description", "OGR_STYLE")]

  sf::st_write(pol, file, layer = layer, driver = "LIBKML", delete_dsn = TRUE, quiet = TRUE)
  invisible(file)
}
