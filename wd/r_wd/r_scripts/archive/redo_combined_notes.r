### atempt 1

#functions

```{r}
# Updated vertical eye legend function to use circle_parameters with proper scaling
f.fig2_legend_eyes_vertical <- function(offset_1, x_space, plot_x, plot_y, text_size, circle_parameters) {
  # Same circle size calculation
  circle_size <- x_space
  
  # Set up vertical positioning
  x1 <- round((offset_1 * plot_x), 1)  # Position for symbols
  x2 <- x1 + x_space                   # Position for text
  
  # Vertical spacing for elements
  y_top <- plot_y * 0.75               # Position for Bud eye
  y_middle <- plot_y * 0.5             # Position for Aletta eye
  y_bottom <- plot_y * 0.25            # Position for Station circle
  
  # Create data objects
  # Bud eye
  new_eye.b <- data.frame(x=x1, y=y_top, image=b$image[1])
  
  # Aletta eye
  new_eye.a <- data.frame(x=x1, y=y_middle) %>% 
    st_as_sf(., coords=c("x", "y"))
  st_crs(new_eye.a) <- st_crs(a)
  
  # Station circle with proper creation
  center_point <- st_point(c(x1, y_bottom))
  center_sf <- st_sf(geometry = st_sfc(center_point))
  st_crs(center_sf) <- st_crs(a)
  new_circle <- st_buffer(center_sf, dist=circle_size)
  
  # Extract styling from circle_parameters
  circle_color1 <- circle_parameters["circle_color1"]
  circle_color2 <- circle_parameters["circle_color2"]
  base <- as.numeric(circle_parameters["base"])
  outline <- as.numeric(circle_parameters["outline"])
  line1 <- circle_parameters["line1"]
  line2 <- circle_parameters["line2"]
  
  # Adjust line thickness for legend if needed
  # Use larger values to make lines more visible in legend
  legend_base <- base * 0.8    # Still visible but not too thick
  legend_outline <- outline * 0.8
  
  # Text labels
  text1 <- "TC Bud Eye"
  text2 <- "TC Aletta Eye"
  text3 <- "200km Diameter\nStation 3.5"
  
  # Create vertical legend with styling from circle_parameters
  p <- ggplot() +
    # Apply white outline circle with adjusted size
    geom_sf(data = new_circle, color = circle_color1, linetype = line1, linewidth = legend_outline, alpha = 0, inherit.aes = FALSE) +
    # Apply orange dashed circle with adjusted size
    geom_sf(data = new_circle, color = circle_color2, linetype = line2, linewidth = legend_base, alpha = 0, inherit.aes = FALSE) +
    
    # Add Aletta eye
    geom_sf(data=new_eye.a, shape=aletta_symbol, col=aletta_color, fill=aletta_color, size=aletta_size) +
    
    # Add Bud eye with controlled size
    geom_image(data=new_eye.b, aes(x, y, image=image), size=0.15) + 
    scale_size_identity() +
    
    # Add text labels
    annotate("text", x=x2, y=y_top, label=text1, hjust=0, size=text_size) +
    annotate("text", x=x2, y=y_middle, label=text2, hjust=0, size=text_size) +
    annotate("text", x=x2, y=y_bottom, label=text3, hjust=0, size=text_size) +
    
    # Ensure text is properly sized and visible
    theme(text = element_text(size=text_size*3)) +   
    theme_void() +
    
    # Use coord_sf with fixed aspect ratio to ensure circle stays round
    coord_sf(xlim=c(0, plot_x), ylim=c(0, plot_y), expand=FALSE, datum=NA)
  
  return(p)
}

# Completely revised function for creating and positioning the eye legend overlay
f.create_eye_legend_overlay <- function(main_width, main_height, circle_parameters) {
  # Create the vertical eye legend first
  legend_vertical <- f.fig2_legend_eyes_vertical(
    offset_1 = 0.1,       # Position within legend viewport
    x_space = 0.2,        # Spacing between symbol and text (as fraction of width)
    plot_x = 1,           # Use unit coordinates (0-1)
    plot_y = 1,           # Use unit coordinates (0-1)
    text_size = 3,        # Base text size
    circle_parameters = circle_parameters
  )
  
  # Convert to grob
  legend_grob <- ggplotGrob(legend_vertical)
  
  # Create a blank canvas with explicit dimensions
  overlay_plot <- ggplot() +
    theme_void() +
    theme(
      plot.background = element_blank(),
      panel.background = element_blank(),
      plot.margin = margin(0, 0, 0, 0)
    ) +
    # Make sure coordinate system is properly defined
    scale_x_continuous(limits = c(0, main_width), expand = c(0,0)) +
    scale_y_continuous(limits = c(0, main_height), expand = c(0,0))
  
  # Position the legend at the bottom right with exact coordinates
  final_overlay <- overlay_plot + 
    annotation_custom(
      grob = legend_grob,
      xmin = main_width * 0.75,
      xmax = main_width,
      ymin = 0, 
      ymax = main_height * 0.2  # Take up bottom 20% of height
    )
  
  return(final_overlay)
}

# Updated function to save the PDFs with better positioning
f.save_separate_pdfs <- function(params) {
  # Helper function for default values
  "%||%" <- function(x, y) {
    if (is.null(x)) y else x
  }
  
  # Extract parameters with defaults
  main_pdf_name <- params$main_pdf_name %||% "figure_two_main.pdf"
  overlay_pdf_name <- params$overlay_pdf_name %||% "figure_two_overlay.pdf"
  test_pdf_name <- params$test_pdf_name %||% "figure_two_vertical_eye_legend.pdf"
  
  main_width <- params$main_width %||% 8.5
  main_height <- params$main_height %||% 11
  test_width <- params$test_width %||% 3
  test_height <- params$test_height %||% 4
  
  rel_widths <- params$rel_widths %||% c(14, 1)
  
  # Eye legend parameters - simplified for better scaling
  eye_plot_x <- params$eye_plot_x %||% 1
  eye_plot_y <- params$eye_plot_y %||% 1
  eye_text_size <- params$eye_text_size %||% 4
  
  # Circle styling parameters
  circle_parameters <- params$circle_parameters
  
  # 1. Save main figure with color legend
  pdf(main_pdf_name, height=main_height, width=main_width)
  main_plot <- plot_grid(p, legend, ncol=2, rel_widths=rel_widths, align="h")
  print(main_plot)
  dev.off()
  
  # 2. Create a standalone test legend for reference
  standalone_legend <- f.fig2_legend_eyes_vertical(
    offset_1 = 0.2,
    x_space = 0.2,
    plot_x = 1,
    plot_y = 1,
    text_size = eye_text_size,
    circle_parameters = circle_parameters
  )
  
  pdf(test_pdf_name, width=test_width, height=test_height)
  print(standalone_legend)
  dev.off()
  
  # 3. Save overlay with eye legend on blank background
  # Create the overlay with the eye legend properly positioned
  overlay_plot <- f.create_eye_legend_overlay(
    main_width = main_width, 
    main_height = main_height,
    circle_parameters = circle_parameters
  )
  
  # Save overlay with exact dimensions matching main figure
  pdf(overlay_pdf_name, height=main_height, width=main_width)
  print(overlay_plot)
  dev.off()
  
  message("PDFs created successfully. Now run the terminal command to combine them.")
}


```


```{r}
# Updated parameter list for figure 2 with simpler positioning parameters
figure_two_params <- list(
  main_pdf_name = "figure_two_main.pdf", 
  overlay_pdf_name = "figure_two_overlay.pdf",
  test_pdf_name = "figure_two_vertical_eye_legend_test.pdf", 
  main_width = 8.5, 
  main_height = 11, 
  test_width = 3, 
  test_height = 4,
  rel_widths = c(14, 1),
  eye_text_size = 4,
  circle_parameters = circle_parameters
)
```

### new overlay apporach:
functiosn
```{r}
# New function to create a vertical eye legend with 90-degree rotated text
f.create_eye_legend_vertical <- function(text_size, circle_parameters, aletta_size_legend = NULL) {
  # Use provided aletta_size_legend or default to global aletta_size
  if (is.null(aletta_size_legend)) {
    aletta_size_legend <- aletta_size * 3  # Default to 3x size for legend
  }
  
  # Extract styling from circle_parameters
  circle_color1 <- as.character(circle_parameters["circle_color1"])
  circle_color2 <- as.character(circle_parameters["circle_color2"])
  base <- as.numeric(circle_parameters["base"])
  outline <- as.numeric(circle_parameters["outline"])
  line1 <- as.character(circle_parameters["line1"])
  line2 <- as.character(circle_parameters["line2"])
  
  # Create data for positioning elements
  # Everything is in normalized coordinates (0-1)
  elements_df <- data.frame(
    y = c(0.8, 0.5, 0.2),           # Vertical positioning (top to bottom)
    x = rep(0.3, 3),                # All elements aligned on left side
    label = c("TC Bud Eye", "TC Aletta Eye", "200km Diameter\nStation 3.5"),
    element = c("bud", "aletta", "station")
  )
  
  # Create Bud eye data
  bud_df <- data.frame(
    x = elements_df$x[1],
    y = elements_df$y[1],
    image = b$image[1]
  )
  
  # Create Aletta eye
  aletta_point <- st_point(c(elements_df$x[2], elements_df$y[2]))
  aletta_sf <- st_sf(geometry = st_sfc(aletta_point))
  st_crs(aletta_sf) <- st_crs(a)
  
  # Create station circle
  station_point <- st_point(c(elements_df$x[3], elements_df$y[3]))
  station_sf <- st_sf(geometry = st_sfc(station_point))
  st_crs(station_sf) <- st_crs(a)
  circle_size <- 0.05  # Size of circle in normalized units
  station_circle <- st_buffer(station_sf, dist = circle_size)
  
  # Create the plot
  p <- ggplot() +
    # Add station circle with white outline
    geom_sf(data = station_circle, color = circle_color1, 
            linetype = line1, linewidth = outline * 0.3, 
            fill = NA, alpha = 1, inherit.aes = FALSE) +
    
    # Add station circle with orange dashed line
    geom_sf(data = station_circle, color = circle_color2, 
            linetype = line2, linewidth = base * 0.3, 
            fill = NA, alpha = 1, inherit.aes = FALSE) +
    
    # Add Aletta eye symbol with specified size
    geom_sf(data = aletta_sf, shape = aletta_symbol, 
            color = aletta_color, fill = aletta_color, 
            size = aletta_size_legend, inherit.aes = FALSE) +
    
    # Add Bud eye hurricane image
    geom_image(data = bud_df, aes(x, y, image = image), size = 0.1) +
    scale_size_identity() +
    
    # Add text labels with 90-degree rotation (vertical text)
    geom_text(data = elements_df, 
              aes(x = x + 0.15, y = y, label = label),
              angle = 270, hjust = 0.5, vjust = 0, 
              size = text_size) +
    
    # Clean theme with transparent background
    theme_void() +
    theme(
      plot.background = element_rect(fill = NA, color = NA),
      panel.background = element_rect(fill = NA, color = NA)
    ) +
    
    # Set coordinate system with fixed aspect ratio
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE)
  
  return(p)
}

# New function to create an overlay that matches main figure layout
f.create_layout_matched_overlay <- function(main_width, main_height, rel_widths, circle_parameters, aletta_size_legend = NULL) {
  # Calculate column widths based on rel_widths
  total_width <- sum(rel_widths)
  left_width <- (rel_widths[1] / total_width) * main_width
  right_width <- (rel_widths[2] / total_width) * main_width
  
  # Create an empty plot for the left column (where the main figure panels would be)
  empty_plot <- ggplot() + 
    theme_void() + 
    theme(
      plot.background = element_rect(fill = NA, color = NA),
      panel.background = element_rect(fill = NA, color = NA)
    ) +
    xlim(0, 1) + ylim(0, 1)
  
  # Create the eye legend for the right column
  eye_legend <- f.create_eye_legend_vertical(
    text_size = 4,
    circle_parameters = circle_parameters,
    aletta_size_legend = aletta_size_legend
  )
  
  # Combine empty left column with eye legend in right column
  # This will match the layout structure of the main figure
  overlay_plot <- plot_grid(
    empty_plot, eye_legend,
    ncol = 2,
    rel_widths = rel_widths,
    align = "h"
  )
  
  # Return the complete layout-matched overlay
  return(overlay_plot)
}

# New function to save PDFs with matched layouts
f.save_matched_layout_pdfs <- function(params) {
  # Helper function for default values
  "%||%" <- function(x, y) {
    if (is.null(x)) y else x
  }
  
  # Extract parameters with defaults
  main_pdf_name <- params$main_pdf_name %||% "figure_two_main.pdf"
  overlay_pdf_name <- params$overlay_pdf_name %||% "figure_two_overlay.pdf"
  test_pdf_name <- params$test_pdf_name %||% "figure_two_eye_legend_test.pdf"
  combined_pdf_name <- params$combined_pdf_name %||% "figure_two_combined.pdf"
  
  main_width <- params$main_width %||% 8.5
  main_height <- params$main_height %||% 11
  test_width <- params$test_width %||% 3
  test_height <- params$test_height %||% 4
  
  rel_widths <- params$rel_widths %||% c(14, 1)
  
  # Circle styling parameters
  circle_parameters <- params$circle_parameters
  
  # Use a larger aletta_size for the legend if specified
  aletta_size_legend <- params$aletta_size_legend %||% (aletta_size * 3)
  
  # 1. Save main figure with color legend
  pdf(main_pdf_name, height = main_height, width = main_width)
  main_plot <- plot_grid(p, legend, ncol = 2, rel_widths = rel_widths, align = "h")
  print(main_plot)
  dev.off()
  
  # 2. Create a standalone test version of the eye legend
  test_legend <- f.create_eye_legend_vertical(
    text_size = 4,
    circle_parameters = circle_parameters,
    aletta_size_legend = aletta_size_legend
  )
  
  pdf(test_pdf_name, width = test_width, height = test_height)
  print(test_legend)
  dev.off()
  
  # 3. Create and save the overlay with matched layout structure
  overlay_plot <- f.create_layout_matched_overlay(
    main_width = main_width,
    main_height = main_height,
    rel_widths = rel_widths,
    circle_parameters = circle_parameters,
    aletta_size_legend = aletta_size_legend
  )
  
  # Save overlay with transparent background
  pdf(overlay_pdf_name, height = main_height, width = main_width)
  par(bg = "transparent")
  print(overlay_plot)
  dev.off()
  
  message("PDFs created successfully with matched layouts.")
  
  # Return file paths for use with combine function
  return(list(
    main = main_pdf_name,
    overlay = overlay_pdf_name,
    test = test_pdf_name,
    combined = combined_pdf_name
  ))
}

# Function to combine PDFs (same as before)
f.combine_pdfs_matched <- function(files = NULL) {
  # Default files if not provided
  if (is.null(files)) {
    files <- list(
      main = "figure_two_main.pdf",
      overlay = "figure_two_overlay.pdf",
      combined = "figure_two_combined.pdf"
    )
  }
  
  # Check if pdftk is available
  has_pdftk <- system("which pdftk", ignore.stdout = TRUE, ignore.stderr = TRUE) == 0
  
  # Check if ghostscript is available
  has_gs <- system("which gs", ignore.stdout = TRUE, ignore.stderr = TRUE) == 0
  
  # Combine PDFs using available tool
  if (has_pdftk) {
    # Use stamp mode for proper overlay
    cmd <- sprintf("pdftk %s stamp %s output %s", 
                   files$main, files$overlay, files$combined)
    result <- system(cmd)
    if (result == 0) {
      message("PDFs combined using pdftk stamp mode")
    } else {
      # Fallback to different pdftk command if stamp fails
      cmd <- sprintf("pdftk %s background %s output %s", 
                     files$overlay, files$main, files$combined)
      result <- system(cmd)
      if (result == 0) {
        message("PDFs combined using pdftk background mode")
      } else {
        message("pdftk commands failed, trying alternative method")
        return(FALSE)
      }
    }
  } else if (has_gs) {
    # Use ghostscript to overlay PDFs
    cmd <- sprintf("gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=%s %s %s", 
                   files$combined, files$main, files$overlay)
    result <- system(cmd)
    if (result == 0) {
      message("PDFs combined using Ghostscript")
    } else {
      message("Ghostscript command failed")
      return(FALSE)
    }
  } else {
    message("Neither pdftk nor Ghostscript found. Install one of these tools for PDF combining")
    return(FALSE)
  }
  
  # Check if successful
  if (file.exists(files$combined)) {
    message("Combined PDF created successfully: ", files$combined)
    return(TRUE)
  } else {
    message("Error: Combined PDF was not created. Check terminal output for errors.")
    return(FALSE)
  }
}

# One-step function to create and combine PDFs with matched layouts
f.create_figure_two_matched <- function(params = NULL) {
  # Use default parameters if none provided
  if (is.null(params)) {
    params <- list(
      main_pdf_name = "figure_two_main.pdf",
      overlay_pdf_name = "figure_two_overlay.pdf",
      test_pdf_name = "figure_two_eye_legend_test.pdf",
      combined_pdf_name = "figure_two_combined.pdf",
      main_width = 8.5,
      main_height = 11,
      test_width = 3,
      test_height = 4,
      rel_widths = c(14, 1),
      circle_parameters = circle_parameters,
      aletta_size_legend = aletta_size * 3
    )
  }
  
  # Create individual PDFs with matched layouts
  files <- f.save_matched_layout_pdfs(params)
  
  # Combine PDFs
  success <- f.combine_pdfs_matched(files)
  
  if (success) {
    message("Figure Two created successfully with matched layout approach!")
  } else {
    message("Warning: Figure Two creation might be incomplete.")
  }
  
  return(success)
}

# Example usage:
f.create_figure_two_matched()
```

input
```{r}
setwd(wd)
setwd(fig)

params <- list(
  # Output filenames
  main_pdf_name = "figure_two_main.pdf",
  overlay_pdf_name = "figure_two_overlay.pdf",
  test_pdf_name = "figure_two_eye_legend_test.pdf",
  combined_pdf_name = "figure_two_combined.pdf",
  
  # PDF dimensions
  main_width = 8.5,
  main_height = 11,
  test_width = 3,
  test_height = 4,
  
  # Layout proportions
  rel_widths = c(14, 1),  # Ratio of main figure to legend column
  
  # Circle styling parameters (should match what you've already defined)
  circle_parameters = circle_parameters,  # Use your existing named vector
  
  # Aletta symbol size for the legend
  aletta_size_legend = 6  # Adjust this value to match your desired triangle size
)

files <- f.save_matched_layout_pdfs(params)
f.combine_pdfs_matched(files)
```
running
```{r}

```


