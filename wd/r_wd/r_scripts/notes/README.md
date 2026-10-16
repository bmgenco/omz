---
title: "previous_code"
output: html_document
date: '2025-03-14'
---
This first code block defines all the functions needed for creating duplicate structure PDFs:

```{r define_functions}
# Function to create an overlay with duplicate structure that preserves axis text
f.create_duplicate_structure_overlay <- function(
  plot_obj, 
  legend_obj, 
  rel_widths = c(14, 1),
  legend_xmin = NULL,  # Will be calculated if NULL
  legend_xmax = 1,
  legend_ymin = 0,
  legend_ymax = 0.25,
  preserve_axis_text = TRUE,  # New parameter to control axis text preservation
  create_diagnostic = FALSE
) {
  # Helper function for null values
  "%||%" <- function(x, y) {
    if (is.null(x)) y else x
  }
  
  # Calculate legend_xmin if not provided
  legend_xmin <- legend_xmin %||% (1 - (rel_widths[2] / sum(rel_widths)))
  
  # Create a version of the main/supplement plot
  # Selectively make elements transparent, potentially preserving axis text
  invisible_plot <- plot_obj + 
    theme(
      panel.background = element_rect(fill = NA, color = NA),
      plot.background = element_rect(fill = NA, color = NA),
      panel.grid = element_blank(),
      # Conditionally preserve axis text
      axis.text = if (preserve_axis_text) element_text() else element_text(color = NA),
      axis.title = element_text(color = NA),
      axis.ticks = if (preserve_axis_text) element_line() else element_line(color = NA),
      legend.position = "none",
      plot.margin = margin(0, 0, 0, 0)
    )
  
  # Create an empty placeholder for the top of the right column
  # This is where the original color legend would be
  empty_legend_top <- ggplot() + 
    theme_void() +
    theme(
      panel.background = element_rect(fill = NA, color = NA),
      plot.background = element_rect(fill = NA, color = NA)
    )
  
  # Create the main structure with invisible plot and empty top legend
  overlay_structure <- plot_grid(
    invisible_plot, 
    empty_legend_top,
    ncol = 2, 
    rel_widths = rel_widths, 
    align = "h"
  )
  
  # If diagnostic mode is requested, create a grid overlay instead
  if (create_diagnostic) {
    # Create a grid overlay to visualize coordinate space
    grid_plot <- ggplot() +
      # Add vertical lines at 0.1 intervals
      geom_vline(xintercept = seq(0, 1, by = 0.1), color = "blue", alpha = 0.5) +
      # Add horizontal lines at 0.1 intervals
      geom_hline(yintercept = seq(0, 1, by = 0.1), color = "red", alpha = 0.5) +
      # Add text labels for coordinates
      annotate("text", x = seq(0, 1, by = 0.1), y = 0.02, 
               label = paste0(seq(0, 1, by = 0.1)), 
               color = "blue", angle = 90, hjust = 0) +
      annotate("text", x = 0.02, y = seq(0, 1, by = 0.1), 
               label = paste0(seq(0, 1, by = 0.1)), 
               color = "red", hjust = 0) +
      # Highlight the legend position rectangle
      annotate("rect", 
               xmin = legend_xmin, xmax = legend_xmax,
               ymin = legend_ymin, ymax = legend_ymax,
               alpha = 0.2, fill = "green") +
      # Add a note about the highlighted area
      annotate("text", x = (legend_xmin + legend_xmax)/2, y = (legend_ymin + legend_ymax)/2,
               label = "Legend\nPosition", color = "darkgreen") +
      # Clean theme with labels
      theme_minimal() +
      labs(title = "Coordinate Grid Overlay",
           subtitle = paste0("Legend position: xmin=", legend_xmin, 
                            ", xmax=", legend_xmax,
                            ", ymin=", legend_ymin,
                            ", ymax=", legend_ymax))
    
    # Create a version of the overlay with grid lines for diagnostic purposes
    diagnostic_overlay <- overlay_structure + 
      annotation_custom(
        grob = ggplotGrob(grid_plot),
        xmin = 0, xmax = 1, ymin = 0, ymax = 1
      )
    
    return(diagnostic_overlay)
  }
  
  # Add the legend at the specified position
  final_overlay <- overlay_structure + 
    annotation_custom(
      grob = legend_obj,  # Use provided legend object
      xmin = legend_xmin,
      xmax = legend_xmax,
      ymin = legend_ymin,
      ymax = legend_ymax
    )
  
  return(final_overlay)
}

# Function to save PDFs with duplicate structure
f.save_duplicate_structure_pdfs <- function(
  plot_obj,
  legend_obj,
  plot_type = "main",
  main_pdf_name = NULL,
  overlay_pdf_name = NULL,
  combined_pdf_name = NULL,
  main_width = 8.5,
  main_height = 11,
  rel_widths = c(14, 1),
  legend_xmin = NULL,
  legend_xmax = 1,
  legend_ymin = 0,
  legend_ymax = 0.25,
  preserve_axis_text = TRUE,  # Added parameter for axis text preservation
  create_diagnostic = FALSE
) {
  # Helper function for null values
  "%||%" <- function(x, y) {
    if (is.null(x)) y else x
  }
  
  # Set default filenames based on plot_type
  main_pdf_name <- main_pdf_name %||% paste0("figure_", plot_type, "_main.pdf")
  overlay_pdf_name <- overlay_pdf_name %||% paste0("figure_", plot_type, "_overlay.pdf")
  combined_pdf_name <- combined_pdf_name %||% paste0("figure_", plot_type, "_combined.pdf")
  
  # If diagnostic mode, modify overlay filename
  if (create_diagnostic) {
    overlay_pdf_name <- sub("\\.pdf$", "_diagnostic.pdf", overlay_pdf_name)
  }
  
  # 1. Save main figure with color legend
  pdf(main_pdf_name, height = main_height, width = main_width)
  main_plot <- plot_grid(plot_obj, legend_obj, ncol = 2, rel_widths = rel_widths, align = "h")
  print(main_plot)
  dev.off()
  
  # 2. Create and save overlay with duplicate structure
  overlay_plot <- f.create_duplicate_structure_overlay(
    plot_obj = plot_obj,
    legend_obj = legend_obj,
    rel_widths = rel_widths,
    legend_xmin = legend_xmin,
    legend_xmax = legend_xmax,
    legend_ymin = legend_ymin,
    legend_ymax = legend_ymax,
    preserve_axis_text = preserve_axis_text,  # Pass axis text preservation parameter
    create_diagnostic = create_diagnostic
  )
  
  pdf(overlay_pdf_name, height = main_height, width = main_width)
  par(bg = "transparent")  # Set transparent background
  print(overlay_plot)
  dev.off()
  
  message(paste0(plot_type, " PDFs created successfully."))
  
  # Return file paths
  result <- list(
    main = main_pdf_name,
    overlay = overlay_pdf_name,
    combined = combined_pdf_name
  )
  
  # Add the diagnostic file path if created
  if (create_diagnostic) {
    result$diagnostic <- overlay_pdf_name
  }
  
  return(result)
}

# Function to combine PDFs with error handling and multiple methods
f.combine_duplicate_pdfs <- function(
  files,
  pdftk_options = NULL,
  gs_options = NULL
) {
  # Helper function for null values
  "%||%" <- function(x, y) {
    if (is.null(x)) y else x
  }
  
  # Default options for pdftk and gs if not provided
  pdftk_options <- pdftk_options %||% ""  # No extra options by default
  gs_options <- gs_options %||% ""  # No extra options by default
  
  # Check if pdftk is available
  has_pdftk <- system("which pdftk", ignore.stdout = TRUE, ignore.stderr = TRUE) == 0
  
  # Check if ghostscript is available
  has_gs <- system("which gs", ignore.stdout = TRUE, ignore.stderr = TRUE) == 0
  
  success <- FALSE
  
  # Combine PDFs using available tool
  if (has_pdftk) {
    # Try stamp mode with optional parameters
    cmd <- sprintf("pdftk %s stamp %s %s output %s", 
                   files$main, files$overlay, pdftk_options, files$combined)
    result <- system(cmd)
    
    if (result == 0) {
      message("PDFs combined using pdftk stamp mode")
      success <- TRUE
    } else {
      # Try alternative modes if stamp fails
      alt_modes <- list(
        list(name = "background", cmd = sprintf("pdftk %s background %s %s output %s", 
                                               files$main, files$overlay, pdftk_options, files$combined)),
        list(name = "multistamp", cmd = sprintf("pdftk %s multistamp %s %s output %s", 
                                               files$main, files$overlay, pdftk_options, files$combined))
      )
      
      for (mode in alt_modes) {
        result <- system(mode$cmd)
        if (result == 0) {
          message(paste0("PDFs combined using pdftk ", mode$name, " mode"))
          success <- TRUE
          break
        }
      }
      
      if (!success) {
        message("All pdftk modes failed")
      }
    }
  } else if (has_gs) {
    # Use ghostscript as fallback with optional parameters
    cmd <- sprintf("gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite %s -sOutputFile=%s %s %s", 
                   gs_options, files$combined, files$main, files$overlay)
    result <- system(cmd)
    
    if (result == 0) {
      message("PDFs combined using Ghostscript")
      success <- TRUE
    } else {
      message("Ghostscript command failed")
    }
  } else {
    message("Neither pdftk nor Ghostscript found. Please install one of these tools.")
  }
  
  # Verify the combined file exists
  if (success && file.exists(files$combined)) {
    message("Combined PDF created successfully: ", files$combined)
    return(TRUE)
  } else {
    message("Error: Combined PDF was not created successfully.")
    return(FALSE)
  }
}

# Function to optimize/compress PDF after combining with high-quality settings
f.optimize_pdf <- function(
  input_pdf,
  output_pdf = NULL,
  quality = "prepress",   # Options: screen, ebook, printer, prepress
  compatibility = "1.5",  # PDF compatibility level
  embed_fonts = TRUE,
  compress_images = TRUE
) {
  # Helper function for null values
  "%||%" <- function(x, y) {
    if (is.null(x)) y else x
  }
  
  # Use input filename with _optimized suffix if no output name provided
  if (is.null(output_pdf)) {
    output_pdf <- sub("\\.pdf$", "_optimized.pdf", input_pdf)
  }
  
  # Check if ghostscript is available
  has_gs <- system("which gs", ignore.stdout = TRUE, ignore.stderr = TRUE) == 0
  
  if (!has_gs) {
    message("Ghostscript not found. Cannot optimize PDF.")
    return(FALSE)
  }
  
  # Build options for Ghostscript
  gs_options <- c(
    sprintf("-dPDFSETTINGS=/%s", quality),
    sprintf("-dCompatibilityLevel=%s", compatibility),
    "-dPrinted=false"  # Optimize for screen viewing
  )
  
  # Add font embedding if requested
  if (embed_fonts) {
    gs_options <- c(gs_options,
                    "-dEmbedAllFonts=true",
                    "-dSubsetFonts=true",
                    "-dMaxSubsetPct=100")
  }
  
  # Add image compression if requested
  if (compress_images) {
    gs_options <- c(gs_options,
                    "-dColorImageDownsampleType=/Bicubic",
                    "-dGrayImageDownsampleType=/Bicubic",
                    "-dMonoImageDownsampleType=/Bicubic")
  }
  
  # Combine options into a single string
  gs_options_str <- paste(gs_options, collapse = " ")
  
  # Ghostscript command for optimization
  cmd <- sprintf("gs -sDEVICE=pdfwrite -dNOPAUSE -dQUIET -dBATCH %s -sOutputFile=%s %s",
                gs_options_str, output_pdf, input_pdf)
  
  result <- system(cmd)
  
  if (result == 0) {
    message("PDF optimized successfully: ", output_pdf)
    
    # Get file sizes for reporting
    original_size <- file.size(input_pdf) / (1024 * 1024)  # Convert to MB
    optimized_size <- file.size(output_pdf) / (1024 * 1024)
    reduction <- (1 - optimized_size/original_size) * 100
    
    message(sprintf("File size reduced from %.2f MB to %.2f MB (%.1f%% reduction)",
                   original_size, optimized_size, reduction))
    
    return(TRUE)
  } else {
    message("PDF optimization failed")
    return(FALSE)
  }
}

# Main function to run the duplicate structure test
# With fully parameterized options
f.run_duplicate_structure_test <- function(
  # Required parameters
  plot_obj,
  legend_obj,
  
  # Identification parameters
  plot_type = "two",
  
  # File parameters
  main_pdf_name = NULL,
  overlay_pdf_name = NULL,
  combined_pdf_name = NULL,
  optimized_pdf_name = NULL,
  
  # Dimension parameters
  main_width = 8.5,
  main_height = 11,
  rel_widths = c(14, 1),
  
  # Legend positioning parameters
  legend_xmin = NULL,
  legend_xmax = 1,
  legend_ymin = 0,
  legend_ymax = 0.25,
  
  # Display parameters
  preserve_axis_text = TRUE,  # Added parameter for axis text preservation
  
  # PDF combination parameters
  pdftk_options = NULL,
  gs_options = NULL,
  
  # Optimization parameters
  optimize_pdf = TRUE,
  quality = "prepress",    # Options: screen, ebook, printer, prepress
  compatibility = "1.5",   # PDF compatibility level
  embed_fonts = TRUE,
  compress_images = TRUE,
  
  # Diagnostic parameters
  create_diagnostic = FALSE,
  
  # Cleanup parameters
  delete_intermediate = FALSE
) {
  # Helper function for null values
  "%||%" <- function(x, y) {
    if (is.null(x)) y else x
  }
  
  # Set default filenames based on plot_type
  combined_pdf_name <- combined_pdf_name %||% paste0("figure_", plot_type, "_combined.pdf")
  optimized_pdf_name <- optimized_pdf_name %||% paste0("figure_", plot_type, "_optimized.pdf")
  
  # Create PDFs with duplicate structure
  files <- f.save_duplicate_structure_pdfs(
    plot_obj = plot_obj,
    legend_obj = legend_obj,
    plot_type = plot_type,
    main_pdf_name = main_pdf_name,
    overlay_pdf_name = overlay_pdf_name,
    combined_pdf_name = combined_pdf_name,
    main_width = main_width,
    main_height = main_height,
    rel_widths = rel_widths,
    legend_xmin = legend_xmin,
    legend_xmax = legend_xmax,
    legend_ymin = legend_ymin,
    legend_ymax = legend_ymax,
    preserve_axis_text = preserve_axis_text,  # Pass axis text preservation parameter
    create_diagnostic = create_diagnostic
  )
  
  # Skip PDF combination if only creating diagnostic
  if (create_diagnostic) {
    message("Diagnostic overlay created. Skipping PDF combination.")
    message("To see coordinate grid, check: ", files$diagnostic)
    return(TRUE)
  }
  
  # Combine PDFs
  success <- f.combine_duplicate_pdfs(
    files = files,
    pdftk_options = pdftk_options,
    gs_options = gs_options
  )
  
  # Optimize the combined PDF if requested
  if (success && optimize_pdf) {
    f.optimize_pdf(
      input_pdf = files$combined, 
      output_pdf = optimized_pdf_name,
      quality = quality,
      compatibility = compatibility,
      embed_fonts = embed_fonts,
      compress_images = compress_images
    )
    # Update the combined PDF reference to the optimized version
    files$final <- optimized_pdf_name
  } else {
    files$final <- files$combined
  }
  
  # Delete intermediate files if requested
  if (success && delete_intermediate) {
    file.remove(files$main)
    file.remove(files$overlay)
    if (optimize_pdf) {
      file.remove(files$combined)  # Remove the unoptimized combined file
    }
    message("Intermediate files deleted.")
  }
  
  if (success) {
    message(paste0(plot_type, " figure completed successfully! Check: ", files$final))
  } else {
    message(paste0(plot_type, " figure encountered issues."))
  }
  
  return(success)
}
```



Integarte eyes legnd with maste plot and other legend
```{r}
# Set the working directory to the figures folder
setwd(wd)
setwd(fig)

# Create the horizontal eye legend
legend_horizontal <- f.fig2_legend_eyes_horizontal(
  plot_width = 800,
  plot_height = 150,  # Made shorter since it's a horizontal layout
  text_size = 12,     # Adjusted for better readability
  circle_parameters = circle_parameters
)

# Generate the final PDF with both legends
f.run_two_legend_test(
  plot_obj = p,
  color_legend_obj = legend,
  eye_legend_obj = legend_horizontal,
  plot_type = "two",
  main_pdf_name = "figure_two_main.pdf",
  overlay_pdf_name = "figure_two_eye_overlay.pdf",
  combined_pdf_name = "figure_two_both_legends.pdf",
  optimized_pdf_name = "figure_two_both_legends_optimized.pdf",
  main_width = 8.5,
  main_height = 11,
  rel_widths = c(14, 1),
  eye_legend_xmin = 0.75,   # Right side of the figure 
  eye_legend_xmax = 0.99,   # Not quite to the edge
  eye_legend_ymin = -0.1,   # Position slightly below the plot
  eye_legend_ymax = 0.15,   # About 1/4 of the way up
  preserve_axis_text = TRUE,
  optimize_pdf = TRUE,
  quality = "prepress"
)
```

