## Refined eye legend function that uses the original hurricane image
## Updated eye legend function that uses linewidth instead of size for lines
f.refined_eye_legend <- function(plot_width=250, plot_height=400) {
  # Create a data frame for Aletta eye (triangle)
  aletta_df <- data.frame(x = plot_width*0.3, y = plot_height*0.5)
  
  # Create a data frame for Bud hurricane icon
  bud_df <- data.frame(
    x = plot_width*0.3, 
    y = plot_height*0.75, 
    image = b$image[1]  # Use the original hurricane image
  )
  
  # Create a data frame for the circle points
  theta <- seq(0, 2*pi, length.out=100)
  circle_size <- 20
  circle_df <- data.frame(
    x = plot_width*0.3 + circle_size*cos(theta),
    y = plot_height*0.25 + circle_size*sin(theta),
    group = rep(1, 100)
  )
  
  # Create a data frame for the outer circle (white)
  outer_circle_df <- data.frame(
    x = plot_width*0.3 + (circle_size+2)*cos(theta),
    y = plot_height*0.25 + (circle_size+2)*sin(theta),
    group = rep(1, 100)
  )
  
  # Create the plot
  p <- ggplot() +
    # Add the Aletta eye (triangle)
    geom_point(
      data = aletta_df, 
      aes(x = x, y = y),
      shape = aletta_symbol,  # Use original symbol
      color = aletta_color,   # Use original color
      fill = aletta_color,    # Use original fill
      size = aletta_size * 2  # Size is still valid for points
    ) +
    
    # Add Bud hurricane image
    geom_image(
      data = bud_df,
      aes(x = x, y = y, image = image),
      size = 0.1  # Size is valid for geom_image
    ) +
    
    # Add Station 3.5 circle - dashed orange line
    # Changed from size to linewidth
    geom_path(
      data = circle_df,
      aes(x = x, y = y, group = group),
      color = circle_color2,  # Use original orange color
      linetype = "dashed",
      linewidth = base * 0.4  # Changed from size to linewidth
    ) +
    
    # Add outer white circle
    # Changed from size to linewidth
    geom_path(
      data = outer_circle_df,
      aes(x = x, y = y, group = group),
      color = circle_color1,  # Use original white color
      linewidth = outline * 0.4  # Changed from size to linewidth
    ) +
    
    # Add text labels
    annotate(
      "text", 
      x = plot_width*0.5, 
      y = plot_height*0.75, 
      label = "TC Bud Eye", 
      hjust = 0,
      size = 3
    ) +
    
    annotate(
      "text", 
      x = plot_width*0.5, 
      y = plot_height*0.5, 
      label = "TC Aletta Eye", 
      hjust = 0,
      size = 3
    ) +
    
    annotate(
      "text", 
      x = plot_width*0.5, 
      y = plot_height*0.25, 
      label = "200km Diameter\nStation 3.5", 
      hjust = 0,
      size = 3
    ) +
    
    # Clean theme
    theme_void() +
    
    # Set plot dimensions
    xlim(0, plot_width) +
    ylim(0, plot_height) +
    
    # Make sure no background elements show
    theme(
      panel.background = element_blank(),
      plot.background = element_blank()
    )
  
  return(p)
}

# Usage remains the same:
# refined_legend <- f.refined_eye_legend()
# 
# setwd(wd)
# setwd(fig)
# pdf("figure_two_refined_eye_legend.pdf", width=3, height=4)
# print(refined_legend)
# dev.off()

# Test the refined legend
refined_legend <- f.refined_eye_legend()

# Save as a test
setwd(wd)
setwd(fig)
pdf("figure_two_refined_eye_legend.pdf", width=3, height=4)
print(refined_legend)
dev.off()

# Create the main figure with just the color legend
t<--390
legend <- f.fig2_legend(combined, v_rgb, height=bar_height, width=1, t=t, b=0, r=0, l=-115, 
                        title.textsize=10, marks.textsize=8)

# Save main figure with color legend
pdf("figure_two_main_with_color_legend.pdf", height=11, width=8.5)
main_plot <- plot_grid(p, legend, ncol=2, rel_widths=c(14,1), align="h")
print(main_plot)
dev.off()

# Try a direct approach with both legends combined
pdf("figure_two_with_both_legends.pdf", height=11, width=8.5)
# First combine the two legends
combined_legends <- tryCatch({
  plot_grid(
    legend,
    refined_legend,
    ncol = 1,
    rel_heights = c(4, 1),
    align = "v"
  )
}, error = function(e) {
  # If combining fails, return NULL
  message("Error combining legends: ", e$message)
  return(NULL)
})

# If combination succeeded, create the full figure
if (!is.null(combined_legends)) {
  full_plot <- plot_grid(p, combined_legends, ncol=2, rel_widths=c(14,1), align="h")
  print(full_plot)
} else {
  # Otherwise just print the main plot with color legend
  print(main_plot)
}
dev.off()

# If the above doesn't work, create a version of the eye legend
# that can be manually pasted into your layout

# Create a standalone eye legend with proper dimensions for manual placement
# pdf("figure_two_eye_legend_for_manual_placement.pdf", width=1.5, height=3)
# standalone_legend <- f.refined_eye_legend(plot_width=150, plot_height=300)
# print(standalone_legend)
# dev.off()


#### new ####
## Create a function specifically designed to produce a vertically stacked legend
## This will create the legend as a single element rather than trying to combine two plots
f.create_combined_legend <- function() {
  # First, create the color scale legend using your existing function
  t <- -390
  color_legend <- f.fig2_legend(combined, v_rgb, height=bar_height, width=1, t=t, b=0, r=0, l=-115, 
                                title.textsize=10, marks.textsize=8)
  
  # Get the pixel dimensions of the color legend
  # (This is just an estimate since we can't directly measure the plot dimensions)
  color_legend_height <- 11 * 72 * 0.8  # 80% of 11 inches in points
  color_legend_width <- 1.5 * 72        # 1.5 inches in points
  
  # Create a vertical eye legend with matching width
  # Create a data frame for Aletta eye (triangle)
  plot_width <- color_legend_width
  plot_height <- 11 * 72 * 0.2          # 20% of 11 inches in points
  
  # Generate the eye legend at the bottom
  symbol_x <- plot_width * 0.3
  label_x <- plot_width * 0.5
  spacing <- plot_height / 4
  
  # Y positions from bottom to top
  circle_y <- spacing * 1
  aletta_y <- spacing * 2
  bud_y <- spacing * 3
  
  # Create data frames for each element
  aletta_df <- data.frame(x = symbol_x, y = aletta_y)
  
  bud_df <- data.frame(
    x = symbol_x, 
    y = bud_y, 
    image = b$image[1]
  )
  
  # Create the circle points
  theta <- seq(0, 2*pi, length.out=100)
  circle_size <- spacing * 0.4
  circle_df <- data.frame(
    x = symbol_x + circle_size*cos(theta),
    y = circle_y + circle_size*sin(theta),
    group = rep(1, 100)
  )
  
  # Create the outer circle (white outline)
  outer_circle_df <- data.frame(
    x = symbol_x + (circle_size+2)*cos(theta),
    y = circle_y + (circle_size+2)*sin(theta),
    group = rep(1, 100)
  )
  
  # Create eye legend
  eye_legend <- ggplot() +
    # Add the Aletta eye (triangle)
    geom_point(
      data = aletta_df, 
      aes(x = x, y = y),
      shape = aletta_symbol, 
      color = aletta_color, 
      fill = aletta_color, 
      size = aletta_size * 2
    ) +
    
    # Add Bud hurricane image
    geom_image(
      data = bud_df,
      aes(x = x, y = y, image = image),
      size = 0.1
    ) +
    
    # Add Station 3.5 circle - dashed orange line
    geom_path(
      data = circle_df,
      aes(x = x, y = y, group = group),
      color = circle_color2,
      linetype = "dashed",
      linewidth = base * 0.4
    ) +
    
    # Add outer white circle
    geom_path(
      data = outer_circle_df,
      aes(x = x, y = y, group = group),
      color = circle_color1,
      linewidth = outline * 0.4
    ) +
    
    # Add text labels with appropriate spacing
    annotate(
      "text", 
      x = label_x, 
      y = bud_y, 
      label = "TC Bud Eye", 
      hjust = 0,
      size = 3
    ) +
    
    annotate(
      "text", 
      x = label_x, 
      y = aletta_y, 
      label = "TC Aletta Eye", 
      hjust = 0,
      size = 3
    ) +
    
    annotate(
      "text", 
      x = label_x, 
      y = circle_y, 
      label = "200km Diameter\nStation 3.5", 
      hjust = 0,
      size = 3
    ) +
    
    # Clean theme
    theme_void() +
    
    # Set plot dimensions to match the width of the color legend
    xlim(0, plot_width) +
    ylim(0, plot_height) +
    
    # Make sure no background elements show
    theme(
      panel.background = element_blank(),
      plot.background = element_blank(),
      plot.margin = margin(0, 0, 0, 0)
    )
  
  # Combine the two legends with appropriate heights
  combined <- plot_grid(
    color_legend,
    eye_legend,
    ncol = 1,
    rel_heights = c(4, 1),  # Adjust this ratio as needed
    align = "v"
  )
  
  return(combined)
}

# Use the function to create a combined legend
combined_legend <- f.create_combined_legend()

# Save the combined legend separately for testing
setwd(wd)
setwd(fig)
pdf("figure_two_combined_legend_test.pdf", width=2, height=11)
print(combined_legend)
dev.off()

# Create the full figure with the combined legend
pdf("figure_two_with_properly_combined_legends.pdf", height=11, width=8.5)
full_plot <- plot_grid(p, combined_legend, ncol=2, rel_widths=c(14,1), align="h")
print(full_plot)
dev.off()

# Alternative implementation with explicit sizing for each element
# This version creates each component with very precise control
f.create_explicit_legend <- function() {
  # Create separate plots for each element
  
  # 1. Color scale
  t <- -390
  color_legend <- f.fig2_legend(combined, v_rgb, height=bar_height, width=1, t=t, b=0, r=0, l=-115, 
                                title.textsize=10, marks.textsize=8)
  
  # 2. Bud symbol with label
  bud_plot <- ggplot() +
    geom_image(
      data = data.frame(x=0.3, y=0.5, image=b$image[1]),
      aes(x=x, y=y, image=image),
      size = 0.1
    ) +
    annotate("text", x=0.6, y=0.5, label="TC Bud Eye", hjust=0, size=3) +
    theme_void() +
    xlim(0, 2) + ylim(0, 1) +
    theme(plot.margin = margin(5, 0, 5, 0))
  
  # 3. Aletta symbol with label
  aletta_plot <- ggplot() +
    geom_point(
      data = data.frame(x=0.3, y=0.5),
      aes(x=x, y=y),
      shape = aletta_symbol,
      color = aletta_color,
      fill = aletta_color,
      size = aletta_size * 2
    ) +
    annotate("text", x=0.6, y=0.5, label="TC Aletta Eye", hjust=0, size=3) +
    theme_void() +
    xlim(0, 2) + ylim(0, 1) +
    theme(plot.margin = margin(5, 0, 5, 0))
  
  # 4. Circle with label
  circle_plot <- ggplot() +
    # Dashed circle
    geom_path(
      data = data.frame(
        x = 0.3 + 0.2*cos(seq(0, 2*pi, length.out=100)),
        y = 0.5 + 0.2*sin(seq(0, 2*pi, length.out=100)),
        group = 1
      ),
      aes(x=x, y=y, group=group),
      color = circle_color2,
      linetype = "dashed",
      linewidth = base * 0.4
    ) +
    # White outline
    geom_path(
      data = data.frame(
        x = 0.3 + 0.22*cos(seq(0, 2*pi, length.out=100)),
        y = 0.5 + 0.22*sin(seq(0, 2*pi, length.out=100)),
        group = 1
      ),
      aes(x=x, y=y, group=group),
      color = circle_color1,
      linewidth = outline * 0.4
    ) +
    annotate("text", x=0.6, y=0.5, label="200km Diameter\nStation 3.5", hjust=0, size=3) +
    theme_void() +
    xlim(0, 2) + ylim(0, 1) +
    theme(plot.margin = margin(5, 0, 5, 0))
  
  # Combine the eye legend elements vertically
  eye_legend <- plot_grid(
    bud_plot,
    aletta_plot,
    circle_plot,
    ncol = 1,
    align = "v"
  )
  
  # Now combine with the color legend
  combined <- plot_grid(
    color_legend,
    eye_legend,
    ncol = 1,
    rel_heights = c(4, 1),
    align = "v"
  )
  
  return(combined)
}

# Try the explicit legend approach
explicit_legend <- f.create_explicit_legend()

pdf("figure_two_with_explicit_legend.pdf", height=11, width=8.5)
plot_grid(p, explicit_legend, ncol=2, rel_widths=c(14,1), align="h")
dev.off()







### newer ####
###################################################
# Create a clean, simple solution for combining legends
###################################################

# First, create a function that only generates the eye legend 
# without any color legend components
f.eye_legend_only <- function() {
  # Create a simple plot with the three elements
  
  # Set basic sizing parameters
  plot_width <- 2 * 72  # 2 inches in points
  plot_height <- 3 * 72  # 3 inches in points
  symbol_x <- 0.3 * plot_width
  text_x <- 0.55 * plot_width
  
  # Vertical positions - evenly spaced
  y1 <- 0.75 * plot_height  # Bud position
  y2 <- 0.5 * plot_height   # Aletta position
  y3 <- 0.25 * plot_height  # Circle position
  circle_size <- 0.1 * plot_height
  
  # Create the plot
  p <- ggplot() +
    # Bud hurricane image
    geom_image(
      data = data.frame(x = symbol_x, y = y1, image = b$image[1]),
      aes(x = x, y = y, image = image),
      size = 0.1
    ) +
    
    # Aletta triangle
    geom_point(
      data = data.frame(x = symbol_x, y = y2),
      aes(x = x, y = y),
      shape = aletta_symbol,
      color = aletta_color,
      fill = aletta_color,
      size = aletta_size * 2
    ) +
    
    # Station circle - dashed orange
    geom_path(
      data = data.frame(
        x = symbol_x + circle_size*cos(seq(0, 2*pi, length.out=100)),
        y = y3 + circle_size*sin(seq(0, 2*pi, length.out=100)),
        group = 1
      ),
      aes(x = x, y = y, group = group),
      color = circle_color2,
      linetype = "dashed",
      linewidth = base * 0.4
    ) +
    
    # Station circle - white outline
    geom_path(
      data = data.frame(
        x = symbol_x + (circle_size*1.1)*cos(seq(0, 2*pi, length.out=100)),
        y = y3 + (circle_size*1.1)*sin(seq(0, 2*pi, length.out=100)),
        group = 1
      ),
      aes(x = x, y = y, group = group),
      color = circle_color1,
      linewidth = outline * 0.4
    ) +
    
    # Text labels
    annotate("text", x = text_x, y = y1, label = "TC Bud Eye", hjust = 0, size = 3) +
    annotate("text", x = text_x, y = y2, label = "TC Aletta Eye", hjust = 0, size = 3) +
    annotate("text", x = text_x, y = y3, label = "200km Diameter\nStation 3.5", hjust = 0, size = 3) +
    
    # Clean theme and proper dimensions
    theme_void() +
    xlim(0, plot_width) +
    ylim(0, plot_height) +
    theme(
      panel.background = element_blank(),
      plot.background = element_blank(),
      plot.margin = margin(0, 0, 0, 0)
    )
  
  return(p)
}

#################################################
# Main implementation - 2-step approach
#################################################

# STEP 1: Create and save the main figure with color legend only
setwd(wd)
setwd(fig)

# First get the color legend
t <- -390
color_legend <- f.fig2_legend(combined, v_rgb, height=bar_height, width=1, t=t, b=0, r=0, l=-115, 
                              title.textsize=10, marks.textsize=8)

# Save the main figure with just the color legend
pdf("figure_two_main_with_color_legend.pdf", height=11, width=8.5)
main_plot <- plot_grid(p, color_legend, ncol=2, rel_widths=c(14,1), align="h")
print(main_plot)
dev.off()

# STEP 2: Create and save just the eye legend separately
# This produces a standalone legend with proper dimensions
pdf("figure_two_eye_legend_final.pdf", width=2, height=3)
eye_legend <- f.eye_legend_only()
print(eye_legend)
dev.off()

# Optional: Try a direct combined approach if desired
# This may or may not work depending on memory constraints
tryCatch({
  # Run garbage collection to free memory
  gc()
  
  # Create the eye legend
  eye_legend <- f.eye_legend_only()
  
  # Try to combine with the color legend in a specific layout
  # First create a blank padding element
  blank_plot <- ggplot() + theme_void()
  
  # Create a layout that puts the eye legend at the bottom right
  layout <- rbind(
    c(1, 2),
    c(1, 2),
    c(1, 2),
    c(1, 2),
    c(1, 3)
  )
  
  # Use this layout to arrange the plots
  final_combined <- grid.arrange(
    p, color_legend, eye_legend,
    layout_matrix = layout
  )
  
  # Save the result
  pdf("figure_two_final_combined_layout.pdf", height=11, width=8.5)
  print(final_combined)
  dev.off()
}, error = function(e) {
  message("Error creating combined layout: ", e$message)
  message("Please use the separate PDF files to manually combine the legends.")
})