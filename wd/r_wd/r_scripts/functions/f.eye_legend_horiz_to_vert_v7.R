#Function to create a horizontal or vertical eye legend - Version 7
# Dynamically adjusts dimensions based on orientation
f.eye_legend_horiz_to_vert_v7 <- function(
    # Basic dimensions with dynamic orientation handling
  long_side_dimension = 800,   # Longer side of the plot
  short_side_dimension = 100,  # Shorter side of the plot
  is_vertical = FALSE,         # Control orientation
  
  # Spacing parameters
  left_margin = 40,         # Space at left edge
  right_margin = 40,        # Space at right edge
  icon_text_gap = 20,       # Space between icon and its text
  element_gap = 60,         # Space between text and next icon
  
  # Flexible spacing ratios
  margin_ratio = NULL,      # Optional ratio-based margin calculation
  inter_element_ratio = NULL, # Optional ratio-based element spacing
  
  # Sizing parameters
  text_size = 10,           # Base text size
  circle_size = 20,         # Radius of the station circle
  aletta_scale = 2.5,       # Scale factor for Aletta symbol
  bud_scale = 0.18,         # Scale factor for Bud hurricane image
  icon_scale_factor = 1,    # Global scaling for all icons
  
  # Circle styling parameters
  circle_parameters = NULL, # Circle styling parameters
  orange_circle_dash_length = 1  # Parameter to control dash spacing
) {
  # Helper function for null values
  "%||%" <- function(x, y) {
    if (is.null(x)) y else x
  }
  
  # Dynamically set plot dimensions based on orientation
  if (is_vertical) {
    plot_width <- short_side_dimension
    plot_height <- long_side_dimension
  } else {
    plot_width <- long_side_dimension
    plot_height <- short_side_dimension
  }
  
  # Calculate margins using ratio if provided
  if (!is.null(margin_ratio)) {
    total_margin_width <- plot_width * margin_ratio
    left_margin <- total_margin_width / 2
    right_margin <- total_margin_width / 2
  }
  
  # Calculate element gap using ratio if provided
  if (!is.null(inter_element_ratio)) {
    element_gap <- plot_width * inter_element_ratio
  }
  
  # Adjust scaling for icons
  adjusted_circle_size <- circle_size * icon_scale_factor
  adjusted_aletta_scale <- aletta_scale * icon_scale_factor
  adjusted_bud_scale <- bud_scale * icon_scale_factor
  
  # Calculate usable width (accounting for margins)
  usable_width <- plot_width - left_margin - right_margin
  
  # Extract styling from circle_parameters
  if (!is.null(circle_parameters)) {
    circle_color1 <- circle_parameters["circle_color1"] %||% "white"
    circle_color2 <- circle_parameters["circle_color2"] %||% "#E65832"
    base <- as.numeric(circle_parameters["base"] %||% "5.0")
    outline <- as.numeric(circle_parameters["outline"] %||% "7.5")
  } else {
    # Default values if circle_parameters is NULL
    circle_color1 <- "white"
    circle_color2 <- "#E65832"
    base <- 5.0
    outline <- 7.5
  }
  
  # Set line thicknesses for legend
  legend_outline <- outline * 0.5
  legend_base <- base * 0.5
  
  # Text labels
  text1 <- "TC Bud Eye"
  text2 <- "TC Aletta Eye"
  text3 <- "200km Diameter\nStation 3.5"
  
  # Create base plot
  p <- ggplot() +
    # Clean theme with minimal margins
    theme_void() +
    theme(
      plot.margin = margin(0, 0, 0, 0),
      panel.border = element_blank(),
      panel.spacing = unit(0, "mm")
    ) +
    coord_cartesian(xlim = c(0, plot_width), ylim = c(0, plot_height), expand = FALSE)
  
  # Layout parameters for single-row approach
  x_spacing <- plot_width / 6  # Divide width into 6 equal segments
  y_center <- plot_height / 2
  
  # Horizontal layout
  if (!is_vertical) {
    p <- p +
      # Bud Eye
      geom_image(data = data.frame(x = x_spacing, y = y_center, 
                                   image = b$image[1]), 
                 aes(x, y, image = image), size = adjusted_bud_scale) +
      # Bud Eye Text
      annotate("text", x = x_spacing * 2, y = y_center, label = text1, 
               hjust = 0, size = text_size/2.5) +
      
      # Aletta Eye
      geom_point(data = data.frame(x = x_spacing * 3, y = y_center), 
                 aes(x, y), shape = aletta_symbol, color = aletta_color, 
                 fill = aletta_color, size = aletta_size * adjusted_aletta_scale) +
      # Aletta Eye Text
      annotate("text", x = x_spacing * 4, y = y_center, label = text2, 
               hjust = 0, size = text_size/2.5) +
      
      # Station Circle with white and orange outlines
      geom_path(data = data.frame(
        x = x_spacing * 5 + adjusted_circle_size * cos(seq(0, 2*pi, length.out = 100)),
        y = y_center + adjusted_circle_size * sin(seq(0, 2*pi, length.out = 100))
      ), aes(x, y), color = circle_color1, linewidth = legend_outline, 
      linetype = "solid", alpha = 1, inherit.aes = FALSE) +
      geom_path(data = data.frame(
        x = x_spacing * 5 + adjusted_circle_size * cos(seq(0, 2*pi, length.out = 100)),
        y = y_center + adjusted_circle_size * sin(seq(0, 2*pi, length.out = 100))
      ), aes(x, y), color = circle_color2, linewidth = legend_base, 
      linetype = orange_circle_dash_length, alpha = 1, inherit.aes = FALSE) +
      # Station Circle Text
      annotate("text", x = x_spacing * 6, y = y_center, label = text3, 
               hjust = 0, size = text_size/2.5)
  } else {
    # Vertical layout (single row with rotated text)
    p <- p +
      # Bud Eye
      geom_image(data = data.frame(x = x_spacing, y = y_center, 
                                   image = b$image[1]), 
                 aes(x, y, image = image), size = adjusted_bud_scale) +
      # Bud Eye Text (rotated)
      annotate("text", x = x_spacing * 2, y = y_center, label = text1, 
               hjust = 1, angle = 90, size = text_size/2.5) +
      
      # Aletta Eye
      geom_point(data = data.frame(x = x_spacing * 3, y = y_center), 
                 aes(x, y), shape = aletta_symbol, color = aletta_color, 
                 fill = aletta_color, size = aletta_size * adjusted_aletta_scale) +
      # Aletta Eye Text (rotated)
      annotate("text", x = x_spacing * 4, y = y_center, label = text2, 
               hjust = 1, angle = 90, size = text_size/2.5) +
      
      # Station Circle with white and orange outlines
      geom_path(data = data.frame(
        x = x_spacing * 5 + adjusted_circle_size * cos(seq(0, 2*pi, length.out = 100)),
        y = y_center + adjusted_circle_size * sin(seq(0, 2*pi, length.out = 100))
      ), aes(x, y), color = circle_color1, linewidth = legend_outline, 
      linetype = "solid", alpha = 1, inherit.aes = FALSE) +
      geom_path(data = data.frame(
        x = x_spacing * 5 + adjusted_circle_size * cos(seq(0, 2*pi, length.out = 100)),
        y = y_center + adjusted_circle_size * sin(seq(0, 2*pi, length.out = 100))
      ), aes(x, y), color = circle_color2, linewidth = legend_base, 
      linetype = orange_circle_dash_length, alpha = 1, inherit.aes = FALSE) +
      # Station Circle Text (rotated)
      annotate("text", x = x_spacing * 6, y = y_center, label = text3, 
               hjust = 1, angle = 90, size = text_size/2.5)
  }
  
  return(p)
}
Last edited 6 minutes ago
