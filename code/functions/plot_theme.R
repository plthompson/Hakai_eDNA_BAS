plt_theme <- theme_bw() + theme(
  plot.background = element_blank(),#element_rect(),    # Background of the entire plot

  panel.background = element_rect(),   # Background of plotting area
  panel.border = element_rect(),       # Border around plotting area.
  # fill argument should be NA

  panel.grid = element_blank(),         # All grid lines
  panel.grid.major = element_blank(),   # Major grid lines
  panel.grid.minor = element_blank(),   # Minor grid lines

  panel.grid.major.x = element_blank(), # Vertical major grid lines
  panel.grid.major.y = element_blank(), # Horizontal major grid lines
  panel.grid.minor.x = element_blank(), # Vertical minor grid lines
  panel.grid.minor.y = element_blank(),  # Vertical major grid lines

  strip.background = element_rect(colour=NA, fill=NA),
  strip.text.x = element_text(size = 12),
  strip.text.y = element_text(size = 12),
  axis.title=element_text(size = 12),
  plot.title = element_text(hjust = 0.5),
  legend.background = element_rect(fill= NA)
)
