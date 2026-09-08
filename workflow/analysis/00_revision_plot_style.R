# Central plotting style for PFAS mSystems revision.
# Purple-focused palette; source this file before making figures.

pfas_palette <- list(
  purple_dark   = "#3B0F70",
  purple_main   = "#6A00A8",
  purple_mid    = "#8C2981",
  purple_light  = "#B12A90",
  mauve         = "#CC79A7",
  lavender      = "#D7B5F7",
  grey_dark     = "#333333",
  grey_mid      = "#777777",
  grey_light    = "#D9D9D9",
  black         = "#000000"
)

pfas_discrete <- c(
  pfas_palette$purple_dark,
  pfas_palette$purple_main,
  pfas_palette$purple_mid,
  pfas_palette$purple_light,
  pfas_palette$mauve,
  pfas_palette$lavender,
  pfas_palette$grey_mid
)

pfas_exposure <- c(
  "Reference" = pfas_palette$grey_mid,
  "Low"       = pfas_palette$lavender,
  "Medium"    = pfas_palette$purple_light,
  "High"      = pfas_palette$purple_dark
)

pfas_theme_base <- function(base_size=10) {
  if (!requireNamespace("ggplot2", quietly=TRUE)) {
    stop("ggplot2 is required for pfas_theme_base()")
  }

  ggplot2::theme_classic(base_size=base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(color=pfas_palette$grey_dark),
      axis.text = ggplot2::element_text(color=pfas_palette$grey_dark),
      axis.title = ggplot2::element_text(color=pfas_palette$grey_dark),
      plot.title = ggplot2::element_text(color=pfas_palette$grey_dark, face="bold"),
      legend.title = ggplot2::element_text(color=pfas_palette$grey_dark),
      legend.text = ggplot2::element_text(color=pfas_palette$grey_dark)
    )
}
