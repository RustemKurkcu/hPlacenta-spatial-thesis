# Centralized aesthetics for thesis spatial plots

misi_celltype_colors <- c(
  "vCTB2" = "#E6A0C4",
  "Fibroblast2" = "#FBE735",
  "Fibroblast1" = "#F28239",
  "STB" = "#92D3E3",
  "Hofbauer cells" = "#8C9ECF",
  "EVT2" = "#006400",
  "Endothelial_cells" = "#D42F2E",
  "EVT3" = "#228B22",
  "maternal macrophages" = "#6A3D9A",
  "vCTBp" = "#F4CAE4",
  "unknown2" = "#BDBDBD",
  "unknown3" = "#969696"
)

misi_week_colors <- c(
  "W7" = "#1B9E77",
  "W8-2" = "#D95F02",
  "W9" = "#7570B3",
  "W11" = "#E7298A"
)

theme_thesis_spatial <- function(base_size = 11, legend_position = "right") {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", linewidth = 0.4),
      legend.position = legend_position,
      legend.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      strip.background = ggplot2::element_rect(fill = "grey95", color = "grey70"),
      strip.text = ggplot2::element_text(face = "bold")
    )
}
