# TrackLandscapeR

TrackLandscapeR builds publication-ready, multi-track density landscapes from genomic interval (BED-like) files. It bins intervals along a target chromosome, applies per-track proportional scaling, and renders aligned area plots with consistent x-axes and per-track y-limits.

## Key features
- Single entry point `plot_track_landscape()` to process multiple tracks and plot a faceted landscape.
- Per-track proportional scaling and optional custom y-axis limits.
- Publication-oriented defaults (minimal theme, strip placement, axis labeling in Mb).
- Included `sample_tracks` dataset and a vignette to reproduce the demo figure.

## Installation
```r
# from a local checkout
install.packages("TrackLandscapeR", repos = NULL, type = "source")
```

## Quick start
```r
library(TrackLandscapeR)

data("sample_tracks", package = "TrackLandscapeR")

p <- plot_track_landscape(
  bed_files = sample_tracks,
  target_chr = "Chr1",
  genome_size = 54e6,
  bin_size = 1e5,
  track_colors = c(Copia = "#0072B2", TRASH = "#D92121", Gypsy = "#CC79A7"),
  track_scale_ratios = c(Copia = 2.5, TRASH = 0.03, Gypsy = 2.5),
  custom_y_limits = c(Copia = 150, TRASH = 70, Gypsy = 150),
  save_path = "sample_output.pdf"
)

print(p)
```

## Run on your real data (with dashed lines if desired)
```r
# install from local source if not yet installed
install.packages("/media/liurxx/华为硬盘/30.igv/1.TRASH/all/TrackLandscapeR",
                 repos = NULL, type = "source")

library(TrackLandscapeR)

bed_files <- c(
  Copia = "Chr1.Copia.Glygla.bed",
  TRASH = "Chr1_Glygla_TRASH_CEN60.bed",
  Gypsy = "Chr1.Gypsy.Glygla.bed"
)

p <- plot_track_landscape(
  bed_files = bed_files,
  target_chr = "Chr1",
  genome_size = 53933422,
  bin_size = 100000,
  track_colors = c(Copia="#0072B2", TRASH="#D92121", Gypsy="#CC79A7"),
  track_scale_ratios = c(Copia=2.5, TRASH=0.03, Gypsy=2.5),
  custom_y_limits = c(Copia=150, TRASH=70, Gypsy=150),
  save_path = "/media/liurxx/华为硬盘/30.igv/1.TRASH/all/output_multitrack_scaled_actual.pdf",
  width = 15, height = 6, dpi = 300
)

# optional dashed reference lines
p <- p +
  ggplot2::geom_vline(xintercept = c(10e6, 30e6), linetype = "dashed", color = "gray40", linewidth = 0.4) +
  ggplot2::geom_hline(yintercept = 50, linetype = "dashed", color = "gray40", linewidth = 0.3)

print(p)
ggplot2::ggsave("output_with_dashed_lines.pdf", p, width = 15, height = 6, dpi = 300)
```

## Reproducibility
The vignette (`vignettes/TrackLandscapeR.Rmd`) walks through the workflow, including binning logic, scaling, and figure export, using the bundled sample data.

