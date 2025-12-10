\name{plot_track_landscape}
\alias{plot_track_landscape}
\title{Plot scaled multi-track genomic landscapes}
\usage{
plot_track_landscape(
  bed_files,
  target_chr,
  genome_size,
  bin_size = 100000,
  track_names = NULL,
  track_colors = NULL,
  track_scale_ratios = NULL,
  custom_y_limits = NULL,
  alpha = 0.9,
  title = NULL,
  save_path = NULL,
  width = 15,
  height = 6,
  dpi = 300
)
}
\arguments{
\item{bed_files}{Named character vector of BED file paths, or a named list of data.frames with at least \code{chr}, \code{start}, and \code{end} columns.}

\item{target_chr}{Chromosome/contig identifier to plot.}

\item{genome_size}{Length of the chromosome in base pairs.}

\item{bin_size}{Bin width in base pairs (default 100000).}

\item{track_names}{Optional labels matching \code{bed_files}. Defaults to names(bed_files).}

\item{track_colors}{Named vector of colors for each track. If \code{NULL}, a qualitative palette is generated.}

\item{track_scale_ratios}{Named numeric vector of multiplicative factors applied to each track's binned counts. Missing tracks default to 1.}

\item{custom_y_limits}{Named numeric vector giving the maximum y for each track (applied after scaling).}

\item{alpha}{Fill alpha for the area geometry (default 0.9).}

\item{title}{Optional plot title.}

\item{save_path}{Optional file path to save the plot via \code{ggsave}.}

\item{width, height, dpi}{Export size/dpi passed to \code{ggsave} when \code{save_path} is provided.}
}
\value{
ggplot object representing the faceted landscape.
}
\description{
Bins genomic intervals per track along a target chromosome, applies optional proportional scaling, and renders aligned area plots with publication-ready defaults.
}
\examples{
data("sample_tracks", package = "TrackLandscapeR")
p <- plot_track_landscape(
  bed_files = sample_tracks,
  target_chr = "Chr1",
  genome_size = 54e6,
  bin_size = 1e5,
  track_colors = c(Copia = "#0072B2", TRASH = "#D92121", Gypsy = "#CC79A7"),
  track_scale_ratios = c(Copia = 2.5, TRASH = 0.03, Gypsy = 2.5),
  custom_y_limits = c(Copia = 150, TRASH = 70, Gypsy = 150)
)
}

