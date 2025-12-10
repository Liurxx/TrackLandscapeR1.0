\name{sample_tracks}
\alias{sample_tracks}
\docType{data}
\title{Example genomic interval tracks}
\format{
A named list of three data.frames (\code{Copia}, \code{TRASH}, \code{Gypsy}).
Each data.frame contains up to nine BED-like columns: \code{chr}, \code{start},
\code{end}, \code{seq}, \code{score}, \code{strand}, \code{thickStart},
\code{thickEnd}, and \code{rgb}. Only \code{chr}, \code{start}, and \code{end}
are required for binning.
}
\usage{
sample_tracks
}
\description{
Lightweight subsets (first 200 rows) from the three BED files used in the demo.
They enable running \code{plot_track_landscape()} without external downloads.
}
\examples{
data("sample_tracks", package = "TrackLandscapeR")
str(sample_tracks)
}

