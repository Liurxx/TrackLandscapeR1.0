# Creates a lightweight example dataset derived from the three provided BED files.
# The script samples the first 200 entries of each track to keep package size modest.

library(utils)

root_dir <- "/media/liurxx/华为硬盘/30.igv/1.TRASH/all"

read_subset <- function(path, n = 200) {
  tbl <- read.table(path, header = FALSE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE, quote = "")
  cols_to_name <- c("chr", "start", "end", "seq", "score", "strand", "thickStart", "thickEnd", "rgb")
  num_cols <- min(ncol(tbl), length(cols_to_name))
  colnames(tbl)[seq_len(num_cols)] <- cols_to_name[seq_len(num_cols)]
  head(tbl, n)
}

sample_tracks <- list(
  Copia = read_subset(file.path(root_dir, "Chr1.Copia.Glygla.bed")),
  TRASH = read_subset(file.path(root_dir, "Chr1_Glygla_TRASH_CEN60.bed")),
  Gypsy = read_subset(file.path(root_dir, "Chr1.Gypsy.Glygla.bed"))
)

dir.create(file.path(root_dir, "TrackLandscapeR", "data"), showWarnings = FALSE, recursive = TRUE)
save(sample_tracks, file = file.path(root_dir, "TrackLandscapeR", "data", "sample_tracks.rda"), compress = "xz")

