#' Plot scaled multi-track genomic landscapes
#'
#' Bins genomic intervals per track along a target chromosome, applies optional
#' per-track proportional scaling, and renders aligned area plots with
#' publication-ready defaults.
#'
#' @param bed_files Named character vector of BED/bed-like file paths, or a
#'   named list of data.frames containing at least `chr`, `start`, and `end`
#'   columns. Names are used as track labels.
#' @param target_chr Chromosome/contig identifier to plot.
#' @param genome_size Length of the chromosome (numeric, in bp).
#' @param bin_size Bin width in base pairs. Default: 100000.
#' @param track_names Optional vector of labels matching `bed_files`.
#' @param track_colors Named vector of colors for each track. If `NULL`, a
#'   qualitative palette is generated.
#' @param track_scale_ratios Named numeric vector of multiplicative factors
#'   applied to each track's binned counts. Missing tracks default to 1.
#' @param custom_y_limits Named numeric vector giving the maximum visible y for
#'   each track (applied after scaling).
#' @param alpha Fill alpha for the area geometry. Default: 0.9.
#' @param title Plot title. If `NULL`, a title is generated from the chromosome.
#' @param save_path Optional path to save the plot (pdf/png compatible with
#'   `ggsave`).
#' @param width,height,dpi Plot export size/dpi passed to `ggsave`.
#'
#' @return A ggplot object.
#' @export
plot_track_landscape <- function(
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
) {
  if (is.null(track_names)) {
    track_names <- names(bed_files)
  }
  if (is.null(track_names) || length(track_names) != length(bed_files)) {
    stop("`track_names` must be provided and match the length of `bed_files`.")
  }
  if (is.null(names(bed_files))) {
    names(bed_files) <- track_names
  }

  # helper: read path or accept data.frame
  normalize_input <- function(x) {
    if (is.character(x) && length(x) == 1) {
      utils::read.table(x, header = FALSE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE, quote = "")
    } else if (is.data.frame(x)) {
      x
    } else {
      stop("Each element of `bed_files` must be a file path or a data.frame.")
    }
  }

  # helper: ensure standard column names
  standardize_cols <- function(df) {
    cols_to_name <- c("chr", "start", "end", "seq", "score", "strand", "thickStart", "thickEnd", "rgb")
    num_cols <- min(ncol(df), length(cols_to_name))
    colnames(df)[seq_len(num_cols)] <- cols_to_name[seq_len(num_cols)]
    df
  }

  # helper: bin one track
  bin_track <- function(df, track_label) {
    df <- standardize_cols(df)
    if (!all(c("chr", "start", "end") %in% colnames(df))) {
      stop("Input must contain at least chr, start, end columns.")
    }
    df <- df[df$chr == target_chr, , drop = FALSE]
    if (!nrow(df)) {
      return(NULL)
    }
    df$start <- suppressWarnings(as.numeric(df$start))
    df$end <- suppressWarnings(as.numeric(df$end))
    df$mid <- (df$start + df$end) / 2
    df <- df[!is.na(df$mid), , drop = FALSE]

    max_pos <- max(df$mid, na.rm = TRUE)
    current_max <- if (max_pos > genome_size) max_pos else genome_size
    num_bins <- ceiling(current_max / bin_size)
    breaks <- seq(0, num_bins * bin_size, by = bin_size)
    bin_counts <- hist(df$mid, breaks = breaks, plot = FALSE)

    out <- data.frame(
      bin_start = bin_counts$breaks[-length(bin_counts$breaks)],
      freq = bin_counts$counts,
      track = track_label,
      stringsAsFactors = FALSE
    )
    out[out$bin_start < genome_size, , drop = FALSE]
  }

  data_list <- mapply(function(x, lab) {
    tab <- normalize_input(x)
    bin_track(tab, lab)
  }, bed_files, track_names, SIMPLIFY = FALSE)

  plot_data <- dplyr::bind_rows(data_list)
  if (!nrow(plot_data)) {
    stop("No data available after filtering for target_chr.")
  }
  plot_data$track <- factor(plot_data$track, levels = track_names)

  # apply scaling
  if (!is.null(track_scale_ratios)) {
    for (nm in names(track_scale_ratios)) {
      ratio <- track_scale_ratios[[nm]]
      plot_data$freq[plot_data$track == nm] <- plot_data$freq[plot_data$track == nm] * ratio
    }
  }

  if (is.null(track_colors)) {
    track_colors <- stats::setNames(scales::hue_pal()(length(track_names)), track_names)
  }

  dummy_data <- data.frame()
  if (!is.null(custom_y_limits)) {
    for (nm in names(custom_y_limits)) {
      dummy_data <- rbind(dummy_data, data.frame(
        bin_start = 0,
        freq = custom_y_limits[[nm]],
        track = nm,
        stringsAsFactors = FALSE
      ))
    }
    dummy_data$track <- factor(dummy_data$track, levels = track_names)
  }

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = bin_start, y = freq)) +
    ggplot2::geom_area(ggplot2::aes(fill = track), alpha = alpha) +
    ggplot2::scale_fill_manual(values = track_colors) +
    {if (nrow(dummy_data) > 0) ggplot2::geom_blank(data = dummy_data, ggplot2::aes(y = freq))} +
    ggplot2::facet_grid(track ~ ., scales = "free_y", switch = "y") +
    ggplot2::scale_x_continuous(
      labels = scales::label_number(scale = 1e-6, suffix = " Mb"),
      limits = c(0, genome_size),
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(expand = c(0.05, 0)) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      legend.position = "none",
      axis.line.x = ggplot2::element_line(color = "black", linewidth = 0.5),
      axis.line.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 16),
      strip.background = ggplot2::element_blank(),
      strip.text.y.left = ggplot2::element_text(angle = 0, face = "bold", size = 12, margin = ggplot2::margin(r = 10)),
      strip.placement = "outside",
      panel.spacing = grid::unit(0.5, "lines")
    ) +
    ggplot2::labs(
      title = if (is.null(title)) paste0(target_chr, " landscape (scaled)") else title,
      x = "Genomic position",
      y = "Density (scaled)"
    )

  if (!is.null(save_path)) {
    ggplot2::ggsave(save_path, plot = p, width = width, height = height, dpi = dpi)
  }

  p
}

