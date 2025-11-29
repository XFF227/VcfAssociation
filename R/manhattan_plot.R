#' manhattan_plot
#' @title Create a Manhattan plot from GWAS results
#'
#' @description
#' This function generates a Manhattan plot from a GWAS result table, such as the
#' output of gwas_single(). It displays genomic position on the x-axis and
#' -log10(p-value) on the y-axis. Genome-wide and suggestive significance thresholds
#' can be visualized, and top hits can be labeled or highlighted.
#'
#' @param df Data frame of GWAS results.
#' @param chr_col Column name for chromosome (default "CHROM").
#' @param pos_col Column name for position (default "POS").
#' @param p_col Column name for p-values (default "p").
#' @param genome_wide Optional genome-wide significance threshold.
#' @param suggestive Optional suggestive threshold.
#' @param highlight_snps Optional vector of variant IDs to highlight.
#' @param snp_col Column name for variant IDs.
#' @param annotate_top_n Number of top hits to label.
#' @param base_size Base font size for the plot theme (default 16).
#' @param title_size Plot title font size (default 18).
#' @param axis_title_size Axis titles font size (default 15).
#' @param axis_text_size Axis text font size (default 12).
#' @param point_size Point size for all SNPs (default 1.2).
#' @param point_alpha Point transparency for all SNPs (default 0.9).
#' @param highlight_point_size Point size for highlighted SNPs (default 2.5).
#' @param suggestive_line_size Line width for suggestive threshold (default 0.8).
#' @param genome_line_size Line width for genome-wide threshold (default 1).
#' @param label_size Font size for annotated SNP labels (default 4).
#'
#' @import ggplot2
#' @importFrom dplyr as_tibble group_by mutate ungroup distinct arrange pull summarise n left_join
#' @importFrom stats setNames
#' @return A ggplot object representing the Manhattan plot.
#' @export
manhattan_plot <- function(df,
                           chr_col = "CHROM",
                           pos_col = "POS",
                           p_col   = "p",
                           genome_wide = 5e-8,
                           suggestive  = 1e-5,
                           highlight_snps = NULL,
                           snp_col = NULL,
                           annotate_top_n = 0,
                           base_size          = 16,
                           title_size         = 18,
                           axis_title_size    = 15,
                           axis_text_size     = 12,
                           point_size         = 1.2,
                           point_alpha        = 0.9,
                           highlight_point_size = 2.5,
                           suggestive_line_size = 0.8,
                           genome_line_size     = 1,
                           label_size           = 4) {
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please install.packages('dplyr')")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install.packages('ggplot2')")
  
  dplyr   <- asNamespace("dplyr")
  ggplot2 <- asNamespace("ggplot2")
  
  dat <- df
  req  <- c(chr_col, pos_col, p_col)
  miss <- setdiff(req, names(dat))
  if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))
  
  to_chr_label <- function(x) {
    x <- as.character(x)
    x <- sub("^chr", "", x, ignore.case = TRUE)
    x[x %in% c("X","x")]  <- "23"
    x[x %in% c("Y","y")]  <- "24"
    x[x %in% c("M","MT","m","mt")] <- "25"
    x
  }
  dat[[chr_col]] <- to_chr_label(dat[[chr_col]])
  
  dat <- dat[is.finite(dat[[p_col]]) & dat[[p_col]] > 0 & is.finite(dat[[pos_col]]), , drop = FALSE]
  if (!nrow(dat)) stop("No valid rows with finite positive p-values.")
  
  dat$logp <- -log10(dat[[p_col]])
  suppressWarnings({ chr_num <- as.integer(dat[[chr_col]]) })
  dat$chr_num <- ifelse(is.na(chr_num), Inf, chr_num)
  
  dat <- dplyr::as_tibble(dat) |>
    dplyr$group_by(.data[[chr_col]]) |>
    dplyr$mutate(chr_len = max(.data[[pos_col]], na.rm = TRUE)) |>
    dplyr$ungroup()
  
  chr_order <- dat |>
    dplyr$distinct(.data[[chr_col]], chr_num) |>
    dplyr$arrange(.data$chr_num, .data[[chr_col]]) |>
    dplyr$pull(.data[[chr_col]])
  
  dat[[chr_col]] <- factor(dat[[chr_col]], levels = unique(chr_order), ordered = TRUE)
  
  chr_table <- dat |>
    dplyr$group_by(.data[[chr_col]]) |>
    dplyr$summarise(chr_len = max(.data[[pos_col]], na.rm = TRUE), .groups = "drop") |>
    dplyr$mutate(offset = c(0, cumsum(chr_len)[-dplyr$n()]))
  
  dat <- dat |>
    dplyr$left_join(chr_table, by = setNames(chr_col, chr_col)) |>
    dplyr$mutate(pos_cum = .data[[pos_col]] + .data$offset)
  
  axis_df <- chr_table |>
    dplyr$mutate(center = offset + chr_len / 2)
  
  dat$chr_index   <- as.integer(dat[[chr_col]])
  dat$color_group <- dat$chr_index %% 2
  
  p <- ggplot2$ggplot(
    dat,
    ggplot2$aes(x = .data$pos_cum, y = .data$logp, color = factor(.data$color_group))
  ) +
    ggplot2$geom_point(size = point_size, alpha = point_alpha, na.rm = TRUE) +
    ggplot2$scale_color_manual(values = c("#4063D8", "#389826"), guide = "none") +
    ggplot2$scale_x_continuous(
      label  = levels(dat[[chr_col]]),
      breaks = axis_df$center,
      expand = ggplot2$expansion(mult = c(0.01, 0.02))
    ) +
    ggplot2$labs(
      x     = "Chromosome",
      y     = expression(-log[10](p)),
      title = "Manhattan Plot"
    ) +
    ggplot2$theme_minimal(base_size = base_size) +
    ggplot2$theme(
      plot.title      = ggplot2$element_text(hjust = 0.5, face = "bold", size = title_size),
      axis.title.x    = ggplot2$element_text(size = axis_title_size),
      axis.title.y    = ggplot2$element_text(size = axis_title_size),
      axis.text.x     = ggplot2$element_text(size = axis_text_size),
      axis.text.y     = ggplot2$element_text(size = axis_text_size),
      panel.grid.major.x = ggplot2$element_blank(),
      panel.grid.minor.x = ggplot2$element_blank()
    )
  
  if (!is.null(suggestive) && is.finite(suggestive) && suggestive > 0) {
    p <- p + ggplot2$geom_hline(
      yintercept = -log10(suggestive),
      linetype   = "dashed",
      size       = suggestive_line_size
    )
  }
  if (!is.null(genome_wide) && is.finite(genome_wide) && genome_wide > 0) {
    p <- p + ggplot2$geom_hline(
      yintercept = -log10(genome_wide),
      color      = "red",
      size       = genome_line_size
    )
  }
  
  if (!is.null(highlight_snps) && !is.null(snp_col) && snp_col %in% names(dat)) {
    dat$highlight <- dat[[snp_col]] %in% highlight_snps
    p <- p + ggplot2$geom_point(
      data = dat[dat$highlight, , drop = FALSE],
      ggplot2$aes(x = .data$pos_cum, y = .data$logp),
      size = highlight_point_size,
      color = "orange",
      inherit.aes = FALSE
    )
  }
  
  if (annotate_top_n > 0) {
    if (is.null(snp_col) || !(snp_col %in% names(dat))) {
      dat$.__tmp_label__ <- paste0(dat[[chr_col]], ":", dat[[pos_col]])
      lab_col <- ".__tmp_label__"
    } else {
      lab_col <- snp_col
    }
    top_idx <- order(dat[[p_col]], decreasing = FALSE)[seq_len(min(annotate_top_n, nrow(dat)))]
    ann <- dat[top_idx, , drop = FALSE]
    p <- p + ggplot2$geom_text(
      data = ann,
      ggplot2$aes(x = .data$pos_cum, y = .data$logp, label = .data[[lab_col]]),
      size = label_size,
      vjust = -0.4,
      check_overlap = TRUE,
      inherit.aes   = FALSE
    )
  }
  
  p
}
