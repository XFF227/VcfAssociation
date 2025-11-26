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
#'
#' @return A ggplot object representing the Manhattan plot.
#'
#' @examples
#' vcf_path <- system.file("extdata", "toy.vcf",
#'                         package = "VcfAssociation", mustWork = TRUE)
#' vcf <- read_vcf(vcf_path)
#' ph <- generate_phenotype(vcf_path,
#'                          chrom = "chr12", pos = 11161,
#'                          model = "carrier")
#' tmp <- tempfile(fileext = ".csv")
#' write.csv(ph, tmp, row.names = FALSE)
#' ph_list <- read_phenotypes(tmp, id_col = "sample", genotypes = vcf$genotypes)
#' gwas_res <- gwas_single(ph_list, pheno_col = "phenotype")
#' manhattan_plot(gwas_res)
#' @seealso gwas_single
#' @references
#' **ggplot2**: Wickham, H. (2016). *ggplot2: Elegant Graphics for Data Analysis.*
#'  Springer-Verlag New York. <https://ggplot2.tidyverse.org>
#'
#' **dplyr**: Wickham, H., François, R., Henry, L., & Müller, K. (2023).
#'  *dplyr: A Grammar of Data Manipulation.* R package version 1.x.
#'  <https://CRAN.R-project.org/package=dplyr>
#' @export
manhattan_plot <- function(df,
                           chr_col = "CHROM",
                           pos_col = "POS",
                           p_col   = "p",
                           genome_wide = 5e-8,
                           suggestive  = 1e-5,
                           highlight_snps = NULL,
                           snp_col = NULL,
                           annotate_top_n = 0) {
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please install.packages('dplyr')")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install.packages('ggplot2')")
  
  dplyr <- asNamespace("dplyr"); ggplot2 <- asNamespace("ggplot2")
  
  dat <- df
  # Basic column checks
  req <- c(chr_col, pos_col, p_col)
  miss <- setdiff(req, names(dat))
  if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))
  
  # Clean and standardize chromosome values (remove 'chr' prefix; map X/Y/MT)
  to_chr_label <- function(x) {
    x <- as.character(x)
    x <- sub("^chr", "", x, ignore.case = TRUE)
    x[x %in% c("X","x")]  <- "23"
    x[x %in% c("Y","y")]  <- "24"
    x[x %in% c("M","MT","m","mt")] <- "25"
    x
  }
  dat[[chr_col]] <- to_chr_label(dat[[chr_col]])
  
  # Keep only rows with finite p and valid positions
  dat <- dat[is.finite(dat[[p_col]]) & dat[[p_col]] > 0 & is.finite(dat[[pos_col]]), , drop = FALSE]
  if (!nrow(dat)) stop("No valid rows with finite positive p-values.")
  
  # Compute -log10(p); guard against p==0 (already filtered)
  dat$logp <- -log10(dat[[p_col]])
  
  # Order chromosomes naturally (numeric, then as factor for plotting)
  # Unknown/Non-numeric chromosomes will be placed after numeric ones
  suppressWarnings({
    chr_num <- as.integer(dat[[chr_col]])
  })
  # If NA for some rows, keep original label but push to end ordered by label
  dat$chr_num <- ifelse(is.na(chr_num), Inf, chr_num)
  
  # Build cumulative position across chromosomes
  # For non-numeric chroms (Inf), order by label to keep deterministic
  dat <- dplyr::as_tibble(dat) |>
    dplyr$group_by(.data[[chr_col]]) |>
    dplyr$mutate(chr_len = max(.data[[pos_col]], na.rm = TRUE)) |>
    dplyr$ungroup()
  
  # Chromosome order
  chr_order <- dat |>
    dplyr$distinct(.data[[chr_col]], chr_num) |>
    dplyr$arrange(.data$chr_num, .data[[chr_col]]) |>
    dplyr$pull(.data[[chr_col]])
  
  dat[[chr_col]] <- factor(dat[[chr_col]], levels = unique(chr_order), ordered = TRUE)
  
  # Compute cumulative offsets
  chr_table <- dat |>
    dplyr$group_by(.data[[chr_col]]) |>
    dplyr$summarise(chr_len = max(.data[[pos_col]], na.rm = TRUE), .groups = "drop") |>
    dplyr$mutate(offset = c(0, cumsum(chr_len)[-dplyr$n()]))
  
  dat <- dat |>
    dplyr$left_join(chr_table, by = setNames(chr_col, chr_col)) |>
    dplyr$mutate(pos_cum = .data[[pos_col]] + .data$offset)
  
  # X-axis tick positions at chromosome centers
  axis_df <- chr_table |>
    dplyr$mutate(center = offset + chr_len / 2)
  
  # Alternating colors by chromosome
  dat$chr_index <- as.integer(dat[[chr_col]])
  dat$color_group <- dat$chr_index %% 2
  
  p <- ggplot2$ggplot(dat, ggplot2$aes(x = .data$pos_cum, y = .data$logp, color = factor(.data$color_group))) +
    ggplot2$geom_point(size = 0.6, alpha = 0.8, na.rm = TRUE) +
    ggplot2$scale_color_manual(values = c("#4063D8", "#389826"), guide = "none") +
    ggplot2$scale_x_continuous(
      label = levels(dat[[chr_col]]),
      breaks = axis_df$center,
      expand = ggplot2$expansion(mult = c(0.01, 0.02))
    ) +
    ggplot2$labs(x = "Chromosome", y = expression(-log[10](p)), title = "Manhattan Plot") +
    ggplot2$theme_minimal(base_size = 12) +
    ggplot2$theme(
      panel.grid.major.x = ggplot2$element_blank(),
      panel.grid.minor.x = ggplot2$element_blank()
    )
  
  # Threshold lines
  if (!is.null(suggestive) && is.finite(suggestive) && suggestive > 0) {
    p <- p + ggplot2$geom_hline(yintercept = -log10(suggestive), linetype = "dashed")
  }
  if (!is.null(genome_wide) && is.finite(genome_wide) && genome_wide > 0) {
    p <- p + ggplot2$geom_hline(yintercept = -log10(genome_wide), color = "red")
  }
  
  # Highlight specific SNPs if requested
  if (!is.null(highlight_snps) && !is.null(snp_col) && snp_col %in% names(dat)) {
    dat$highlight <- dat[[snp_col]] %in% highlight_snps
    p <- p + ggplot2$geom_point(
      data = dat[dat$highlight, , drop = FALSE],
      ggplot2$aes(x = .data$pos_cum, y = .data$logp),
      size = 1.6, color = "orange", inherit.aes = FALSE
    )
  }
  
  # Optional annotation of top-N hits by p-value
  if (annotate_top_n > 0) {
    if (is.null(snp_col) || !(snp_col %in% names(dat))) {
      # fabricate simple labels if SNP IDs are not provided
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
      size = 3, vjust = -0.3, check_overlap = TRUE, inherit.aes = FALSE
    )
  }
  
  p
}
