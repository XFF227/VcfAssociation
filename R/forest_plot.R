#' forest_plot
#' @title Create a forest plot of effect sizes
#'
#' @description
#' This function visualizes model results from build_model() as a forest plot.
#' Each row represents a variant or region, showing its estimated effect size
#' and confidence interval. It helps interpret the strength and direction
#' of associations across multiple loci.
#'
#' @param res Data frame returned by build_model().
#' @param label_col Column used for y-axis labels.
#' @param facet_col Optional column to group results.
#' @param show_errors Whether to display errors for failed fits.
#' @param annotate_p Whether to label p-values next to estimates.
#' @param base_size Base font size for the plot theme (default 16).
#' @param title_size Plot title font size (default 18).
#' @param axis_title_size Axis titles font size (default 14).
#' @param axis_text_size Axis text font size (default 12).
#' @param point_size Point size for valid estimates (default 3).
#' @param errorbar_size Line width for confidence intervals (default 0.8).
#' @param zero_line_size Line width for the vertical reference line at 0 (default 0.9).
#' @param error_point_size Point size for invalid rows (default 3).
#' @param error_label_size Text size for error messages (default 3.6).
#' @param p_label_size Text size for p-value labels (default 3.6).
#' @import ggplot2
#' @importFrom dplyr as_tibble tibble group_by mutate ungroup
#' @return A ggplot object showing the forest plot.
#' @seealso build_model
#' @export
forest_plot <- function(res,
                        label_col        = "table_id",
                        facet_col        = "outcome",
                        show_errors      = TRUE,
                        annotate_p       = TRUE,
                        base_size        = 16,
                        title_size       = 18,
                        axis_title_size  = 14,
                        axis_text_size   = 12,
                        point_size       = 3,
                        errorbar_size    = 0.8,
                        zero_line_size   = 0.9,
                        error_point_size = 3,
                        error_label_size = 3.6,
                        p_label_size     = 3.6) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install.packages('ggplot2')")
  if (!requireNamespace("dplyr", quietly = TRUE))   stop("Please install.packages('dplyr')")
  if (!is.data.frame(res)) stop("`res` must be a data.frame produced by build_model().")
  
  dplyr   <- asNamespace("dplyr")
  ggplot2 <- asNamespace("ggplot2")
  
  df <- dplyr::as_tibble(res)
  
  ## required columns
  need <- c("model", "beta", "p", label_col)
  miss <- setdiff(need, names(df))
  if (length(miss)) stop("Missing columns in `res`: ", paste(miss, collapse = ", "))
  
  ## logistic-like models
  is_logit <- df$model %in% c("logistic", "firth")
  
  if (!"OR" %in% names(df)) df$OR <- NA_real_
  df$OR[is_logit & is.na(df$OR) & is.finite(df$beta)] <- exp(df$beta)
  
  get_col <- function(nm) if (nm %in% names(df)) df[[nm]] else rep(NA_real_, nrow(df))
  lin_ci_lo <- get_col("ci_lo"); lin_ci_hi <- get_col("ci_hi")
  log_or_lo <- get_col("OR_lo"); log_or_hi <- get_col("OR_hi")
  
  effect <- df$beta
  lo <- lin_ci_lo
  hi <- lin_ci_hi
  
  effect[is_logit] <- df$beta[is_logit]  # log(OR)
  lo[is_logit]     <- ifelse(is.na(log_or_lo[is_logit]),
                             NA_real_,
                             log(log_or_lo[is_logit]))
  hi[is_logit]     <- ifelse(is.na(log_or_hi[is_logit]),
                             NA_real_,
                             log(log_or_hi[is_logit]))
  
  has_error_col  <- "error" %in% names(df)
  invalid_reason <- if (has_error_col) df$error else NA_character_
  valid          <- is.finite(effect)
  
  ord_key <- ifelse(is.finite(df$p), df$p, Inf)
  ylab    <- df[[label_col]]
  
  df_plot <- dplyr::tibble(
    label  = ylab,
    effect = effect,
    lo     = lo,
    hi     = hi,
    p      = df$p,
    model  = df$model,
    facet  = if (!is.null(facet_col) && facet_col %in% names(df)) df[[facet_col]] else factor(""),
    valid  = valid,
    reason = invalid_reason
  )
  
  df_plot <- df_plot |>
    dplyr::group_by(.data$facet) |>
    dplyr::mutate(label = factor(.data$label,
                                 levels = .data$label[order(ord_key)],
                                 ordered = TRUE)) |>
    dplyr::ungroup()
  
  good <- df_plot[df_plot$valid, , drop = FALSE]
  bad  <- df_plot[!df_plot$valid, , drop = FALSE]
  
  ## base plot
  p <- ggplot2$ggplot() +
    ggplot2$geom_vline(xintercept = 0,
                       linetype   = "dashed",
                       color      = "#F16916",
                       linewidth  = zero_line_size) +
    ggplot2$theme_minimal(base_size = base_size) +
    ggplot2$theme(
      plot.title   = ggplot2$element_text(hjust = 0.5,
                                          face  = "bold",
                                          size  = title_size),
      axis.title.x = ggplot2$element_text(size = axis_title_size),
      axis.title.y = ggplot2$element_blank(),
      axis.text.x  = ggplot2$element_text(size = axis_text_size),
      axis.text.y  = ggplot2$element_text(size = axis_text_size)
    ) +
    ggplot2$labs(
      x     = "Effect (beta or log(OR))",
      y     = NULL,
      title = "Forest plot"
    )
  
  ## valid rows
  if (nrow(good) > 0) {
    p <- p +
      ggplot2$geom_errorbar(
        data        = good,
        ggplot2$aes(y = .data$label,
                    xmin = .data$lo,
                    xmax = .data$hi),
        width       = 0.2,
        alpha       = 0.7,
        na.rm       = TRUE,
        orientation = "y",
        color       = "#1C65A3",
        linewidth   = errorbar_size
      ) +
      ggplot2$geom_point(
        data  = good,
        ggplot2$aes(y = .data$label, x = .data$effect),
        size  = point_size,
        alpha = 0.9,
        color = "#1C65A3"
      )
  }
  
  ## invalid rows
  if (nrow(bad) > 0) {
    p <- p +
      ggplot2$geom_point(
        data  = bad,
        ggplot2$aes(y = .data$label, x = 0),
        size  = error_point_size,
        shape = 1,
        color = "grey30"
      )
    
    if (show_errors) {
      short_reason <- bad$reason
      short_reason[is.na(short_reason) | short_reason == ""] <- "invalid / no information"
      short_reason <- ifelse(
        nchar(short_reason) > 28,
        paste0(substr(short_reason, 1, 25), "..."),
        short_reason
      )
      p <- p +
        ggplot2$geom_text(
          data = transform(bad, reason_short = short_reason),
          ggplot2$aes(y = .data$label, x = 0, label = .data$reason_short),
          nudge_x = 0.02,
          hjust   = 0,
          size    = error_label_size,
          color   = "grey30"
        )
    }
  }
  
  ## p-value annotation
  if (annotate_p && nrow(df_plot) > 0) {
    p_lab <- ifelse(
      is.finite(df_plot$p),
      formatC(df_plot$p, format = "e", digits = 2),
      "NA"
    )
    p <- p +
      ggplot2$geom_text(
        data  = transform(df_plot, p_lab = paste0("p=", p_lab)),
        ggplot2$aes(y = .data$label, x = Inf, label = .data$p_lab),
        hjust = 1.05,
        size  = p_label_size,
        color = "grey20"
      ) +
      ggplot2$coord_cartesian(clip = "off") +
      ggplot2$theme(plot.margin = ggplot2$margin(5.5, 40, 5.5, 5.5))
  }
  
  if (!is.null(facet_col) && facet_col %in% names(res)) {
    p <- p + ggplot2$facet_wrap(~ facet, scales = "free_y", ncol = 1)
  }
  
  p
}
