#' forest_plot
#' @title Forest plot for build_model() results (robust to NA / no-variation cases)
#' @description
#' Draw a forest plot from the tidy summary returned by \code{build_model()}.
#' Handles both linear and (Firth) logistic models in one figure by plotting
#' **effects on a common scale**:
#' - Linear models: use \code{beta} and its CI.
#' - Logistic/Firth: use \code{log(OR)} and \code{log(CI)} so that 0 is the
#'   null line for all panels.
#'
#' Rows where the dosage term is missing (e.g., no genotype variation) or a fit
#' failed will be shown as hollow markers with a short error label (optional).
#'
#' @param res Data frame returned by \code{build_model()} (one row per result).
#' @param label_col Column used for y-axis labels (default: "table_id").
#' @param facet_col Optional column to facet by (e.g., "outcome" or "model"). Use \code{NULL} for no facet.
#' @param show_errors Logical; whether to display short error text for invalid rows. Default TRUE.
#' @param annotate_p Logical; whether to print p-values at the right side. Default TRUE.
#' @return A \code{ggplot} object.
#' @export
forest_plot <- function(res,
                               label_col = "table_id",
                               facet_col = "outcome",
                               show_errors = TRUE,
                               annotate_p  = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install.packages('ggplot2')")
  if (!requireNamespace("dplyr", quietly = TRUE))   stop("Please install.packages('dplyr')")
  if (!is.data.frame(res)) stop("`res` must be a data.frame produced by build_model().")
  
  dplyr <- asNamespace("dplyr"); ggplot2 <- asNamespace("ggplot2")
  
  df <- dplyr::as_tibble(res)
  
  # Ensure required columns exist (tolerate missing CI columns)
  need <- c("model", "beta", "p", label_col)
  miss <- setdiff(need, names(df))
  if (length(miss)) stop("Missing columns in `res`: ", paste(miss, collapse = ", "))
  
  # Identify logistic-like models
  is_logit <- df$model %in% c("logistic", "firth")
  
  # If OR not provided but beta is, derive OR for display
  if (!"OR" %in% names(df)) df$OR <- NA_real_
  df$OR[is_logit & is.na(df$OR) & is.finite(df$beta)] <- exp(df$beta)
  
  # Build a unified "effect" scale: beta for linear; log(OR) for logistic/firth
  # Also build corresponding CI columns on the same scale.
  # If CI is missing, keep NA (plot will show points without errorbars).
  get_col <- function(nm) if (nm %in% names(df)) df[[nm]] else rep(NA_real_, nrow(df))
  
  lin_ci_lo <- get_col("ci_lo"); lin_ci_hi <- get_col("ci_hi")
  log_or_lo <- get_col("OR_lo"); log_or_hi <- get_col("OR_hi")
  
  effect <- df$beta
  lo <- lin_ci_lo
  hi <- lin_ci_hi
  effect[is_logit] <- df$beta[is_logit]                        # log(OR) == beta
  lo[is_logit]     <- ifelse(is.na(log_or_lo[is_logit]), NA_real_, log(log_or_lo[is_logit]))
  hi[is_logit]     <- ifelse(is.na(log_or_hi[is_logit]), NA_real_, log(log_or_hi[is_logit]))
  
  # Valid row: effect is finite; Invalid row: NA/Inf or explicit error message
  has_error_col <- "error" %in% names(df)
  invalid_reason <- if (has_error_col) df$error else NA_character_
  valid <- is.finite(effect)
  
  # Ordering: by p ascending (NA last)
  ord_key <- ifelse(is.finite(df$p), df$p, Inf)
  # y as factor with chosen label
  ylab <- df[[label_col]]
  df_plot <- dplyr::tibble(
    label = ylab,
    effect = effect,
    lo = lo,
    hi = hi,
    p = df$p,
    model = df$model,
    facet = if (!is.null(facet_col) && facet_col %in% names(df)) df[[facet_col]] else factor(""),
    valid = valid,
    reason = invalid_reason
  )
  
  # Reorder labels within each facet by p-value (smallest on top)
  df_plot <- df_plot |>
    dplyr$group_by(.data$facet) |>
    dplyr$mutate(label = factor(.data$label, levels = .data$label[order(ord_key)], ordered = TRUE)) |>
    dplyr$ungroup()
  
  # Split valid / invalid for separate geoms
  good <- df_plot[df_plot$valid, , drop = FALSE]
  bad  <- df_plot[!df_plot$valid, , drop = FALSE]
  
  # Compose base plot (common 0 reference line: beta=0 or log(OR)=0)
  p <- ggplot2$ggplot() +
    ggplot2$geom_vline(xintercept = 0, linetype = "dashed", color = "#F16916") +
    ggplot2$theme_minimal(base_size = 12) +
    ggplot2$labs(x = "Effect (beta or log(OR))", y = NULL, title = "Forest plot")
  
  # Valid rows: point + CI
  if (nrow(good) > 0) {
    p <- p +
      ggplot2::geom_errorbar(
        color = "#1C65A3",
        data = good,
        ggplot2::aes(y = .data$label, xmin = .data$lo, xmax = .data$hi),
        width = 0.2,                
        alpha = 0.7,
        na.rm = TRUE,
        orientation = "y"          
      ) +
      ggplot2::geom_point(
        color = "#1C65A3",
        data = good,
        ggplot2::aes(y = .data$label, x = .data$effect),
        size = 2.2,
        alpha = 0.9
      )
  }
  # Invalid rows: hollow point at 0 with optional short error text
  if (nrow(bad) > 0) {
    p <- p +
      ggplot2$geom_point(
        data = bad,
        ggplot2$aes(y = .data$label, x = 0),
        size = 2.2, shape = 1, color = "grey30"
      )
    if (show_errors) {
      # shorten long messages
      short_reason <- bad$reason
      short_reason[is.na(short_reason) | short_reason == ""] <- "invalid / no information"
      short_reason <- ifelse(nchar(short_reason) > 28,
                             paste0(substr(short_reason, 1, 25), "..."),
                             short_reason)
      p <- p + ggplot2$geom_text(
        data = transform(bad, reason_short = short_reason),
        ggplot2$aes(y = .data$label, x = 0, label = .data$reason_short),
        nudge_x = 0.02, hjust = 0, size = 3.2, color = "grey30"
      )
    }
  }
  
  # Optional p-value annotation (on the right margin)
  if (annotate_p && nrow(df_plot) > 0) {
    p_lab <- ifelse(is.finite(df_plot$p),
                    formatC(df_plot$p, format = "e", digits = 2),
                    "NA")
    p <- p + ggplot2$geom_text(
      data = transform(df_plot, p_lab = paste0("p=", p_lab)),
      ggplot2$aes(y = .data$label, x = Inf, label = .data$p_lab),
      hjust = 1.05, size = 3.2, color = "grey20"
    ) +
      ggplot2$coord_cartesian(clip = "off") +
      ggplot2$theme(plot.margin = ggplot2$margin(5.5, 40, 5.5, 5.5))
  }
  
  # Facet if requested
  if (!is.null(facet_col) && facet_col %in% names(res)) {
    p <- p + ggplot2$facet_wrap(~ facet, scales = "free_y", ncol = 1)
  }
  
  p
}
