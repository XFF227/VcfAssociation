#' assoc_viz
#' @title Forest plot for association results
#' @description Plot either odds ratios (logistic models) or beta with CI (linear).
#' @param results Data frame with columns:
#'   For OR scale: label, OR, OR_lo, OR_hi
#'   For beta scale: label, beta, ci_lo, ci_hi
#' @param effect_scale "OR" or "beta".
#' @param sort_by Optional column name to order labels.
#' @return ggplot object
#' @export
#' @examples
#' \dontrun{
#' p1 <- assoc_viz(res_or, effect_scale="OR", sort_by="OR")
#' p2 <- assoc_viz(res_b, effect_scale="beta", sort_by="beta")
#' }
assoc_viz <- function(results, effect_scale = c("OR","beta"), sort_by = NULL) {
  effect_scale <- match.arg(effect_scale)
  df <- results
  
  if (!is.null(sort_by) && sort_by %in% names(df)) {
    df <- df[order(df[[sort_by]]), , drop = FALSE]
  }
  
  if (effect_scale == "OR") {
    need <- c("label","OR","OR_lo","OR_hi")
    if (!all(need %in% names(df))) stop("Need columns: label, OR, OR_lo, OR_hi.")
    p <- ggplot2::ggplot(df, ggplot2::aes(y = factor(label, levels = df$label), x = OR)) +
      ggplot2::geom_point() +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = OR_lo, xmax = OR_hi), height = 0.2) +
      ggplot2::geom_vline(xintercept = 1, linetype = 2) +
      ggplot2::labs(x = "Odds Ratio (95% CI)", y = NULL, title = "Association Forest Plot") +
      ggplot2::theme_bw()
    return(p)
  }
  
  need <- c("label","beta","ci_lo","ci_hi")
  if (!all(need %in% names(df))) stop("Need columns: label, beta, ci_lo, ci_hi.")
  ggplot2::ggplot(df, ggplot2::aes(y = factor(label, levels = df$label), x = beta)) +
    ggplot2::geom_point() +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = ci_lo, xmax = ci_hi), height = 0.2) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::labs(x = "Effect (95% CI)", y = NULL, title = "Association Forest Plot") +
    ggplot2::theme_bw()
}
