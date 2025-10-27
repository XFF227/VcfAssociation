#' Manhattan plot (toy)
#' @param res data.frame(CHR, POS, P)
#' @export
plot_manhattan <- function(res){
  res$CHR <- as.factor(res$CHR)
  ggplot2::ggplot(res, ggplot2::aes(x=POS, y=-log10(as.numeric(P)), col=CHR)) +
    ggplot2::geom_point(size=1) + ggplot2::theme_bw() +
    ggplot2::labs(x="Position", y="-log10(P)", title="Toy Manhattan")
}

#' QQ plot (toy)
#' @export
plot_qq <- function(res){
  p <- sort(as.numeric(res$P)); exp <- stats::ppoints(length(p))
  df <- data.frame(exp=-log10(exp), obs=-log10(p))
  ggplot2::ggplot(df, ggplot2::aes(exp, obs)) + ggplot2::geom_point() +
    ggplot2::geom_abline(slope=1, intercept=0, linetype=2) +
    ggplot2::theme_bw() + ggplot2::labs(x="Expected -log10 P", y="Observed -log10 P")
}

#' Forest plot for one gene (toy)
#' @export
plot_forest <- function(tbl){
  ggplot2::ggplot(tbl, ggplot2::aes(x=term, y=beta, ymin=beta-1.96*se, ymax=beta+1.96*se))+
    ggplot2::geom_pointrange() + ggplot2::coord_flip() + ggplot2::theme_bw() +
    ggplot2::labs(x="", y="Effect (log-odds)")
}
