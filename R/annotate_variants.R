#' Minimal annotator (toy)
#' @param df output of read_vcf()
#' @return df with dummy gene/impact columns (for vignette/tests)
#' @export
annotate_variants <- function(df){
  df$GENE <- paste0("GENE", ((df$POS %% 3)+1))
  df$IMPACT <- c("LOW","MODERATE","HIGH")[(df$POS %% 3)+1]
  df
}
