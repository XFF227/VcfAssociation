#' Simple gene-level logistic regression (toy)
#' @param var_df annotated variant table
#' @param pheno data.frame(SAMPLE, PHENO[, AGE, SEX])
#' @param gene character
#' @return data.frame(term, beta, se, p)
#' @export
assoc_gene <- function(var_df, pheno, gene){
  x <- subset(var_df, GENE == gene)
  carr <- aggregate((GT != "0/0") ~ SAMPLE, x, mean)
  colnames(carr) <- c("SAMPLE","CARR")
  d <- merge(pheno, carr, by="SAMPLE", all.x=TRUE); d$CARR[is.na(d$CARR)] <- 0
  fit <- stats::glm(PHENO ~ CARR, data=d, family=stats::binomial())
  s <- summary(fit)$coef["CARR", , drop=FALSE]
  data.frame(term="CARR", beta=s[1], se=s[2], p=s[4], row.names=NULL)
}
