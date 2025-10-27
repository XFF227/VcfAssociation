#' assoc_gene
#' @title Associate a variant with a gene-level binary label
#' @description Model: Gene(+/-) ~ Variant + covariates using Logistic/Firth Logistic.
#' @param df Data with gene_flag (0/1), dosage and covariates.
#' @param gene_col Column name of gene flag (0/1).
#' @param dosage Dosage column name.
#' @param covars Optional covariate names.
#' @param model "logistic" or "firth".
#' @return list(tidy, model)
#' @export
#' @examples
#' \dontrun{
#' assoc_gene(df, gene_col="is_driver", dosage="Dosage", covars=c("age","sex"), model="logistic")
#' }
assoc_gene <- function(df, gene_col = "gene_flag", dosage = "Dosage",
                       covars = character(), model = c("logistic","firth")) {
  model <- match.arg(model)
  if (!all(c(gene_col, dosage) %in% names(df))) stop("gene_col or dosage missing.")
  fm <- stats::as.formula(paste0(gene_col, " ~ ", paste(c(dosage, covars), collapse = " + ")))
  
  if (model == "logistic") {
    fit <- stats::glm(fm, data = df, family = stats::binomial())
    co  <- summary(fit)$coef
    est <- co[dosage, "Estimate"]; se <- co[dosage, "Std. Error"]
    p   <- co[dosage, "Pr(>|z|)"]
    ci  <- suppressMessages(stats::confint(fit, parm = dosage))
    tidy <- data.frame(term=dosage, beta=est, se=se, p=p,
                       OR=exp(est), OR_lo=exp(ci[1]), OR_hi=exp(ci[2]), model="logistic")
    return(list(tidy = tidy, model = fit))
  }
  
  if (!requireNamespace("logistf", quietly = TRUE)) stop("Install 'logistf' for Firth logistic.")
  fit <- logistf::logistf(fm, data = df)
  est <- fit$coef[dosage]; se <- sqrt(diag(vcov(fit)))[dosage]
  ci  <- fit$ci[dosage, c("lower","upper")]
  p   <- fit$prob[dosage]
  tidy <- data.frame(term=dosage, beta=est, se=se, p=p,
                     OR=exp(est), OR_lo=ci[1], OR_hi=ci[2], model="firth")
  list(tidy = tidy, model = fit)
}
