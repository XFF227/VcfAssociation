#' assoc_single
#' @title Association for one variant and one phenotype
#' @description Fit Linear/Logistic or Firth Logistic (if \pkg{logistf} installed).
#' @param df Analysis data containing outcome, dosage, optional covars.
#' @param outcome Outcome column name.
#' @param dosage Dosage column name (0/1/2).
#' @param covars Optional covariates.
#' @param model "logistic","linear","firth".
#' @return list(tidy=data.frame, model=fit)
#' @export
#' @examples
#' \dontrun{
#' dat <- merge(pheno_df, subset(v$genotypes, CHROM=="1" & POS==1234)[,c("Sample","Dosage")],
#'              by.x="sample_id", by.y="Sample")
#' assoc_single(dat, outcome="trait", dosage="Dosage", covars=c("age","sex"), model="logistic")
#' }
assoc_single <- function(df, outcome, dosage = "Dosage", covars = character(),
                         model = c("logistic","linear","firth")) {
  model <- match.arg(model)
  if (!all(c(outcome, dosage) %in% names(df))) stop("Outcome or dosage missing in df.")
  
  rhs <- paste(c(dosage, covars), collapse = " + ")
  fm  <- stats::as.formula(paste0(outcome, " ~ ", rhs))
  
  if (model == "linear") {
    fit <- stats::lm(fm, data = df)
    co  <- summary(fit)$coef
    est <- co[dosage, "Estimate"]; se <- co[dosage, "Std. Error"]
    p   <- co[dosage, "Pr(>|t|)"]
    ci  <- stats::confint(fit, parm = dosage)
    tidy <- data.frame(term=dosage, beta=est, se=se, p=p, ci_lo=ci[1], ci_hi=ci[2], model="linear")
    return(list(tidy = tidy, model = fit))
  }
  
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
