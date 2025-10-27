#' gwas_single
#' @title Single-variant GWAS using GLM per site
#' @description Fit GLM per variant with optional covariates; return results and Manhattan plot.
#' @param genotypes Long genotype table from \code{read_vcf()}.
#' @param phenotypes Phenotype data including id and trait.
#' @param id_col Sample id column in \code{phenotypes}.
#' @param pheno_col Phenotype column.
#' @param covars Optional covariate column names.
#' @param family "gaussian" (quantitative) or "binomial" (binary).
#' @param p_adjust Method for \code{p.adjust}.
#' @return list(results=data.frame, manhattan=ggplot)
#' @export
#' @examples
#' \dontrun{
#' g <- v$genotypes
#' ph <- readr::read_csv("pheno.csv")
#' out <- gwas_single(g, ph, id_col="sample_id", pheno_col="trait", covars=c("age","sex"), family="binomial")
#' out$manhattan
#' dplyr::arrange(out$results, padj)[1:10,]
#' }
gwas_single <- function(genotypes, phenotypes, id_col, pheno_col,
                        covars = character(), family = c("gaussian","binomial"),
                        p_adjust = "BH") {
  family <- match.arg(family)
  stopifnot(is.data.frame(genotypes), is.data.frame(phenotypes))
  if (!id_col %in% names(phenotypes)) stop("id_col not in phenotypes")
  if (!pheno_col %in% names(phenotypes)) stop("pheno_col not in phenotypes")
  
  g_wide <- tidyr::pivot_wider(genotypes[, c("Sample","CHROM","POS","Dosage")],
                               names_from = c("CHROM","POS"),
                               values_from = "Dosage")
  colnames(g_wide)[1] <- id_col
  df <- dplyr::left_join(phenotypes, g_wide, by = id_col)
  
  var_cols <- setdiff(colnames(df), c(id_col, pheno_col, covars))
  var_cols <- var_cols[grepl("^\\w+_\\d+$", var_cols)]
  if (!length(var_cols)) stop("No variant dosage columns found.")
  
  fit_one <- function(vc) {
    ftxt <- paste0(pheno_col, " ~ ", vc,
                   if (length(covars)) paste0(" + ", paste(covars, collapse = " + ")) else "")
    fm <- stats::as.formula(ftxt)
    fit <- try(stats::glm(fm, data = df, family = if (family == "binomial") stats::binomial() else stats::gaussian()),
               silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    co <- summary(fit)$coef
    if (!vc %in% rownames(co)) return(NULL)
    est <- co[vc, "Estimate"]; se <- co[vc, "Std. Error"]
    p <- co[vc, if (family == "binomial") "Pr(>|z|)" else "Pr(>|t|)"]
    tibble::tibble(CHROM = sub("_.*$", "", vc),
                   POS   = as.integer(sub("^.*_", "", vc)),
                   beta  = est, se = se, p = p)
  }
  
  results <- dplyr::bind_rows(lapply(var_cols, fit_one))
  if (!nrow(results)) stop("No models converged.")
  results$padj <- p.adjust(results$p, method = p_adjust)
  results$neglog10p <- -log10(results$p)
  
  manh <- ggplot2::ggplot(results, ggplot2::aes(x = POS, y = neglog10p, color = CHROM)) +
    ggplot2::geom_point(size = 0.6) +
    ggplot2::facet_wrap(~CHROM, scales = "free_x") +
    ggplot2::labs(x = "Position", y = "-log10(p)", title = "Single-variant GWAS") +
    ggplot2::theme_bw()
  
  list(results = as.data.frame(results), manhattan = manh)
}
