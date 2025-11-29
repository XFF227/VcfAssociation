#' gwas_single
#'
#' @title Perform single-variant association tests
#'
#' @description
#' `gwas_single()` runs a standard single-variant association scan on a merged
#' phenotype–genotype table. Using the unified table produced by
#' [read_phenotypes()], it fits a separate regression model for each variant,
#' with the phenotype as the outcome and genotype dosage as the main predictor.
#' Optional covariates can be supplied and are included in every model.
#'
#' The regression family is chosen automatically from the phenotype column:
#'
#' * if the phenotype is numeric and takes only the values `0` and `1`
#'   (ignoring `NA`), or a two-level factor, a logistic regression is fitted
#'   with `stats::binomial()`;
#' * otherwise a linear model is fitted with `stats::gaussian()`.
#'
#' For each variant, the function extracts the coefficient of the dosage term
#' (`beta`), its standard error (`se`), and the corresponding Wald p-value
#' (`p`). The result is returned as a tidy data frame containing one row per
#' variant, ordered by chromosome and position, and is suitable for downstream
#' visualization with [manhattan_plot()].
#'
#' @param ph A merged data frame or a list as returned by [read_phenotypes()].
#'   If a list is supplied, the merged phenotype–genotype table is taken from
#'   the `$merged` element.
#' @param pheno_col Character scalar giving the name of the phenotype column.
#'   Default is `"phenotype"`.
#' @param covars Optional character vector of covariate column names. These
#'   columns must all be present in the merged table. If `NULL` (the default),
#'   only the dosage term is included in the model.
#' @param id_col Character scalar giving the sample ID column name. This is
#'   only used to validate that an ID column is present; the column itself is
#'   not included in the regression.
#'   
#' @importFrom dplyr group_by group_modify ungroup arrange tibble
#' @importFrom stats as.formula glm binomial gaussian
#' @return
#' A tibble with one row per variant and the following columns:
#'
#' \describe{
#'   \item{CHROM}{Chromosome.}
#'   \item{POS}{Base-pair position.}
#'   \item{REF}{Reference allele.}
#'   \item{ALT}{Alternate allele.}
#'   \item{beta}{Estimated regression coefficient for genotype dosage.}
#'   \item{se}{Standard error of `beta`.}
#'   \item{p}{Two-sided Wald p-value for `beta`.}
#' }
#'
#' Variants for which the model cannot be fitted (for example because all
#' phenotypes are missing, there is no variation in dosage, or the fit fails)
#' receive `NA` in the `beta`, `se`, and `p` columns.
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
#' head(gwas_res)
#'
#' @seealso [manhattan_plot()]
#'
#' @references
#' \strong{stats}: R Core Team. (2025).
#' *R: A Language and Environment for Statistical Computing.*
#' R Foundation for Statistical Computing, Vienna, Austria.
#' <https://www.R-project.org/>
#'
#' \strong{dplyr}: Wickham, H., François, R., Henry, L., & Müller, K. (2023).
#' *dplyr: A Grammar of Data Manipulation.* R package version 1.x.
#' <https://CRAN.R-project.org/package=dplyr>
#'
#' @export
gwas_single <- function(ph,
                        pheno_col = "phenotype",
                        covars = NULL,
                        id_col = "sample") {
  merged <- .gwas_standardize_input(ph, id_col = id_col)
  
  if (!pheno_col %in% names(merged)) {
    stop("`pheno_col` not found in merged table.")
  }
  
  if (!is.null(covars)) {
    missing_cov <- setdiff(covars, names(merged))
    if (length(missing_cov) > 0) {
      stop("Missing covariate(s) in merged table: ",
           paste(missing_cov, collapse = ", "))
    }
  }
  
  if (!"Dosage" %in% names(merged)) {
    stop("Merged table must contain a 'Dosage' column produced by read_vcf().")
  }
  
  fam <- .gwas_choose_family(merged[[pheno_col]])
  
  fml <- .gwas_build_formula(
    outcome   = pheno_col,
    dosage    = "Dosage",
    covariate = covars
  )
  
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Please install.packages(\"dplyr\") to use gwas_single().")
  }
  
  dplyr::group_by(merged, .data$CHROM, .data$POS, .data$REF, .data$ALT) |>
    dplyr::group_modify(~ .gwas_fit_one_variant(
      dat       = .x,
      formula   = fml,
      family    = fam,
      pheno_col = pheno_col,
      dosage    = "Dosage"
    )) |>
    dplyr::ungroup() |>
    dplyr::arrange(.data$CHROM, .data$POS)
}

# ---- internal helpers for gwas_single ---------------------------------------

#' @keywords internal
.gwas_standardize_input <- function(ph, id_col = "sample") {
  merged <- if (is.list(ph) && !is.null(ph$merged)) ph$merged else ph
  if (!is.data.frame(merged)) {
    stop("`ph` must be either a data.frame or a list with $merged from read_phenotypes().")
  }
  
  sample_col <- if ("sample" %in% names(merged)) {
    "sample"
  } else if (id_col %in% names(merged)) {
    id_col
  } else {
    NA_character_
  }
  
  if (is.na(sample_col)) {
    stop("Missing sample ID column: neither 'sample' nor '", id_col,
         "' is present in the merged table.")
  }
  
  merged
}

#' @keywords internal
.gwas_choose_family <- function(y) {
  if (is.factor(y) && nlevels(y) == 2L) {
    stats::binomial()
  } else if (is.numeric(y) && all(y %in% c(0, 1), na.rm = TRUE)) {
    stats::binomial()
  } else {
    stats::gaussian()
  }
}

#' @keywords internal
.gwas_build_formula <- function(outcome, dosage = "Dosage", covariate = NULL) {
  rhs <- c(dosage, covariate)
  rhs <- rhs[!is.na(rhs) & nzchar(rhs)]
  if (length(rhs) == 0L) {
    stop("Formula must contain at least a dosage term.")
  }
  stats::as.formula(paste(outcome, "~", paste(rhs, collapse = " + ")))
}

#' @keywords internal
.gwas_fit_one_variant <- function(dat,
                                  formula,
                                  family,
                                  pheno_col,
                                  dosage = "Dosage") {
  if (all(is.na(dat[[pheno_col]])) ||
      length(unique(dat[[dosage]][!is.na(dat[[dosage]])])) < 2L) {
    return(dplyr::tibble(beta = NA_real_, se = NA_real_, p = NA_real_))
  }
  
  fit <- try(stats::glm(formula, data = dat, family = family), silent = TRUE)
  if (inherits(fit, "try-error")) {
    return(dplyr::tibble(beta = NA_real_, se = NA_real_, p = NA_real_))
  }
  
  sm <- summary(fit)$coefficients
  if (!(dosage %in% rownames(sm))) {
    return(dplyr::tibble(beta = NA_real_, se = NA_real_, p = NA_real_))
  }
  
  stat_col <- if (family$family == "binomial") "Pr(>|z|)" else "Pr(>|t|)"
  
  dplyr::tibble(
    beta = unname(sm[dosage, "Estimate"]),
    se   = unname(sm[dosage, "Std. Error"]),
    p    = unname(sm[dosage, stat_col])
  )
}