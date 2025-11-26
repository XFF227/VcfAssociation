#' gwas_single
#' @title Perform single-variant association tests
#'
#' @description
#' This function runs the core GWAS step in the pipeline. Using the merged
#' phenotype–genotype table from read_phenotypes(), it performs a regression
#' test for each variant. The phenotype is modeled as a function of genotype
#' dosage and optional covariates.
#'
#' Binary phenotypes use logistic regression; continuous phenotypes use linear
#' regression. The output includes per-variant effect size, standard error, and
#' p-value, suitable for visualization with manhattan_plot().
#'
#' @param ph Merged data frame or list from read_phenotypes().
#' @param pheno_col Name of the phenotype column.
#' @param covars Optional character vector of covariate names.
#' @param id_col Sample ID column name.
#'
#' @return A data frame with CHROM, POS, beta, standard error, and p-value.
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
#' @seealso manhattan_plot
#' @references
#' **stats**: R Core Team. (2025). *R: A Language and Environment for Statistical Computing.*
#'  R Foundation for Statistical Computing, Vienna, Austria.
#'  <https://www.R-project.org/>
#'
#' **dplyr**: Wickham, H., François, R., Henry, L., & Müller, K. (2023).
#'  *dplyr: A Grammar of Data Manipulation.* R package version 1.x.
#'  <https://CRAN.R-project.org/package=dplyr>
#' @export
gwas_single <- function(ph,
                        pheno_col = "phenotype",
                        covars = NULL,
                        id_col = "sample") {
  # Standardize input to a merged data.frame
  merged <- if (is.list(ph) && !is.null(ph$merged)) ph$merged else ph
  if (!is.data.frame(merged)) stop("`ph` must be either a data.frame or a list with $merged from read_phenotypes().")
  
  # Resolve the sample column flexibly: prefer "sample", else use id_col
  sample_col <- if ("sample" %in% names(merged)) "sample" else if (id_col %in% names(merged)) id_col else NA_character_
  if (is.na(sample_col)) {
    stop("Missing sample ID column: neither 'sample' nor '", id_col, "' is present in the merged table.")
  }
  
  # Check required columns
  req <- c("CHROM", "POS", sample_col, "Dosage", pheno_col)
  miss <- setdiff(req, names(merged))
  if (length(miss)) stop("Missing columns in merged table: ", paste(miss, collapse = ", "))
  
  # Choose regression family (binary -> binomial; else gaussian)
  y <- merged[[pheno_col]]
  fam <- if (is.numeric(y) && all(y %in% c(0, 1), na.rm = TRUE)) stats::binomial() else stats::gaussian()
  
  # Build formula: phenotype ~ Dosage + covariates
  rhs <- c("Dosage", covars)
  fml <- stats::as.formula(paste(pheno_col, "~", paste(rhs, collapse = " + ")))
  
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please install.packages('dplyr')")
  if (!requireNamespace("purrr", quietly = TRUE)) stop("Please install.packages('purrr')")
  
  dplyr::group_by(merged, .data$CHROM, .data$POS) |>
    dplyr::group_modify(~{
      df <- .x
      # Skip if no phenotype info or no dosage variation
      if (all(is.na(df[[pheno_col]])) || length(unique(df$Dosage[!is.na(df$Dosage)])) < 2) {
        return(dplyr::tibble(beta = NA_real_, se = NA_real_, p = NA_real_))
      }
      fit <- try(stats::glm(fml, data = df, family = fam), silent = TRUE)
      if (inherits(fit, "try-error")) {
        return(dplyr::tibble(beta = NA_real_, se = NA_real_, p = NA_real_))
      }
      sm <- summary(fit)$coefficients
      if (!("Dosage" %in% rownames(sm))) {
        return(dplyr::tibble(beta = NA_real_, se = NA_real_, p = NA_real_))
      }
      dplyr::tibble(
        beta = unname(sm["Dosage", "Estimate"]),
        se   = unname(sm["Dosage", "Std. Error"]),
        p    = unname(sm["Dosage", if (fam$family == "binomial") "Pr(>|z|)" else "Pr(>|t|)"])
      )
    }) |>
    dplyr::ungroup() |>
    dplyr::arrange(.data$CHROM, .data$POS)
}
