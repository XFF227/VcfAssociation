#' gwas_single
#' @title Single-variant regression using merged pheno+geno
#' @description Accepts either a data.frame that's already merged OR the output of read_phenotypes().
#' @param ph Either a data.frame (already merged) or the list returned by read_phenotypes().
#' @param pheno_col Name of phenotype column in merged table.
#' @param covars Optional character vector of covariate names in merged table.
#' @param id_col Sample ID column name (should match what you used in read_phenotypes()).
#' @return A data.frame with CHROM, POS, beta, se, and p values for the Dosage effect.
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
