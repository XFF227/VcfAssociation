#' build_model
#' @title Fit gene- or region-level association models
#'
#' @description
#' This function performs multi-variant modeling using data prepared by
#' prepare_list(). It fits linear or logistic (Firth) regression models that link
#' phenotype to one or more genotype dosages and optional covariates.
#' It is used to identify the joint effect of variants within a locus or gene.
#'
#' @param df_list List of data frames prepared by prepare_list().
#' @param outcomes Phenotype column names.
#' @param dosage Name of dosage column.
#' @param covars Optional covariate names.
#' @param model Either "logistic" or "linear".
#'
#' @return A data frame summarizing coefficients, standard errors, and p-values.
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
#' variants <- data.frame(CHROM = "chr12", POS = c(10537, 12180, 12372))
#' df_list <- prepare_list(ph_list, variants, phenotypes = "phenotype")
#' res <- build_model(df_list, outcomes = "phenotype", model = "logistic")
#' head(res)
#' @seealso forest_plot
#' @references
#' **stats**: R Core Team. (2025). *R: A Language and Environment for Statistical Computing.*
#'  R Foundation for Statistical Computing, Vienna, Austria.
#'  <https://www.R-project.org/>
#'
#' **dplyr**: Wickham, H., François, R., Henry, L., & Müller, K. (2023).
#'  *dplyr: A Grammar of Data Manipulation.* R package version 1.x.
#'  <https://CRAN.R-project.org/package=dplyr>
#' @export
build_model <- function(df_list,
                        outcomes,
                        dosage = "Dosage",
                        covars = character(),
                        model  = c("logistic","linear")) {
  model <- match.arg(model)
  if (is.data.frame(df_list)) {
    tables <- list(table1 = df_list)
  } else if (is.list(df_list) && length(df_list) > 0L) {
    ok <- vapply(df_list, is.data.frame, logical(1))
    if (!all(ok)) stop("All elements of `df_list` must be data.frames (e.g., from prepare_list()).")
    tables <- df_list
    if (is.null(names(tables))) names(tables) <- paste0("table", seq_along(tables))
  } else stop("`df_list` must be a data.frame or a non-empty list of data.frames.")
  table_ids <- names(tables)
  if (length(outcomes) == 1L) {
    outcome <- outcomes[1L]
  } else {
    u_out <- unique(outcomes)
    if (length(u_out) != 1L) stop("Joint modelling requires a single outcome name used in all tables.")
    outcome <- u_out[1L]
  }
  need_cols <- c(outcome, dosage, covars)
  miss_any  <- lapply(tables, function(dat) setdiff(need_cols, names(dat)))
  if (any(vapply(miss_any, length, integer(1L)) > 0L)) {
    res_err <- Map(
      f = function(miss, id) {
        data.frame(
          table_id = id, outcome = outcome, model = model, n = NA_integer_,
          term = dosage, beta = NA_real_, se = NA_real_, p = NA_real_,
          ci_lo = NA_real_, ci_hi = NA_real_,
          OR = NA_real_, OR_lo = NA_real_, OR_hi = NA_real_,
          error = paste0("Missing columns: ", paste(miss, collapse = ", ")),
          stringsAsFactors = FALSE)
      },
      miss_any, table_ids)
    out <- do.call(rbind, res_err)
    rownames(out) <- NULL
    return(out)
  }
  n_rows <- vapply(tables, nrow, integer(1L))
  if (length(unique(n_rows)) != 1L)
    stop("All tables in `df_list` must have the same number of rows for joint modelling.")
  base_dat <- tables[[1L]]
  model_dat <- base_dat[, c(outcome, covars), drop = FALSE]
  geno_mat <- lapply(tables, function(dat) dat[[dosage]])
  geno_df  <- as.data.frame(do.call(cbind, geno_mat), stringsAsFactors = FALSE)
  colnames(geno_df) <- table_ids
  model_dat <- cbind(model_dat, geno_df)
  rhs_terms <- c(covars, table_ids)
  rhs <- paste(rhs_terms, collapse = " + ")
  fm  <- stats::as.formula(paste0(outcome, " ~ ", rhs))
  use_cols <- c(outcome, covars, table_ids)
  mf <- stats::model.frame(
    stats::as.formula(paste0("~", paste(use_cols, collapse = "+"))),
    data = model_dat, na.action = stats::na.omit)
  n_obs <- nrow(mf)
  if (n_obs < max(3L, length(rhs_terms) + 2L)) {
    res_err <- lapply(table_ids, function(id) {
      data.frame(
        table_id = id, outcome = outcome, model = model, n = n_obs,
        term = dosage, beta = NA_real_, se = NA_real_, p = NA_real_,
        ci_lo = NA_real_, ci_hi = NA_real_,
        OR = NA_real_, OR_lo = NA_real_, OR_hi = NA_real_,
        error = "Not enough observations after NA omission.",
        stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, res_err)
    rownames(out) <- NULL
    return(out)
  }
  if (model == "linear") {
    fit <- try(stats::lm(fm, data = model_dat), silent = TRUE)
  } else {
    X <- stats::model.matrix(fm, data = model_dat)
    qrX <- qr(X)
    if (qrX$rank < ncol(X)) {
      keep_idx <- qrX$pivot[seq_len(qrX$rank)]
      keep_names <- colnames(X)[keep_idx]
      keep_rhs <- setdiff(keep_names, "(Intercept)")
      if (length(keep_rhs) == 0L) {
        res_err <- lapply(table_ids, function(id) {
          data.frame(
            table_id = id, outcome = outcome, model = model, n = n_obs,
            term = dosage, beta = NA_real_, se = NA_real_, p = NA_real_,
            ci_lo = NA_real_, ci_hi = NA_real_,
            OR = NA_real_, OR_lo = NA_real_, OR_hi = NA_real_,
            error = "All variant terms are collinear after removing NAs.",
            stringsAsFactors = FALSE)
        })
        out <- do.call(rbind, res_err)
        rownames(out) <- NULL
        return(out)
      }
      fm <- stats::as.formula(paste(outcome, "~", paste(keep_rhs, collapse = " + ")))
    }
    fit <- try(logistf::logistf(fm, data = model_dat), silent = TRUE)
  }
  if (inherits(fit, "try-error")) {
    msg <- as.character(fit)
    res_err <- lapply(table_ids, function(id) {
      data.frame(
        table_id = id, outcome = outcome, model = model, n = n_obs,
        term = dosage, beta = NA_real_, se = NA_real_, p = NA_real_,
        ci_lo = NA_real_, ci_hi = NA_real_,
        OR = NA_real_, OR_lo = NA_real_, OR_hi = NA_real_,
        error = paste0("Firth logistic failed: ", msg),
        stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, res_err)
    rownames(out) <- NULL
    return(out)
  }
  if (model == "linear") {
    co <- summary(fit)$coef
    res_list <- lapply(table_ids, function(id) {
      if (!(id %in% rownames(co))) {
        return(data.frame(
          table_id = id, outcome = outcome, model = model, n = n_obs,
          term = dosage, beta = NA_real_, se = NA_real_, p = NA_real_,
          ci_lo = NA_real_, ci_hi = NA_real_,
          OR = NA_real_, OR_lo = NA_real_, OR_hi = NA_real_,
          error = "Dosage term for this table_id not found in joint model.",
          stringsAsFactors = FALSE))
      }
      est <- co[id, "Estimate"]
      se  <- co[id, "Std. Error"]
      p   <- co[id, "Pr(>|t|)"]
      ci  <- try(stats::confint(fit, parm = id), silent = TRUE)
      if (inherits(ci, "try-error")) ci <- c(NA_real_, NA_real_)
      data.frame(
        table_id = id, outcome = outcome, model = model, n = n_obs,
        term = dosage, beta = unname(est), se = unname(se), p = unname(p),
        ci_lo = unname(ci[1]), ci_hi = unname(ci[2]),
        OR = NA_real_, OR_lo = NA_real_, OR_hi = NA_real_,
        error = NA_character_, stringsAsFactors = FALSE)
    })
  } else {
    coefs <- coef(fit)
    vc    <- try(stats::vcov(fit), silent = TRUE)
    res_list <- lapply(table_ids, function(id) {
      if (!(id %in% names(coefs))) {
        return(data.frame(
          table_id = id, outcome = outcome, model = model, n = n_obs,
          term = dosage, beta = NA_real_, se = NA_real_, p = NA_real_,
          ci_lo = NA_real_, ci_hi = NA_real_,
          OR = NA_real_, OR_lo = NA_real_, OR_hi = NA_real_,
          error = "Term dropped or not estimable (collinearity / sparse data).",
          stringsAsFactors = FALSE))
      }
      est <- unname(coefs[id])
      if (inherits(vc, "try-error")) {
        se <- NA_real_
      } else {
        se <- sqrt(vc[id, id])
      }
      if (!is.finite(se) || se <= 0) {
        return(data.frame(
          table_id = id, outcome = outcome, model = model, n = n_obs,
          term = dosage, beta = NA_real_, se = NA_real_, p = NA_real_,
          ci_lo = NA_real_, ci_hi = NA_real_,
          OR = NA_real_, OR_lo = NA_real_, OR_hi = NA_real_,
          error = "Unstable standard error in Firth model.",
          stringsAsFactors = FALSE))
      }
      zval <- est / se
      pval <- 2 * stats::pnorm(-abs(zval))
      zcrit <- stats::qnorm(0.975)
      lo_logit <- est - zcrit * se
      hi_logit <- est + zcrit * se
      OR    <- exp(est)
      OR_lo <- exp(lo_logit)
      OR_hi <- exp(hi_logit)
      data.frame(
        table_id = id, outcome = outcome, model = model, n = n_obs,
        term = dosage, beta = est, se = se, p = pval,
        ci_lo = NA_real_, ci_hi = NA_real_,
        OR = OR, OR_lo = OR_lo, OR_hi = OR_hi,
        error = NA_character_, stringsAsFactors = FALSE)
    })
  }
  res <- do.call(rbind, res_list)
  rownames(res) <- NULL
  res
}
