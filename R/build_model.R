#' build_model
#' @title Fit models for one or multiple phenotypes and return a tidy summary table
#' @description
#' Batch wrapper of single-variant association. Accepts a list input (e.g., the
#' return of `prepare_list()`) and fits Linear or Logistic models of
#' `outcome ~ Dosage + covariates` for one or more phenotypes, returning a single
#' data.frame of summaries.
#'
#' @param df_list A named/unnamed list of data.frames (typically from prepare_list()).
#' @param outcomes Character vector of phenotype column name(s). Length can be 1
#'   (recycled) or match the number of tables in `df_list`.
#' @param dosage Dosage column name. Default "Dosage".
#' @param covars Optional character vector of covariate column names.
#' @param model One of "logistic","linear".
#' @return data.frame with one row per (table × outcome) containing:
#'   \code{table_id, outcome, model, n, term, beta, se, p, ci_lo, ci_hi, OR, OR_lo, OR_hi, error}.
#'   For linear models, OR columns are NA; for logistic, ci_lo/ci_hi are NA.
#' @export
build_model <- function(df_list,
                        outcomes,
                        dosage = "Dosage",
                        covars = character(),
                        model  = c("logistic","linear")) {
  model <- match.arg(model)
  
  # --- Normalize df_list: must be a list of data.frames -----------------------
  if (is.data.frame(df_list)) {
    tables <- list(table1 = df_list)
  } else if (is.list(df_list) && length(df_list) > 0) {
    ok <- vapply(df_list, is.data.frame, logical(1))
    if (!all(ok)) stop("All elements of `df_list` must be data.frames (e.g., from prepare_list()).")
    tables <- df_list
    if (is.null(names(tables))) names(tables) <- paste0("table", seq_along(tables))
  } else {
    stop("`df_list` must be a data.frame or a non-empty list of data.frames.")
  }
  
  # --- Recycle / align outcomes to tables -------------------------------------
  if (length(outcomes) == 1L) {
    outcomes <- rep(outcomes, length(tables))
  } else if (length(outcomes) != length(tables)) {
    stop("Length of `outcomes` must be 1 or equal to the number of tables in `df_list`.")
  }
  
  # --- Internal fitter for one (table, outcome) pair ---------------------------
  fit_one <- function(dat, outcome, dosage, covars, model, table_id) {
    req <- c(outcome, dosage, covars)
    miss <- setdiff(req, names(dat))
    if (length(miss)) {
      return(data.frame(
        table_id = table_id, outcome = outcome, model = model, n = NA_integer_,
        term = dosage, beta = NA_real_, se = NA_real_, p = NA_real_,
        ci_lo = NA_real_, ci_hi = NA_real_,
        OR = NA_real_, OR_lo = NA_real_, OR_hi = NA_real_,
        error = paste0("Missing columns: ", paste(miss, collapse = ", ")),
        stringsAsFactors = FALSE
      ))
    }
    
    # Build formula
    rhs_terms <- c(dosage, covars)
    rhs <- paste(rhs_terms, collapse = " + ")
    fm  <- stats::as.formula(paste0(outcome, " ~ ", rhs))
    
    # Keep complete cases on variables used
    use_cols <- c(outcome, dosage, covars)
    dat2 <- stats::model.frame(stats::as.formula(paste0("~", paste(use_cols, collapse="+"))),
                               data = dat, na.action = stats::na.omit)
    n_obs <- nrow(dat2)
    if (n_obs < max(3, length(rhs_terms) + 2L)) {
      return(data.frame(
        table_id = table_id, outcome = outcome, model = model, n = n_obs,
        term = dosage, beta = NA_real_, se = NA_real_, p = NA_real_,
        ci_lo = NA_real_, ci_hi = NA_real_,
        OR = NA_real_, OR_lo = NA_real_, OR_hi = NA_real_,
        error = "Not enough observations after NA omission.",
        stringsAsFactors = FALSE
      ))
    }
    
    if (model == "linear") {
      fit <- try(stats::lm(fm, data = dat2), silent = TRUE)
      if (inherits(fit, "try-error")) {
        return(data.frame(
          table_id = table_id, outcome = outcome, model = model, n = n_obs,
          term = dosage, beta = NA_real_, se = NA_real_, p = NA_real_,
          ci_lo = NA_real_, ci_hi = NA_real_,
          OR = NA_real_, OR_lo = NA_real_, OR_hi = NA_real_,
          error = as.character(fit),
          stringsAsFactors = FALSE
        ))
      }
      co <- summary(fit)$coef
      if (!(dosage %in% rownames(co))) {
        return(data.frame(
          table_id = table_id, outcome = outcome, model = model, n = n_obs,
          term = dosage, beta = NA_real_, se = NA_real_, p = NA_real_,
          ci_lo = NA_real_, ci_hi = NA_real_,
          OR = NA_real_, OR_lo = NA_real_, OR_hi = NA_real_,
          error = "Dosage term not found in model matrix.",
          stringsAsFactors = FALSE
        ))
      }
      est <- co[dosage, "Estimate"]
      se  <- co[dosage, "Std. Error"]
      p   <- co[dosage, "Pr(>|t|)"]
      ci  <- try(stats::confint(fit, parm = dosage), silent = TRUE)
      if (inherits(ci, "try-error")) ci <- c(NA_real_, NA_real_)
      return(data.frame(
        table_id = table_id, outcome = outcome, model = model, n = n_obs,
        term = dosage, beta = unname(est), se = unname(se), p = unname(p),
        ci_lo = unname(ci[1]), ci_hi = unname(ci[2]),
        OR = NA_real_, OR_lo = NA_real_, OR_hi = NA_real_,
        error = NA_character_, stringsAsFactors = FALSE
      ))
    }
    
    if (model == "logistic") {
      fit <- try(stats::glm(fm, data = dat2, family = stats::binomial()), silent = TRUE)
      if (inherits(fit, "try-error")) {
        return(data.frame(
          table_id = table_id, outcome = outcome, model = model, n = n_obs,
          term = dosage, beta = NA_real_, se = NA_real_, p = NA_real_,
          ci_lo = NA_real_, ci_hi = NA_real_,
          OR = NA_real_, OR_lo = NA_real_, OR_hi = NA_real_,
          error = as.character(fit), stringsAsFactors = FALSE
        ))
      }
      co <- summary(fit)$coef
      if (!(dosage %in% rownames(co))) {
        return(data.frame(
          table_id = table_id, outcome = outcome, model = model, n = n_obs,
          term = dosage, beta = NA_real_, se = NA_real_, p = NA_real_,
          ci_lo = NA_real_, ci_hi = NA_real_,
          OR = NA_real_, OR_lo = NA_real_, OR_hi = NA_real_,
          error = "Dosage term not found in model matrix.",
          stringsAsFactors = FALSE
        ))
      }
      est <- co[dosage, "Estimate"]
      se  <- co[dosage, "Std. Error"]
      p   <- co[dosage, "Pr(>|z|)"]
      z   <- stats::qnorm(0.975)
      lo_logit <- est - z * se
      hi_logit <- est + z * se
      return(data.frame(
        table_id = table_id, outcome = outcome, model = model, n = n_obs,
        term = dosage, beta = unname(est), se = unname(se), p = unname(p),
        ci_lo = NA_real_, ci_hi = NA_real_,
        OR    = exp(unname(est)),
        OR_lo = exp(unname(lo_logit)),
        OR_hi = exp(unname(hi_logit)),
        error = NA_character_, stringsAsFactors = FALSE
      ))
    }
  }
  
  res_list <- Map(
    f = function(tab, out, id) fit_one(tab, out, dosage, covars, model, table_id = id),
    tables, outcomes, names(tables)
  )
  
  res <- do.call(rbind, res_list)
  rownames(res) <- NULL
  res
}

