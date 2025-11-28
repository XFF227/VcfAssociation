# ============================================================
# Helper functions for build_model
# ============================================================

#' @title Create a standardized result row
#' @description Internal utility to create a single standardized data.frame row
#' used throughout `build_model()` for consistent outputs and error handling.
#' @keywords internal
make_row <- function(id,
                     outcome,
                     model,
                     dosage,
                     n,
                     beta   = NA_real_,
                     se     = NA_real_,
                     p      = NA_real_,
                     ci_lo  = NA_real_,
                     ci_hi  = NA_real_,
                     OR     = NA_real_,
                     OR_lo  = NA_real_,
                     OR_hi  = NA_real_,
                     error  = NA_character_) {
  data.frame(
    table_id = id,
    outcome  = outcome,
    model    = model,
    n        = n,
    term     = dosage,
    beta     = beta,
    se       = se,
    p        = p,
    ci_lo    = ci_lo,
    ci_hi    = ci_hi,
    OR       = OR,
    OR_lo    = OR_lo,
    OR_hi    = OR_hi,
    error    = error,
    stringsAsFactors = FALSE
  )
}

#' @title Bind list of data frames without rownames
#' @description Internal utility to rbind a list of data.frames and
#' remove row names for consistent formatting.
#' @keywords internal
bind_rows_list <- function(lst) {
  out <- do.call(rbind, lst)
  rownames(out) <- NULL
  out
}

#' @title Safe Firth logistic regression
#' @description Suppresses warnings and catches errors during
#' Firth logistic regression fitting using `logistf::logistf()`.
#' @param formula Model formula.
#' @param data Data frame.
#' @param maxit Maximum iterations.
#' @keywords internal
firth_fit_safe <- function(formula, data, maxit) {
  suppressWarnings(
    try(
      logistf::logistf(
        formula,
        data      = data,
        plcontrol = logistf::logistpl.control(maxit = maxit)
      ),
      silent = TRUE
    )
  )
}


# ============================================================
# Main function build_model
# ============================================================

#' build_model
#' @title Fit gene- or region-level association models
#'
#' @description
#' `build_model()` fits gene- or region-level association models on the
#' list of data.frames from `prepare_list()`. It merges phenotype, dosage,
#' and optional covariates into a modeling dataset and automatically chooses
#' between linear (`lm`) and Firth-penalized logistic (`logistf`) regression.
#'
#' For linear models, the function reports per-variant coefficients, standard
#' errors, p-values, and confidence intervals. For logistic models, it first
#' screens variants in univariate Firth regressions to remove non-convergent
#' or unstable ones, then fits a joint model on the remaining variants and
#' iteratively drops those causing instability.
#'
#' Results include beta, SE, p, CI, and odds ratios where applicable. Variants
#' excluded due to missing columns, insufficient sample size, or model failure
#' return `NA` estimates with a brief reason in the `error` column.
#'
#' @param df_list List of data frames prepared by prepare_list().
#' @param outcomes Phenotype column names.
#' @param dosage Name of dosage column.
#' @param covars Optional covariate names.
#' @param model Either "logistic" or "linear".
#' @param maxit Maximum iterations for Firth logistic (default 50).
#'
#' @return A data frame summarizing coefficients, standard errors, and p-values.
#' @export
build_model <- function(df_list,
                        outcomes,
                        dosage = "Dosage",
                        covars = character(),
                        model  = c("logistic", "linear"),
                        maxit  = 50) {
  
  model <- match.arg(model)
  
  ## -------- input checks --------
  if (is.data.frame(df_list)) {
    tables <- list(table1 = df_list)
  } else if (is.list(df_list) && length(df_list) > 0L) {
    ok <- vapply(df_list, is.data.frame, logical(1))
    if (!all(ok))
      stop("All elements of `df_list` must be data.frames (e.g., from prepare_list()).")
    tables <- df_list
    if (is.null(names(tables)))
      names(tables) <- paste0("table", seq_along(tables))
  } else {
    stop("`df_list` must be a data.frame or a non-empty list of data.frames.")
  }
  
  table_ids <- names(tables)
  
  if (length(outcomes) == 1L) {
    outcome <- outcomes[1L]
  } else {
    u_out <- unique(outcomes)
    if (length(u_out) != 1L)
      stop("Joint modelling requires a single outcome name.")
    outcome <- u_out[1L]
  }
  
  need_cols <- c(outcome, dosage, covars)
  miss_any  <- lapply(tables, function(dat) setdiff(need_cols, names(dat)))
  if (any(vapply(miss_any, length, integer(1L)) > 0L)) {
    res_err <- Map(
      f = function(miss, id) {
        if (length(miss) > 0L)
          make_row(id, outcome, model, dosage, n = NA_integer_, error = "Missing columns; removed.")
      },
      miss_any, table_ids
    )
    return(bind_rows_list(res_err))
  }
  
  n_rows <- vapply(tables, nrow, integer(1L))
  if (length(unique(n_rows)) != 1L)
    stop("All tables in `df_list` must have the same number of rows.")
  
  base_dat  <- tables[[1L]]
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
    data = model_dat,
    na.action = stats::na.omit
  )
  n_obs <- nrow(mf)
  
  if (n_obs < max(3L, length(rhs_terms) + 2L)) {
    res_err <- lapply(table_ids, function(id)
      make_row(id, outcome, model, dosage, n = n_obs, error = "Too few samples; removed."))
    return(bind_rows_list(res_err))
  }
  
  ## -------- linear branch --------
  if (model == "linear") {
    
    fit <- try(stats::lm(fm, data = model_dat), silent = TRUE)
    
    if (inherits(fit, "try-error")) {
      res_err <- lapply(table_ids, function(id)
        make_row(id, outcome, model, dosage, n = n_obs, error = "Linear model failed."))
      return(bind_rows_list(res_err))
    }
    
    co <- summary(fit)$coef
    
    res_list <- lapply(table_ids, function(id) {
      if (!(id %in% rownames(co))) {
        return(make_row(id, outcome, model, dosage, n = n_obs, error = "Linear term missing."))
      }
      est <- co[id, "Estimate"]
      se  <- co[id, "Std. Error"]
      p   <- co[id, "Pr(>|t|)"]
      ci  <- try(stats::confint(fit, parm = id), silent = TRUE)
      if (inherits(ci, "try-error")) ci <- c(NA_real_, NA_real_)
      make_row(id, outcome, model, dosage, n_obs, est, se, p, ci[1], ci[2])
    })
    
    return(bind_rows_list(res_list))
  }
  
  ## -------- logistic (Firth) branch --------
  bad_ids    <- character(0)
  bad_reason <- character(0)
  good_ids   <- character(0)
  
  for (id in table_ids) {
    terms1 <- c(covars, id)
    fm1    <- stats::as.formula(paste(outcome, "~", paste(terms1, collapse = " + ")))
    mf1    <- stats::model.frame(fm1, data = model_dat, na.action = stats::na.omit)
    n1     <- nrow(mf1)
    
    if (n1 < max(3L, length(terms1) + 2L)) {
      bad_ids    <- c(bad_ids, id)
      bad_reason <- c(bad_reason, "Too few samples.")
      next
    }
    
    fit1 <- firth_fit_safe(fm1, mf1, maxit)
    
    if (inherits(fit1, "try-error")) {
      bad_ids    <- c(bad_ids, id)
      bad_reason <- c(bad_reason, "Firth failed.")
      next
    }
    
    vc1 <- try(stats::vcov(fit1), silent = TRUE)
    if (inherits(vc1, "try-error") ||
        !is.finite(vc1[id, id]) || vc1[id, id] <= 0) {
      bad_ids    <- c(bad_ids, id)
      bad_reason <- c(bad_reason, "Unstable SE.")
      next
    }
    
    good_ids <- c(good_ids, id)
  }
  names(bad_reason) <- bad_ids
  
  if (length(good_ids) == 0L) {
    res_err <- lapply(table_ids, function(id) {
      msg <- if (id %in% bad_ids) bad_reason[[id]] else "Unstable Firth."
      make_row(id, outcome, model, dosage, n_obs, error = msg)
    })
    return(bind_rows_list(res_err))
  }
  
  joint_drop <- character(0)
  keep_ids   <- good_ids
  fit_joint  <- NULL
  
  repeat {
    rhs_joint <- paste(c(covars, keep_ids), collapse = " + ")
    fm_joint  <- stats::as.formula(paste(outcome, "~", rhs_joint))
    fit_joint <- firth_fit_safe(fm_joint, model_dat, maxit)
    
    if (inherits(fit_joint, "try-error")) {
      if (length(keep_ids) <= 1L) {
        joint_drop <- c(joint_drop, keep_ids)
        keep_ids   <- character(0)
        break
      }
      drop_var   <- tail(keep_ids, 1)
      joint_drop <- c(joint_drop, drop_var)
      keep_ids   <- setdiff(keep_ids, drop_var)
      next
    }
    
    vcj <- try(stats::vcov(fit_joint), silent = TRUE)
    if (inherits(vcj, "try-error")) {
      if (length(keep_ids) <= 1L) {
        joint_drop <- c(joint_drop, keep_ids)
        keep_ids   <- character(0)
        fit_joint  <- NULL
        break
      }
      drop_var   <- tail(keep_ids, 1)
      joint_drop <- c(joint_drop, drop_var)
      keep_ids   <- setdiff(keep_ids, drop_var)
      next
    }
    
    bad_in_joint <- keep_ids[!is.finite(diag(vcj)[keep_ids]) |
                               diag(vcj)[keep_ids] <= 0]
    if (length(bad_in_joint) == 0L) break
    if (length(keep_ids) <= length(bad_in_joint)) {
      joint_drop <- c(joint_drop, keep_ids)
      keep_ids   <- character(0)
      fit_joint  <- NULL
      break
    }
    joint_drop <- c(joint_drop, bad_in_joint)
    keep_ids   <- setdiff(keep_ids, bad_in_joint)
  }
  
  joint_drop <- unique(joint_drop)
  
  if (length(keep_ids) == 0L || is.null(fit_joint)) {
    res_err <- lapply(table_ids, function(id) {
      msg <- if (id %in% bad_ids) {
        bad_reason[[id]]
      } else if (id %in% joint_drop) {
        "Unstable Joined, removed"
      } else {
        "Unstable Firth."
      }
      make_row(id, outcome, model, dosage, n_obs, error = msg)
    })
    return(bind_rows_list(res_err))
  }
  
  coefs <- coef(fit_joint)
  vc    <- try(stats::vcov(fit_joint), silent = TRUE)
  
  res_list <- lapply(table_ids, function(id) {
    if (id %in% bad_ids)
      return(make_row(id, outcome, model, dosage, n_obs, error = bad_reason[[id]]))
    if (id %in% joint_drop)
      return(make_row(id, outcome, model, dosage, n_obs, error = "Joint Firth failed."))
    if (!(id %in% names(coefs)))
      return(make_row(id, outcome, model, dosage, n_obs, error = "Term dropped; removed."))
    if (inherits(vc, "try-error") ||
        !is.finite(vc[id, id]) || vc[id, id] <= 0)
      return(make_row(id, outcome, model, dosage, n_obs, error = "Unstable SE; removed."))
    
    est <- unname(coefs[id])
    se  <- sqrt(vc[id, id])
    z   <- est / se
    p   <- 2 * stats::pnorm(-abs(z))
    zc  <- stats::qnorm(0.975)
    lo_logit <- est - zc * se
    hi_logit <- est + zc * se
    make_row(id, outcome, model, dosage, n_obs, est, se, p,
             OR = exp(est), OR_lo = exp(lo_logit), OR_hi = exp(hi_logit))
  })
  
  bind_rows_list(res_list)
}