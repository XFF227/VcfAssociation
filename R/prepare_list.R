#' prepare_list
#' @title Prepare df_list input for build_model() from read_phenotypes() output
#' @description Extract one or multiple variant-level data.frames from ph$merged
#' @param ph The list returned by read_phenotypes()
#' @param variants data.frame with CHROM, POS; if NULL, returns one element with all rows
#' @param phenotypes character vector of phenotype column names
#' @param covars optional covariate names
#' @param dosage_col dosage column name (default "Dosage")
#' @param chrom_col chromosome column name (default "CHROM")
#' @param pos_col position column name (default "POS")
#' @param id_col sample id column name (default "sample")
#' @return named list of data.frames
#' @export
prepare_list <- function(ph,
                         variants = NULL,
                         phenotypes,
                         covars = character(),
                         dosage_col = "Dosage",
                         chrom_col = "CHROM",
                         pos_col = "POS",
                         id_col = "sample") {
  if (!is.list(ph) || is.null(ph$merged))
    stop("Input must be the list returned by read_phenotypes().")
  
  dat <- ph$merged
  
  req <- unique(c(id_col, chrom_col, pos_col, dosage_col, phenotypes, covars))
  miss <- setdiff(req, names(dat))
  if (length(miss)) stop("Missing required columns in ph$merged: ", paste(miss, collapse = ", "))
  if (is.null(variants)) {
    return(list(full = dat[, req, drop = FALSE]))
  }
  
  if (!all(c("CHROM","POS") %in% names(variants))) {
    stop("`variants` must contain columns CHROM and POS.")
  }
  
  df_list <- lapply(seq_len(nrow(variants)), function(i) {
    subset(dat,
           dat[[chrom_col]] == variants$CHROM[i] &
             dat[[pos_col]]   == variants$POS[i],
           select = intersect(req, names(dat)))
  })
  
  names(df_list) <- paste0("chr", variants$CHROM, "_pos", variants$POS)
  df_list
}

