#' prepare_list
#' @title Prepare per-variant data sets for joint modeling
#'
#' @description
#' This function splits the merged phenotype–genotype data from read_phenotypes()
#' into smaller data frames, each representing one variant. Each sub-data frame
#' includes genotype dosage, phenotypes, covariates, and sample IDs.
#' The result is a list used by build_model() for gene- or region-level tests.
#'
#' @param ph Output list from read_phenotypes().
#' @param variants Optional data frame specifying CHROM and POS of variants to include.
#' @param phenotypes Character vector of phenotype columns to include.
#' @param covars Optional covariate columns.
#' @param dosage_col Name of dosage column.
#' @param chrom_col Name of chromosome column.
#' @param pos_col Name of position column.
#' @param id_col Name of sample ID column.
#'
#' @return A named list of data frames, one per variant.
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
#' variants <- data.frame(CHROM = "chr12", POS = c(10537, 12180))
#' df_list <- prepare_list(ph_list, variants, phenotypes = "phenotype")
#' names(df_list)
#' @seealso build_model
#' @references
#' **dplyr**: Wickham, H., François, R., Henry, L., & Müller, K. (2023).
#'  *dplyr: A Grammar of Data Manipulation.* R package version 1.x.
#'  <https://CRAN.R-project.org/package=dplyr>
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

