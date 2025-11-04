#' read_phenotypes
#' @title Read phenotype table and merge with genotypes
#' @param pheno_path CSV/TSV path.
#' @param id_col Sample ID column name in phenotype file.
#' @param genotypes Optional genotype long table from \code{read_vcf()}.
#' @param by Named character vector mapping phenotype id to genotype id (default sample->sample).
#' @return list(merged=data.frame, gaps=list(pheno_only, geno_only))
#' @export
#' @examples
#' \dontrun{
#' ph <- read_phenotypes("pheno.csv", id_col="sample_id", genotypes=v$genotypes,
#'                       by=c("sample_id"="sample"))
#' ph$gaps
#' head(ph$merged)
#' }
read_phenotypes <- function(pheno_path, id_col = "sample", genotypes = NULL, by = c("sample" = "sample")) {
  if (!file.exists(pheno_path)) stop("Phenotype not found: ", pheno_path)
  ext <- tolower(tools::file_ext(pheno_path))
  df <- if (ext %in% c("tsv","txt")) readr::read_tsv(pheno_path, show_col_types = FALSE) else readr::read_csv(pheno_path, show_col_types = FALSE)
  if (!id_col %in% names(df)) stop("id_col not found in phenotype file.")
  
  if (is.null(genotypes)) {
    return(list(merged = df, gaps = list(pheno_only = character(), geno_only = character())))
  }
  
  merged <- dplyr::left_join(df, genotypes, by = setNames(names(by), by))
  pheno_ids <- unique(df[[id_col]])
  geno_ids  <- unique(genotypes[[by[[1]]]])
  gaps <- list(pheno_only = setdiff(pheno_ids, geno_ids), geno_only = setdiff(geno_ids, pheno_ids))
  list(merged = merged, gaps = gaps)
}
