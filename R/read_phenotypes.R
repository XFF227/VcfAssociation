#' read_phenotypes
#' @title Read phenotype data and merge with genotype data
#'
#' @description
#' This function reads a phenotype file (such as one produced by generate_phenotype())
#' and merges it with the genotype table returned by read_vcf(). The merge is
#' performed by sample ID. It produces a unified data frame containing
#' both phenotype and genotype information for each sample.
#'
#' It also reports "gaps": samples that appear only in one table, helping
#' identify ID mismatches early in the analysis.
#'
#' @param pheno_path Path to the phenotype file (CSV or TSV).
#' @param id_col Name of the sample ID column in the phenotype file.
#' @param genotypes Optional genotype table from read_vcf().
#' @param by Mapping of sample ID columns between phenotype and genotype tables.
#'
#' @importFrom readr read_csv read_tsv
#' @importFrom dplyr left_join
#' @importFrom tools file_ext
#' @return A list containing:
#' - merged: combined phenotype–genotype table.
#' - gaps: unmatched sample IDs.
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
#' head(ph_list$merged)
#' @seealso gwas_single, prepare_list
#' @references
#' **readr**: Wickham, H., Hester, J., & Bryan, J. (2023).
#'  *readr: Read Rectangular Text Data.* R package version 2.x.
#'  <https://CRAN.R-project.org/package=readr>
#'
#' **dplyr**: Wickham, H., François, R., Henry, L., & Müller, K. (2023).
#'  *dplyr: A Grammar of Data Manipulation.* R package version 1.x.
#'  <https://CRAN.R-project.org/package=dplyr>
#' @export
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
