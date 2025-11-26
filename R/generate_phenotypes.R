#' generate_phenotype
#' @title Generate simple phenotype data from a single variant
#'
#' @description
#' This function is the second step in the workflow. It selects one variant
#' (identified by chromosome and position) from the VCF data and converts its
#' ALT-allele dosage into a phenotype. Depending on the model:
#' - "carrier" gives binary 0/1 for presence of any ALT allele,
#' - "additive" uses the dosage value directly (0, 1, 2).
#'
#' The result can be saved as a CSV file and used later by read_phenotypes().
#'
#' @param vcf_path Path to the input VCF file.
#' @param chrom Chromosome of the variant to extract.
#' @param pos Position (1-based) of the variant to extract.
#' @param model Either "carrier" or "additive".
#' @param out_csv Optional path to save the generated phenotype table as CSV.
#' @param sample_col Name of the sample ID column in the output.
#'
#' @return A data frame with sample IDs and phenotype values.
#'
#' @examples
#' vcf_path <- system.file("extdata", "toy.vcf",
#'                         package = "VcfAssociation", mustWork = TRUE)
#' ph <- generate_phenotype(vcf_path,
#'                          chrom = "chr12", pos = 11161,
#'                          model = "carrier")
#' head(ph)
#'
#' tmp <- tempfile(fileext = ".csv")
#' ph_add <- generate_phenotype(vcf_path,
#'                              chrom = "chr12", pos = 11161,
#'                              model = "additive",
#'                              out_csv = tmp)
#' head(ph_add)
#' @seealso read_vcf, read_phenotypes
#' @references
#' **dplyr**: Wickham, H., François, R., Henry, L., & Müller, K. (2023).
#'  *dplyr: A Grammar of Data Manipulation.* R package version 1.x.
#'  <https://CRAN.R-project.org/package=dplyr>
#'
#' **stats**: R Core Team. (2025). *R: A Language and Environment for Statistical Computing.*
#'  R Foundation for Statistical Computing, Vienna, Austria.
#'  <https://www.R-project.org/>
#' @export
generate_phenotype <- function(vcf_path, chrom, pos, model = c("carrier","additive"),
                               out_csv = NULL, sample_col = "sample") {
  model <- match.arg(model)
  dat <- read_vcf(vcf_path)
  g <- dat$genotypes
  sub <- g[g$CHROM == chrom & g$POS == pos, c("sample","Dosage")]
  if (!nrow(sub)) stop("Variant not found: ", chrom, ":", pos)
  
  phenotype <- if (model == "carrier") as.integer(sub$Dosage > 0L) else as.integer(sub$Dosage)
  res <- stats::setNames(data.frame(sub$sample, phenotype, check.names = FALSE),
                         c(sample_col, "phenotype"))
  if (!is.null(out_csv)) readr::write_csv(res, out_csv)
  res
}
