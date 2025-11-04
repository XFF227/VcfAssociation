#' generate_phenotype
#' @title Derive phenotype from a specific variant
#' @description Extract a variant (chrom:pos) from VCF and output carrier/additive phenotype.
#' @param vcf_path Path to VCF.
#' @param chrom Chromosome.
#' @param pos 1-based position.
#' @param model "carrier" (0/1) or "additive" (0/1/2).
#' @param out_csv Optional path to write CSV.
#' @param sample_col Name of sample column in result.
#' @return data.frame with columns \code{sample_col} and \code{phenotype}.
#' @export
#' @examples
#' \dontrun{
#' ph2 <- generate_phenotype("demo.vcf.gz", "1", 123456, model="carrier", out_csv="pheno_from_snp.csv")
#' }
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
