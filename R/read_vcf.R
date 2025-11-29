#' read_vcf
#' @title Read and parse annotated VCF files for association analysis
#'
#' @description
#' This function is the first step in the pipeline. It reads an annotated VCF file
#' and extracts two key tables:
#' 1) a variant table with basic variant information (CHROM, POS, REF, ALT),
#' 2) a genotype table in long format with ALT-allele dosage for each sample.
#'
#' These tables form the foundation for all downstream steps, such as
#' generating phenotypes and merging genotype–phenotype data.
#'
#' @param vcf_file Path to the annotated VCF file (.vcf or .vcf.gz).
#' @param samples Optional character vector of sample IDs to retain.
#' @param genome Optional genome build name (not used internally, kept for compatibility).
#'
#' @importFrom vcfR read.vcfR extract.gt
#' @return A list containing:
#' - variants: data frame with CHROM, POS, REF, and ALT.
#' - genotypes: long-format data frame with per-sample Dosage values.
#'
#' @examples
#' vcf_path <- system.file("extdata", "toy.vcf",
#'                         package = "VcfAssociation", mustWork = TRUE)
#' vcf <- read_vcf(vcf_path)
#' head(vcf$variants)
#' head(vcf$genotypes)
#' @seealso generate_phenotype, read_phenotypes
#' @references
#' **vcfR**: Knaus, Brian J.; Grünwald, Niklaus J. (2017).
#'  vcfR: a package to manipulate and visualize variant call format data in R.
#'  *Molecular Ecology Resources* 17(1): 44–53.
#'  <http://dx.doi.org/10.1111/1755-0998.12549>
#' @export
read_vcf <- function(vcf_file, samples = NULL, genome = NULL) {
  if (!file.exists(vcf_file)) stop("VCF file not found: ", vcf_file)
  
  vcf <- vcfR::read.vcfR(vcf_file, verbose = FALSE)
  
  fix <- vcf@fix
  variants <- data.frame(
    CHROM = as.character(fix[, "CHROM"]),
    POS   = as.integer(fix[, "POS"]),
    REF   = as.character(fix[, "REF"]),
    ALT   = as.character(fix[, "ALT"]),
    stringsAsFactors = FALSE
  )
  
  GT <- vcfR::extract.gt(vcf, element = "GT", as.numeric = FALSE)
  if (is.null(GT)) stop("FORMAT/GT not found in VCF.")
  
  all_samps <- colnames(GT)
  if (!is.null(samples)) {
    miss <- setdiff(samples, all_samps)
    if (length(miss)) stop("Samples not in VCF: ", paste(miss, collapse = ", "))
    GT <- GT[, samples, drop = FALSE]
    all_samps <- samples
  }
  
  dosage_from_gt <- function(gt_str) {
    if (is.na(gt_str) || gt_str == "." || gt_str == "./.") return(NA_real_)
    gt_str <- gsub("\\|", "/", gt_str)
    parts <- strsplit(gt_str, "/", fixed = TRUE)[[1]]
    if (length(parts) != 2) return(NA_real_)
    suppressWarnings(nums <- as.integer(parts))
    if (any(is.na(nums))) return(NA_real_)
    sum(nums == 1L)
  }
  
  dos_mat <- apply(GT, 2, function(col_gt) vapply(col_gt, dosage_from_gt, numeric(1)))
  if (is.null(dim(dos_mat))) {
    # Single-sample edge case: make it a matrix with one column
    dos_mat <- matrix(dos_mat, ncol = 1, dimnames = list(NULL, all_samps))
  } else {
    colnames(dos_mat) <- all_samps
  }
  
  locus_key <- paste0(variants$CHROM, ":", variants$POS, "_", variants$REF, ">", variants$ALT)
  rownames(dos_mat) <- locus_key
  
  g <- as.data.frame(as.table(dos_mat), stringsAsFactors = FALSE)
  colnames(g) <- c("VariantKey", "sample", "Dosage")
  
  key_df <- data.frame(
    VariantKey = locus_key,
    CHROM = variants$CHROM,
    POS   = variants$POS,
    REF   = variants$REF,
    ALT   = variants$ALT,
    stringsAsFactors = FALSE
  )
  
  g <- merge(g, key_df, by = "VariantKey", all.x = TRUE, sort = FALSE)
  g <- g[, c("CHROM", "POS", "REF", "ALT", "sample", "Dosage")]
  
  list(variants = variants, genotypes = g)
}
