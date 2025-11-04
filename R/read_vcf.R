#' Read a VCF and return a genotype and variants tables
#'
#' @description
#' Minimal VCF reader for GWAS pipelines using the vcfR package (no Bioconductor).
#' Produces:
#'   (1) variant table with CHROM, POS, REF, ALT (no QUAL/FILTER/INFO),
#'   (2) long-format genotype table with per-sample ALT(1) dosage.
#'
#' @param vcf_file Path to .vcf or .vcf.gz
#' @param samples Optional character vector of sample IDs to keep (subset)
#' @param genome  Ignored (kept for API compatibility)
#'
#' @return list(variants, genotypes)
#'   - variants: data.frame(CHROM, POS, REF, ALT)
#'   - genotypes: data.frame(CHROM, POS, REF, ALT, sample, Dosage)
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
