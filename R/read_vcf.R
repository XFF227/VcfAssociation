#' Read a VCF and return only GWAS-essential tables
#'
#' @description
#' Minimal VCF reader for GWAS pipelines.
#' Produces:
#'   (1) variant table with CHROM, POS, REF, ALT (no QUAL/FILTER/INFO),
#'   (2) long-format genotype table with per-sample ALT(1) dosage.
#'
#' @param vcf_file Path to .vcf or .vcf.gz
#' @param samples Optional character vector of sample IDs to keep (subset)
#' @param genome  Optional genome string for VariantAnnotation::readVcf (e.g., "hg38")
#'
#' @return list(variants, genotypes)
#'   - variants: data.frame(CHROM, POS, REF, ALT)
#'   - genotypes: data.frame(CHROM, POS, REF, ALT, sample, Dosage)
#' @export
read_vcf <- function(vcf_file, samples = NULL, genome = NULL) {
  
  if (!file.exists(vcf_file)) stop("VCF file not found: ", vcf_file)
  
  # Read VCF with minimal overhead
  vcf <- if (is.null(genome)) {
    VariantAnnotation::readVcf(vcf_file)
  } else {
    VariantAnnotation::readVcf(vcf_file, genome = genome)
  }
  
  # Subset samples if requested
  all_samps <- rownames(SummarizedExperiment::colData(vcf))
  if (!is.null(samples)) {
    miss <- setdiff(samples, all_samps)
    if (length(miss)) stop("Samples not in VCF: ", paste(miss, collapse = ", "))
    vcf <- vcf[, samples, drop = FALSE]
  }
  
  # ---------- Variant table (GWAS essentials only) ----------
  rr <- SummarizedExperiment::rowRanges(vcf)
  
  # First ALT only (character), safe from CharacterList
  alt_list <- VariantAnnotation::alt(vcf)
  ALT1 <- vapply(
    alt_list,
    function(x) if (length(x) == 0) NA_character_ else as.character(x[[1]]),
    character(1)
  )
  
  variants <- data.frame(
    CHROM = as.character(GenomicRanges::seqnames(rr)),
    POS   = as.integer(GenomicRanges::start(rr)),
    REF   = as.character(VariantAnnotation::ref(vcf)),
    ALT   = ALT1,
    stringsAsFactors = FALSE
  )
  
  # ---------- Genotype dosage (GT -> 0/1/2 for ALT(1)) ----------
  GT <- VariantAnnotation::geno(vcf)$GT
  if (is.null(GT)) stop("FORMAT/GT not found in VCF.")
  
  dosage_from_gt <- function(gt_str) {
    if (is.na(gt_str) || gt_str == ".") return(NA_real_)
    gt_str <- sub("\\|", "/", gt_str)             
    a <- strsplit(gt_str, "/", fixed = TRUE)[[1]]
    if (length(a) != 2) return(NA_real_)
    suppressWarnings(nums <- as.integer(a))
    if (any(is.na(nums))) return(NA_real_)
    sum(nums == 1)
  }

  dos_mat <- apply(GT, 2, function(col_gt) vapply(col_gt, dosage_from_gt, numeric(1)))
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
  
  return(list(variants = variants, genotypes = g))
}
