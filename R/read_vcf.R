#' read_vcf
#' @title Read a VCF file into tidy data frames
#' @description Read a VCF (.vcf or .vcf.gz) and return `variants` (one row per site)
#'   and `genotypes` (long format: variant x sample) with computed dosage.
#' @param vcf_path Path to VCF file.
#' @param samples Optional character vector of sample IDs to subset.
#' @param info_fields Optional character vector of INFO tags to keep.
#' @return A list with \code{$variants} and \code{$genotypes} data.frames.
#' @export
#' @examples
#' \dontrun{
#' v <- read_vcf("demo.vcf.gz")
#' head(v$variants)
#' head(v$genotypes)
#' }
read_vcf <- function(vcf_path, samples = NULL, info_fields = NULL) {
  stopifnot(length(vcf_path) == 1L, is.character(vcf_path))
  if (!file.exists(vcf_path)) stop("VCF not found: ", vcf_path)
  
  vcf <- VariantAnnotation::readVcf(vcf_path)
  
  if (!is.null(samples)) {
    # Try to infer sample column name from colData if present
    cd <- SummarizedExperiment::colData(vcf)
    present <- colnames(vcf)
    keep <- intersect(samples, present)
    if (!length(keep)) stop("Requested samples not in VCF.")
    vcf <- vcf[, keep, drop = FALSE]
  }
  
  vr <- VariantAnnotation::rowRanges(vcf)
  info_df <- as.data.frame(VariantAnnotation::info(vcf))
  if (!is.null(info_fields)) {
    info_df <- info_df[, intersect(info_fields, colnames(info_df)), drop = FALSE]
  }
  
  variants <- data.frame(
    CHROM = as.character(GenomeInfoDb::seqnames(vr)),
    POS   = as.integer(GenomicRanges::start(vr)),
    ID    = if ("ID" %in% names(vr)) as.character(vr$ID) else NA_character_,
    REF   = as.character(Biostrings::DNAStringSet(GenomicRanges::mcols(vr)$REF)),
    ALT   = vapply(GenomicRanges::mcols(vr)$ALT, function(x) paste(as.character(x), collapse = ","), ""),
    stringsAsFactors = FALSE
  )
  if (ncol(info_df)) variants <- cbind(variants, info_df)
  
  gt <- VariantAnnotation::geno(vcf)$GT
  if (is.null(gt)) stop("GENO/GT missing in VCF.")
  
  gt_df <- tibble::as_tibble(gt, .name_repair = "minimal")
  gt_df$.__rowid__ <- seq_len(nrow(gt_df))
  
  long <- utils::stack(as.data.frame(gt_df[, setdiff(colnames(gt_df), ".__rowid__"), drop = FALSE]))
  names(long) <- c("GT", "Sample")
  long$Row <- rep(gt_df$.__rowid__, times = ncol(gt))
  
  variants$Row <- seq_len(nrow(variants))
  genotypes <- dplyr::left_join(long, variants, by = "Row")
  genotypes$Dosage <- ifelse(grepl("1[\\|/]1", genotypes$GT), 2L,
                             ifelse(grepl("0[\\|/]1|1[\\|/]0", genotypes$GT), 1L,
                                    ifelse(grepl("\\.", genotypes$GT), NA_integer_, 0L)))
  genotypes <- genotypes[, c("Sample","CHROM","POS","ID","REF","ALT","GT","Dosage")]
  
  list(variants = variants[, setdiff(names(variants), "Row"), drop = FALSE],
       genotypes = as.data.frame(genotypes))
}
