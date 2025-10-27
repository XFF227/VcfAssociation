#' annotate_variant
#' @title Annotate variants by gene regions (offline table or TxDb)
#' @description Use a bed-like table (CHROM, START, END, GENE) or TxDb exons to map positions to genes.
#' @param variants Data frame from \code{read_vcf()}$variants.
#' @param annot_tbl Optional data.frame with CHROM, START, END, GENE.
#' @param txdb Optional TxDb for exon overlaps.
#' @param orgdb Optional OrgDb to map Entrez IDs to SYMBOL.
#' @return Data frame with gene annotation columns.
#' @export
#' @examples
#' \dontrun{
#' anno <- annotate_variant(v$variants, annot_tbl = my_gene_bed)
#' }
annotate_variant <- function(variants, annot_tbl = NULL, txdb = NULL, orgdb = NULL) {
  stopifnot(is.data.frame(variants), all(c("CHROM","POS") %in% names(variants)))
  
  out <- variants
  
  if (is.data.frame(annot_tbl)) {
    need <- c("CHROM","START","END","GENE")
    if (!all(need %in% names(annot_tbl))) stop("annot_tbl must have: CHROM, START, END, GENE.")
    annot_tbl$CHROM <- as.character(annot_tbl$CHROM)
    
    q <- tibble::tibble(CHROM = out$CHROM, POS = out$POS, idx = seq_len(nrow(out)))
    hits <- dplyr::inner_join(q, annot_tbl, by = "CHROM")
    hits <- dplyr::filter(hits, POS >= START & POS <= END)
    
    out$GENE <- NA_character_
    if (nrow(hits)) {
      ag <- dplyr::summarise(dplyr::group_by(hits, idx), GENE = paste(unique(GENE), collapse = ";"))
      out$GENE[ag$idx] <- ag$GENE
    }
  } else if (!is.null(txdb)) {
    gr <- GenomicRanges::GRanges(seqnames = out$CHROM, ranges = IRanges::IRanges(out$POS, out$POS))
    ex <- GenomicFeatures::exonsBy(txdb, by = "gene")
    ov <- GenomicRanges::findOverlaps(gr, ex, ignore.strand = TRUE)
    
    out$GENE <- NA_character_
    if (length(ov)) {
      idx <- S4Vectors::queryHits(ov)
      gid <- names(ex)[S4Vectors::subjectHits(ov)]
      coll <- tapply(gid, idx, function(x) paste(unique(x), collapse = ";"))
      out$GENE[as.integer(names(coll))] <- unname(coll)
    }
    if (!is.null(orgdb)) {
      map <- function(g) {
        if (is.na(g)) return(NA_character_)
        ids <- strsplit(g, ";")[[1]]
        sy <- AnnotationDbi::mapIds(orgdb, keys = ids, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
        paste(na.omit(unname(sy)), collapse = ";")
      }
      out$SYMBOL <- vapply(out$GENE, map, FUN.VALUE = character(1))
    }
  } else {
    out$GENE <- NA_character_
  }
  
  out
}
