#' Read a VCF (plain or bgz) into a tidy table
#' @param path Path to .vcf or .vcf.gz
#' @return data.frame with CHROM, POS, REF, ALT, SAMPLE, GT, DP, AD
#' @examples
#' v <- read_vcf(system.file("extdata","toy.vcf", package="VcfAssociation"))
#' head(v)
#' @export
read_vcf <- function(path){
  v <- vcfR::read.vcfR(path, verbose = FALSE)
  fix <- vcfR::getFIX(v)
  gt  <- vcfR::extract_GT(v)
  d <- data.frame(CHROM=fix[,"CHROM"], POS=as.integer(fix[,"POS"]),
                  REF=fix[,"REF"], ALT=fix[,"ALT"], check.names=FALSE)

  sm <- colnames(v@gt)[-1]
  out <- do.call(rbind, lapply(seq_along(sm), function(i){
    data.frame(d, SAMPLE=sm[i], GT=gt[,i], stringsAsFactors=FALSE)
  }))
  out
}