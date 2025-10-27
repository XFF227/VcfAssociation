#' Toy GWAS scan over positions (presence/absence)
#' @param var_df annotated variant table
#' @param pheno data.frame(SAMPLE, PHENO)
#' @return data.frame(CHR, POS, P, BETA)
#' @export
assoc_gwas <- function(var_df, pheno){
  pos_list <- unique(var_df[,c("CHROM","POS")])
  res <- lapply(seq_len(nrow(pos_list)), function(i){
    sl <- subset(var_df, CHROM==pos_list$CHROM[i] & POS==pos_list$POS[i])
    carr <- aggregate((GT!="0/0")~SAMPLE, sl, mean); colnames(carr) <- c("SAMPLE","CARR")
    d <- merge(pheno, carr, by="SAMPLE", all.x=TRUE); d$CARR[is.na(d$CARR)] <- 0
    fit <- stats::glm(PHENO ~ CARR, data=d, family=stats::binomial())
    c(CHR=pos_list$CHROM[i], POS=pos_list$POS[i],
      P=summary(fit)$coef["CARR","Pr(>|z|)"], BETA=coef(fit)["CARR"])
  })
  data.frame(do.call(rbind, res), row.names=NULL)
}
