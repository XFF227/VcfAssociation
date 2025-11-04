# tests/testthat/test-gwas_single.R
test_that("gwas_single returns per-variant association rows and expected fields", {
  vcf_path <- system.file("extdata", "toy.vcf", package = "VcfAssociation", mustWork = TRUE)
  
  geno_data <- read_vcf(vcf_path)$genotypes
  pheno_df  <- generate_phenotype(vcf_path, chrom = "chr12", pos = 11161, model = "carrier")
  
  tmpfile <- tempfile(fileext = ".csv")
  write.csv(pheno_df, tmpfile, row.names = FALSE)
  
  ph <- read_phenotypes(tmpfile, id_col = "sample", genotypes = geno_data)
  
  result <- gwas_single(ph, pheno_col = "phenotype", covars = NULL, id_col = "sample")

  expect_true(is.data.frame(result))
  expect_true(nrow(result) >= 1)
  expect_true(all(c("CHROM", "POS", "beta", "se", "p") %in% names(result)))
  expect_true(!any(duplicated(result[c("CHROM", "POS")])))
  idx <- which(result$CHROM == "chr12" & result$POS == 11161)
  expect_true(length(idx) == 1)
  r <- result[idx, , drop = FALSE]
  expect_true(is.numeric(r$beta[1]) && !is.na(r$beta[1]) && is.finite(r$beta[1]))
  expect_true(is.numeric(r$se[1])   && !is.na(r$se[1])   && is.finite(r$se[1]))
  expect_true(is.numeric(r$p[1])    && !is.na(r$p[1])    && r$p[1] >= 0 && r$p[1] <= 1)
})
