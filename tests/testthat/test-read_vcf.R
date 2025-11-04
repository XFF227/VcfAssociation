# tests/testthat/test-read_vcf.R
test_that("read_vcf reads VCF file and returns expected structure", {
  # Use the toy VCF file included in inst/extdata
  vcf_path <- system.file("extdata", "toy.vcf", package = "VcfAssociation", mustWork = TRUE)
  result <- read_vcf(vcf_path)
  
  # The result should be a list with named elements 'variants' and 'genotypes'
  expect_type(result, "list")
  expect_named(result, c("variants", "genotypes"))
  
  # Both 'variants' and 'genotypes' should be data frames
  expect_true(is.data.frame(result$variants))
  expect_true(is.data.frame(result$genotypes))
  
  # The variants data frame should contain the basic variant columns
  expect_true(all(c("CHROM", "POS", "REF", "ALT") %in% names(result$variants)))
  
  # The genotypes data frame should contain Sample and Dosage columns (long format)
  expect_true(all(c("sample", "Dosage") %in% names(result$genotypes)))
  
  # The toy VCF should include variant chr12:11161 – verify it appears in the variants table
  pos11161 <- subset(result$variants, CHROM == "chr12" & POS == 11161)
  expect_equal(nrow(pos11161), 1)  # exactly one variant at chr12:11161
  
  # Also check that genotypes for chr12:11161 are present for the expected number of samples
  geno_11161 <- subset(result$genotypes, CHROM == "chr12" & POS == 11161)
  unique_samples <- length(unique(result$genotypes$sample))
  expect_equal(nrow(geno_11161), unique_samples)  # one row per sample (including NA genotypes)
  
  # Dosage values should be 0, 1, 2 or NA_real_
  expect_true(all(geno_11161$Dosage %in% c(0, 1, 2, NA)))
})
