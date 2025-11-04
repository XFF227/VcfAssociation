# tests/testthat/test-read_phenotypes.R
test_that("read_phenotypes correctly merges phenotype and genotype data", {
  vcf_path <- system.file("extdata", "toy.vcf", package = "VcfAssociation", mustWork = TRUE)
  geno_data <- read_vcf(vcf_path)$genotypes  # get genotype data frame from toy VCF
  
  # Create a toy phenotype table from the carrier phenotype of chr12:11161
  pheno_df <- generate_phenotype(vcf_path, chrom = "chr12", pos = 11161, model = "carrier")
  
  # Write phenotype data to a temporary CSV (simulate an external phenotype file)
  tmpfile <- tempfile(fileext = ".csv")
  write.csv(pheno_df, tmpfile, row.names = FALSE)
  
  # Read and merge phenotypes with genotype data
  ph <- read_phenotypes(tmpfile, id_col = "sample", genotypes = geno_data)
  
  # The result should be a list with 'merged' data frame and 'gaps' list
  expect_type(ph, "list")
  expect_equal(names(ph), c("merged", "gaps"))
  expect_true(is.data.frame(ph$merged))
  expect_true(is.list(ph$gaps))
  
  # The merged data frame should include phenotype and Dosage columns (and others like CHROM, POS)
  expect_true("phenotype" %in% names(ph$merged))
  expect_true("Dosage" %in% names(ph$merged))
  
  # All rows in merged should have a phenotype (NA if phenotype missing for some genotypes)
  expect_equal(
    nrow(ph$merged),
    nrow(subset(geno_data, sample %in% ph$merged$sample))
  )
})
