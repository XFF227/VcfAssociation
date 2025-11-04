# tests/testthat/test-generate_phenotype.R
context("generate_phenotype")

test_that("generate_phenotype extracts carrier/additive phenotypes correctly", {
  vcf_path <- system.file("extdata", "toy.vcf", package = "VcfAssociation", mustWork = TRUE)
  
  # Generate phenotype for a known variant (chr12:11161) in both carrier and additive modes
  ph_carrier <- generate_phenotype(vcf_path, chrom = "chr12", pos = 11161, model = "carrier")
  ph_additive <- generate_phenotype(vcf_path, chrom = "chr12", pos = 11161, model = "additive")
  
  # The result should be a data frame with sample identifier and phenotype value
  expect_true(is.data.frame(ph_carrier))
  expect_true(is.data.frame(ph_additive))
  expect_true(all(c("sample", "phenotype") %in% names(ph_carrier)))
  expect_true(all(c("sample", "phenotype") %in% names(ph_additive)))
  
  # Phenotype values in carrier model should be 0/1 (integer), and in additive model 0/1/2
  expect_true(all(ph_carrier$phenotype %in% c(0, 1, NA)))
  expect_true(all(ph_additive$phenotype %in% c(0, 1, 2, NA)))
  
  # Ensure that there is at least one carrier (phenotype == 1) and at least one non-carrier (0)
  expect_true(1 %in% ph_carrier$phenotype, info = "Should identify at least one carrier")
  expect_true(0 %in% ph_carrier$phenotype, info = "Should identify at least one non-carrier")
  
  # In additive model, check that the maximum phenotype is 2 (homozygous alt present)
  expect_equal(max(ph_additive$phenotype, na.rm = TRUE), 2)
  
  # Missing genotypes in VCF should result in NA phenotypes for those samples
  expect_true(any(is.na(ph_carrier$phenotype)) || any(is.na(ph_additive$phenotype)))
})
