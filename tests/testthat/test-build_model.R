# tests/testthat/test-build_model.R
context("build_model")

test_that("build_model fits models and returns a summary data frame", {
  vcf_path <- system.file("extdata", "toy.vcf", package = "VcfAssociation", mustWork = TRUE)
  # Prepare a small df_list with two variants for testing
  geno_data <- read_vcf(vcf_path)$genotypes
  pheno_df <- generate_phenotype(vcf_path, chrom = "chr12", pos = 11161, model = "carrier")
  tmpfile <- tempfile(fileext = ".csv")
  write.csv(pheno_df, tmpfile, row.names = FALSE)
  ph <- read_phenotypes(tmpfile, id_col = "sample", genotypes = geno_data)
  # Select two variant positions for analysis
  variants_df <- unique(ph$merged[, c("CHROM", "POS")])[1:10, ]
  df_list_small <- prepare_list(ph, variants = variants_df, phenotypes = "phenotype")
  
  # Fit logistic models for the selected variants
  res_df <- build_model(df_list_small, outcomes = "phenotype", covars = character(), model = "logistic")
  
  # The result should be a data frame with one row per variant × outcome
  expect_true(is.data.frame(res_df))
  expect_equal(nrow(res_df), length(df_list_small) * 1)  # two variants * one outcome
  # Expected columns in summary (table_id identifies variant, outcome name, model type, effect estimates)
  expect_true(all(c("table_id", "outcome", "model", "beta", "se", "p") %in% names(res_df)))
  expect_setequal(res_df$table_id, names(df_list_small))
  expect_equal(unique(res_df$outcome), "phenotype")
  expect_equal(unique(res_df$model), "logistic")
  
  # p-values should be between 0 and 1
  expect_true(all(res_df$p >= 0 & res_df$p <= 1, na.rm = TRUE))
  # Check that each result row corresponds to the correct variant and outcome
  for (i in seq_len(nrow(res_df))) {
    variant_name <- res_df$table_id[i]
    # The variant name should be present in the df_list and outcome should match
    expect_true(variant_name %in% names(df_list_small))
    expect_equal(res_df$outcome[i], "phenotype")
  }
})
