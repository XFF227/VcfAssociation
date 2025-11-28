test_that("build_model (logistic) fits Firth models and returns a well-formed summary", {
  vcf_path <- system.file(
    "extdata", "toy.vcf",
    package = "VcfAssociation",
    mustWork = TRUE
  )
  
  geno_data <- read_vcf(vcf_path)$genotypes
  pheno_df  <- generate_phenotype(
    vcf_path,
    chrom = "chr12",
    pos   = 11161,
    model = "carrier"
  )
  
  tmpfile <- tempfile(fileext = ".csv")
  write.csv(pheno_df, tmpfile, row.names = FALSE)
  ph <- read_phenotypes(tmpfile, id_col = "sample", genotypes = geno_data)
  
  variants_df   <- unique(ph$merged[, c("CHROM", "POS")])[1:10, ]
  df_list_small <- prepare_list(
    ph,
    variants   = variants_df,
    phenotypes = "phenotype"
  )
  
  res_df <- build_model(
    df_list_small,
    outcomes = "phenotype",
    covars   = character(),
    model    = "logistic"
  )
  
  ## basic structure
  expect_true(is.data.frame(res_df))
  expect_equal(nrow(res_df), length(df_list_small) * 1L)
  expect_setequal(res_df$table_id, names(df_list_small))
  
  ## expected columns
  expect_true(all(c(
    "table_id", "outcome", "model", "n",
    "term", "beta", "se", "p",
    "ci_lo", "ci_hi",
    "OR", "OR_lo", "OR_hi",
    "error"
  ) %in% names(res_df)))
  
  ## labels
  expect_equal(unique(res_df$outcome), "phenotype")
  expect_equal(unique(res_df$model), "logistic")
  
  ## n is non-negative when not NA (can be 0 in 'too few' cases)
  expect_true(all(is.na(res_df$n) | res_df$n >= 0))
  
  ## p-values within [0, 1] when non-NA
  expect_true(all(res_df$p[!is.na(res_df$p)] >= 0))
  expect_true(all(res_df$p[!is.na(res_df$p)] <= 1))
  
  ## ORs positive when non-NA
  expect_true(all(res_df$OR[!is.na(res_df$OR)] > 0))
  expect_true(all(res_df$OR_lo[!is.na(res_df$OR_lo)] > 0))
  expect_true(all(res_df$OR_hi[!is.na(res_df$OR_hi)] > 0))
  
  ## each row corresponds to a valid variant + outcome label
  for (i in seq_len(nrow(res_df))) {
    expect_true(res_df$table_id[i] %in% names(df_list_small))
    expect_equal(res_df$outcome[i], "phenotype")
  }
})

test_that("build_model (linear) fits models and returns a well-formed summary", {
  vcf_path <- system.file(
    "extdata", "toy.vcf",
    package = "VcfAssociation",
    mustWork = TRUE
  )
  
  geno_data <- read_vcf(vcf_path)$genotypes
  pheno_df  <- generate_phenotype(
    vcf_path,
    chrom = "chr12",
    pos   = 11161,
    model = "additive"
  )
  
  tmpfile <- tempfile(fileext = ".csv")
  write.csv(pheno_df, tmpfile, row.names = FALSE)
  ph <- read_phenotypes(tmpfile, id_col = "sample", genotypes = geno_data)
  
  variants_df   <- unique(ph$merged[, c("CHROM", "POS")])[1:10, ]
  df_list_small <- prepare_list(
    ph,
    variants   = variants_df,
    phenotypes = "phenotype"
  )
  
  res_df <- build_model(
    df_list_small,
    outcomes = "phenotype",
    covars   = character(),
    model    = "linear"
  )
  
  ## basic structure
  expect_true(is.data.frame(res_df))
  expect_equal(nrow(res_df), length(df_list_small) * 1L)
  expect_setequal(res_df$table_id, names(df_list_small))
  
  ## expected columns
  expect_true(all(c(
    "table_id", "outcome", "model", "n",
    "term", "beta", "se", "p",
    "ci_lo", "ci_hi",
    "OR", "OR_lo", "OR_hi",
    "error"
  ) %in% names(res_df)))
  
  ## labels
  expect_equal(unique(res_df$outcome), "phenotype")
  expect_equal(unique(res_df$model), "linear")
  
  ## n is non-negative when not NA
  expect_true(all(is.na(res_df$n) | res_df$n >= 0))
  
  ## p-values within [0, 1] when non-NA
  expect_true(all(res_df$p[!is.na(res_df$p)] >= 0))
  expect_true(all(res_df$p[!is.na(res_df$p)] <= 1))
  
  ## OR columns should be NA for linear models
  expect_true(all(is.na(res_df$OR)))
  expect_true(all(is.na(res_df$OR_lo)))
  expect_true(all(is.na(res_df$OR_hi)))
})