# tests/testthat/test-prepare_list.R
context("prepare_list")

test_that("prepare_list splits merged data into list of variant-level data frames", {
  vcf_path <- system.file("extdata", "toy.vcf", package = "VcfAssociation", mustWork = TRUE)
  # Prepare a merged phenotype-genotype table using the toy data
  geno_data <- read_vcf(vcf_path)$genotypes
  pheno_df <- generate_phenotype(vcf_path, chrom = "chr12", pos = 11161, model = "carrier")
  tmpfile <- tempfile(fileext = ".csv")
  write.csv(pheno_df, tmpfile, row.names = FALSE)
  ph <- read_phenotypes(tmpfile, id_col = "sample", genotypes = geno_data)
  
  # Case 1: No specific variants provided (should return one combined list element 'full')
  df_list_all <- prepare_list(ph, variants = NULL, phenotypes = "phenotype", covars = character())
  expect_true(is.list(df_list_all))
  expect_identical(names(df_list_all), "full")
  expect_true(is.data.frame(df_list_all$full))
  # The 'full' data frame should have all rows from ph$merged and include required columns
  expect_equal(nrow(df_list_all$full), nrow(ph$merged))
  expect_true(all(c("CHROM", "POS", "Dosage", "phenotype") %in% names(df_list_all$full)))
  
  # Case 2: Provide specific variants (should split data for each variant)
  # Choose two distinct variant positions from the merged data
  variants_df <- unique(ph$merged[, c("CHROM", "POS")])[1:2, ]
  df_list_sel <- prepare_list(ph, variants = variants_df, phenotypes = "phenotype")
  expect_equal(length(df_list_sel), 2)
  # Names of list elements should correspond to each variant (format "chr<CHROM>_pos<POS>")
  expect_setequal(names(df_list_sel),
                  paste0("chr", variants_df$CHROM, "_pos", variants_df$POS))
  # Each list element data frame should contain only one variant's data
  for (i in seq_len(nrow(variants_df))) {
    v_chr <- variants_df$CHROM[i]; v_pos <- variants_df$POS[i]
    df_i <- df_list_sel[[i]]
    expect_true(all(df_i$CHROM == v_chr) && all(df_i$POS == v_pos))
    expect_true("phenotype" %in% names(df_i) && "Dosage" %in% names(df_i))
    # Each variant subset should have one row per sample for that variant
    expect_equal(nrow(df_i), length(unique(df_i$sample)))
  }
})
