# tests/testthat/test-manhattan_plot.R
test_that("manhattan_plot generates Manhattan plot with highlights and annotations", {
  # Simulate GWAS results for two chromosomes
  set.seed(100)
  df <- data.frame(
    CHROM = rep(c("chr1", "chr2"), each = 25),
    POS   = rep(1:25, times = 2),
    p     = runif(50)
  )
  # Add an SNP identifier column and create some notable p-values
  df$SNP <- paste0("rs", 1:50)
  df$p[1] <- 1e-9   # very significant SNP
  df$p[2] <- 1e-6   # suggestive SNP
  
  # Draw Manhattan plot, highlighting two SNPs and annotating the top hit
  gp <- manhattan_plot(df, chr_col = "CHROM", pos_col = "POS", p_col = "p",
                       highlight_snps = df$SNP[1:2], snp_col = "SNP", annotate_top_n = 1)
  
  # The result should be a ggplot object
  expect_s3_class(gp, "ggplot")
  # Expect the plot to have layers: points, significance lines, highlight points, and annotation text
  expect_equal(length(gp$layers), 5)
})
