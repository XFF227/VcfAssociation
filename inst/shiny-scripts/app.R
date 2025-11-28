library(shiny)
library(DT)
library(VcfAssociation)

ui <- fluidPage(
  titlePanel("VcfAssociation Shiny App"),
  sidebarLayout(
    sidebarPanel(
      h4("1. Data source"),
      radioButtons(
        "dataSource",
        "Choose data source:",
        choices = c("Use example data" = "example",
                    "Upload my own files" = "upload"),
        selected = "example"
      ),
      conditionalPanel(
        condition = "input.dataSource == 'upload'",
        fileInput("vcfFile",  "Upload VCF file",       accept = c(".vcf", ".vcf.gz")),
        fileInput("phenoFile","Upload phenotype file", accept = c(".csv", ".tsv")),
        helpText(
          "Phenotype file must contain a sample ID column named 'sample' ",
          "and at least one trait column (e.g. 'phenotype')."
        )
      ),
      conditionalPanel(
        condition = "input.dataSource == 'example'",
        helpText(
          "Example data: built-in 'toy.vcf' and a generated phenotype ",
          "based on a single variant (e.g. chr12:11161)."
        )
      )
    ),
    mainPanel(
      h4("2. Load data and run GWAS"),
      helpText(
        "Specify the phenotype column and click the button to load data ",
        "and run single-variant association using gwas_single()."
      ),
      textInput("phenoColInput", "Phenotype column", value = "phenotype"),
      actionButton("loadBtn", "Load data and run GWAS"),
      br(), br(),
      tabsetPanel(
        id = "gwasTabs",
        tabPanel(
          "GWAS results table",
          DTOutput("gwasTable")
        ),
        tabPanel(
          "GWAS Manhattan plot",
          plotOutput("manhattanPlot"),
          tags$div(
            "Figure 1. Manhattan plot showing −log₁₀(P) values across chromosomes.The red line indicates the genome-wide significance threshold",
            style = "text-align:center; font-style:italic; color:#555; margin-top:20px;"
          ))
      ),
      br(), br(),
      h4("3. Joint model for selected variants"),
      helpText(
        "After GWAS finishes, select one or more variants (CHROM:POS) ",
        "and run a joint model (multivariable regression)."
      ),
      selectizeInput(
        "variantSelect",
        "Select variants of interest (CHROM:POS):",
        choices  = NULL,
        multiple = TRUE,
        options  = list(maxItems = 20, placeholder = "Run GWAS first.")
      ),
      actionButton("runJointBtn", "Run joint model"),
      br(), br(),
      tabsetPanel(
        id = "jointTabs",
        tabPanel(
          "Joint model results table",
          DTOutput("jointTable")
        ),
        tabPanel(
          "Joint model forest plot",
          plotOutput("forestPlot"),
          tags$div(
            "Figure 2. Forest plot showing estimated effects (β or log(OR)) and 95% confidence intervals for selected variants. Variants labelled “Firth failed” showed non-convergence or separation in Firth logistic regression, and were therefore excluded from the joint model",
            style = "text-align:center; font-style:italic; color:#555; margin-top:20px;"
          ))
      )
    )
  )
)

server <- function(input, output, session) {
  
  rv <- reactiveValues(
    ph_list = NULL,
    res_df  = NULL
  )
  
  get_pheno_col <- reactive({
    if (is.null(input$phenoColInput) || input$phenoColInput == "") {
      "phenotype"
    } else {
      input$phenoColInput
    }
  })
  
  observeEvent(input$loadBtn, {
    pheno_col <- get_pheno_col()
    
    if (identical(input$dataSource, "upload")) {
      vcf_path <- if (!is.null(input$vcfFile)) {
        input$vcfFile$datapath
      } else {
        system.file("extdata", "toy.vcf",
                    package = "VcfAssociation",
                    mustWork = TRUE)
      }
      
      pheno_path <- if (!is.null(input$phenoFile)) {
        input$phenoFile$datapath
      } else {
        pheno_df <- generate_phenotype(
          vcf_path,
          chrom = "chr12",
          pos   = 11161,
          model = "carrier"
        )
        tmp <- tempfile(fileext = ".csv")
        write.csv(pheno_df, tmp, row.names = FALSE)
        tmp
      }
      
    } else {
      vcf_path <- system.file("extdata", "toy.vcf",
                              package = "VcfAssociation",
                              mustWork = TRUE)
      pheno_df <- generate_phenotype(
        vcf_path,
        chrom = "chr12",
        pos   = 11161,
        model = "carrier"
      )
      pheno_path <- tempfile(fileext = ".csv")
      write.csv(pheno_df, pheno_path, row.names = FALSE)
    }
    
    vcf <- read_vcf(vcf_path)
    geno_data <- vcf$genotypes
    
    ph_list <- read_phenotypes(
      pheno_path,
      id_col    = "sample",
      genotypes = geno_data
    )
    
    res_df <- gwas_single(
      ph_list,
      pheno_col = pheno_col
    )
    
    rv$ph_list <- ph_list
    rv$res_df  <- res_df
    
    output$gwasTable <- DT::renderDT({
      if (is.null(res_df) || nrow(res_df) == 0L) return(NULL)
      res_df
    },
    options = list(
      pageLength = 10,
      lengthMenu = c(10, 20, 50, 100)
    ))
    
    output$manhattanPlot <- renderPlot({
      if (is.null(res_df) || nrow(res_df) == 0L) return(NULL)
      manhattan_plot(res_df)
    })
    
    if (!is.null(res_df) && nrow(res_df) > 0L) {
      ordered <- res_df[order(res_df$p), ]
      labels  <- paste0(ordered$CHROM, ":", ordered$POS)
      names(labels) <- labels
      
      updateSelectizeInput(
        session,
        "variantSelect",
        choices  = labels,
        selected = intersect(c("chr12:10537", "chr12:12372", "chr12:12180"), labels),
        server   = TRUE
      )
    } else {
      updateSelectizeInput(
        session,
        "variantSelect",
        choices  = NULL,
        selected = NULL
      )
    }
  })
  
  observeEvent(input$runJointBtn, {
    pheno_col <- get_pheno_col()
    ph_list   <- rv$ph_list
    res_df    <- rv$res_df
    selected  <- input$variantSelect
    
    if (is.null(ph_list) || is.null(res_df) || length(selected) == 0L) {
      output$jointTable  <- DT::renderDT(NULL)
      output$forestPlot  <- renderPlot(NULL)
      return(NULL)
    }
    
    split_vals <- strsplit(selected, ":", fixed = TRUE)
    chroms     <- vapply(split_vals, `[`, character(1L), 1L)
    poses_chr  <- vapply(split_vals, `[`, character(1L), 2L)
    poses      <- suppressWarnings(as.numeric(poses_chr))
    
    variant_df <- data.frame(
      CHROM = chroms,
      POS   = poses,
      stringsAsFactors = FALSE
    )
    
    df_list <- prepare_list(
      ph_list,
      variants   = variant_df,
      phenotypes = pheno_col
    )
    
    merged_df <- ph_list$merged
    y <- merged_df[[pheno_col]]
    uniq_y <- unique(stats::na.omit(y))
    model_type <- if (length(uniq_y) == 2L) "logistic" else "linear"
    
    variant_result <- build_model(
      df_list  = df_list,
      outcomes = pheno_col,
      dosage   = "Dosage",
      covars   = character(),
      model    = model_type
    )
    
    output$jointTable <- DT::renderDT({
      variant_result
    },
    options = list(
      pageLength = 10,
      lengthMenu = c(10, 20, 50, 100)
    ))
    
    output$forestPlot <- renderPlot({
      forest_plot(variant_result)
    })
  })
}

shinyApp(ui = ui, server = server)