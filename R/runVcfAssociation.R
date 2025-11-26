# R/runVcfAssociation.R

#' Launch the VcfAssociation Shiny application
#'
#' This function launches the Shiny application bundled with the
#' \code{VcfAssociation} package. The app provides an interactive interface for
#' exploring VCF association analyses. Users can upload their own VCF and
#' phenotype files or run the analysis using the built-in example data.
#'
#' @return No return value; called for side effects. A Shiny application is
#'   started in the current R session.
#' @importFrom shiny runApp
#' @examples
#' if (interactive()) {
#'   runVcfAssociation()
#' }
#' @export
runVcfAssociation <- function() {
  appDir <- system.file("shiny-scripts", package = "VcfAssociation")
  if (appDir == "") {
    stop(
      "Could not find the Shiny application directory 'shiny-scripts'. ",
      "Please check that the package is installed correctly.",
      call. = FALSE
    )
  }
  shiny::runApp(appDir, display.mode = "normal")
}