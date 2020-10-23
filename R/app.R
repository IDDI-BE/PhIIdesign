
#' @title Launch Shiny app for Phase II clinical trials
#' @description Web application allowing to see designs.
#' Note that this requires the rmarkdown, flexdashboard, shiny, shinyWidgets, shinyBS, ggplot2 and shinipsum packages to be installed
#' @export
#' @examples
#' \dontrun{
#' library(rmarkdown)
#' library(flexdashboard)
#' library(shiny)
#' library(shinyWidgets)
#' library(shinyBS)
#' library(shinipsum)
#' library(ggplot2)
#' app_start()
#' }
app_start <- function(){
  path <- system.file(package = "PhIIdesign", "dashboard", "ph2design.Rmd")
  rmarkdown::run(path)
}
