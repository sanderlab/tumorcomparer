#' Run Shiny App
#' 
#' @param port the port number to use (Default: 3838)
#' 
#' @examples
#' \dontrun{
#'     runShinyApp()
#' }
#' 
#' @concept tumorcomparer
#' @export
run_shiny_app <- function (port=3838) {
  shiny::runApp(system.file('shinyApp', package='tumorcomparer'), port=port)	
}
