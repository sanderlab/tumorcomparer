#' Run Shiny App
#' 
#' @param port the port number to use (Default: 3838)
#' 
#' @examples
#' \dontrun{
#'     runShinyApp()
#' }
#' 
#' @concept PanCanMet
#' @export
runShinyApp <- function (port=3838) {
    shiny::runApp(system.file('shinyApp', package='PanCanMet'), port=port)	
}
