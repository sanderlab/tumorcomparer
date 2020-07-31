#' Starts a stopwatch timer to measure performance
#' 
#' @param gcFirst a boolean whether to run garbage collection before starting watch
#' @param type the time of time to keep track of 
#' 
#' @concept tumorcomparer
#' @export
tic <- function(gcFirst = TRUE, type=c("elapsed", "user.self", "sys.self")) {
  type <- match.arg(type)
  assign(".type", type, envir=baseenv())
  if(gcFirst) gc(FALSE)
  tic <- proc.time()[type]         
  assign(".tic", tic, envir=baseenv())
  invisible(tic)
}

#' Stops a stopwatch timer to measure performance 
#' 
#' @concept tumorcomparer
#' @export
toc <- function() {
  type <- get(".type", envir=baseenv())
  toc <- proc.time()[type]
  tic <- get(".tic", envir=baseenv())
  print(toc - tic)
  invisible(toc)
}