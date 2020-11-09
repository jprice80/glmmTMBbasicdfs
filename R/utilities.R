#====================================== Data Classes ===================================================


#' Extract variable classes
#'
#' @param model
#'
#' @return a \code{\link{glmmTMB}} model object
#'
#' @return a \code{\link{data.frame}}
#' @keywords internal
dataClasses <- function(model){
  temp <- attr(model$modelInfo$terms$cond$fixed, "dataClasses")
  out <- data.frame(terms=names(temp), class=as.vector(temp))
  return(out)
}
