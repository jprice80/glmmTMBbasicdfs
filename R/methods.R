#====================================== Summary ===================================================

#' Summary methods
#'
#' @param object a \code{\link{bdf}} object
#' @export
# summary.bdf<-function (object){
#   c1 <- class(object)[1]
#   c2 <- class(object)[2]
#
#   if(c1 != "bdr"){
#     stop(paste(object, "not an bdf object"))
#   }
#
#   if(c2 == inner_outer){
#
#    print.bdf(object)
#
#   }
# }


#====================================== Summary ===================================================

#' Print methods
#'
#' @param object
#'
#' @return
#' @export
# print.bdf<-function (object){
#   c1 <- class(object)[1]
#   c2 <- class(object)[2]
#
#   if(c1 != "bdr"){
#     stop(paste(object, "not an bdf object"))
#   }
#
#   if(c2 == inner_outer){
#
#     #row.names(aod) <- fixed$terms
#     #class(aod) <- c("bdf", "inner_outer", "data.frame")
#     #attr(out, "heading") <- c("Analysis of variance with inner-outer degrees of freedom \n",
#     #                         paste("Response:", names(model$modelInfo$respCol)))
#
#     print(attr(object, "heading"))
#     print(object)
#     cat("---")
#
#   }
# }

