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
print.bdf <- function (object){
  c1 <- class(object)[1]
  c2 <- class(object)[2]

  if(c1 != "bdf"){
    stop(paste(object, "not an bdf object"))
  }

  if(c2 == "inner_outer"){

    out <- data.frame(object)

    colnames(out) <- c("numDF", "denDF", "F value", "Pr(>F)")


    cat(attr(object, "heading"), "\n")
    print(out)
    #cat("---")

  }
}

