#====================================== Summary ===================================================

#' Summary methods
#'
#' @param object a \code{\link{bdf}} object
#' @export
summary.bdf <- function (object){
  c1 <- class(object)[1]
  c2 <- class(object)[2]

  if(c1 != "bdf"){
    stop(paste(object, "not an bdf object"))
  }

  if(c2 == "inner_outer"){

    out <- inner_outer_summary(object)
    print(out)

  } else if(c2 == "nlme") {
    out <- nlme_summary(object)
    print(out)
  }
}


#====================================== Summary ===================================================

#' Print methods
#'
#' @param object
#'
#' @export
print.bdf <- function (object){
  c1 <- class(object)[1]
  c2 <- class(object)[2]

  if(c1 != "bdf"){
    stop(paste(object, "not an bdf object"))
  }

  if(c2 == "inner_outer" || c2 == "nlme" || c2 == "containment"){

    out <- data.frame(object)

    out$pvalue <- ifelse(out$pvalue == 0.0000, as.character("<0.0001"),
                         as.character(format(out$pvalue, scientific = FALSE)))

    colnames(out) <- c("numDF", "denDF", "F value", "Pr(>F)")

    cat(attr(object, "heading"), "\n")
    print(out)
    #cat("---")

  }
}

