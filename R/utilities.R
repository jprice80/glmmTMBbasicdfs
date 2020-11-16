#====================================== Data Classes ===================================================


#' Extract variable classes
#'
#' @param model
#'
#' @return a \code{\link{data.frame}}
#' @keywords internal
dataClasses <- function(model){
  temp <- attr(model$modelInfo$terms$cond$fixed, "dataClasses")
  out <- data.frame(terms=names(temp), class=as.vector(temp))
  return(out)
}

#====================================== Model Properties ===================================================


#' Define Model Properties for Various Terms
#'
#' @param model
#'
#' @return a \code{\link{data.frame}}
#' @keywords internal
model_properties <- function(model){

  form <- model$modelInfo$allForm$combForm
  fixed <- unique(attr(terms(model), "term.labels"))
  n_random <- length(model$modelInfo$reStruc$condReStruc)
  prop_out <- data.frame(terms=fixed, slopeterm=NA, covar=NA)

  # Define random slope terms
  for(i in 1:n_random) {
    random_effect <- model$modelInfo$reStruc$condReStruc[i]
    covar <- names(model$modelInfo$reStruc$condReStruc[[i]]$blockCode)
    random_effect_name <- names(random_effect[1])
    lhs <- trimws(strsplit(as.character(random_effect_name), "\\| ")[[1]][1])
    rhs <- trimws(strsplit(as.character(random_effect_name), "\\| ")[[1]][2])

    lhs_terms <- strsplit(as.character(lhs), "\\+ ")[[1]]
    rslope <- (any(lhs_terms == 0) == TRUE)

    if(rslope == TRUE){
      for(j in 1:length(lhs_terms)){
        st <- trimws(lhs_terms[j])

        if(st != "0"){
          currow <- data.frame(terms = rhs, slopeterm = paste(st,"+ 0"), covar = covar)
          prop_out <- rbind(prop_out, currow)
        }
      }
    } else {
      for(j in 1:length(lhs_terms)){
        currow <- data.frame(terms = rhs, slopeterm =  trimws(lhs_terms[j]), covar = covar)
        prop_out <- rbind(prop_out, currow)
      }
    }
  }

  # Determine if we have a random int, slope, or int_and_slope
  prop_out$termtype <- as.character(NA)
  for(i in 1:nrow(prop_out)) {

    slope_term <- prop_out$slopeterm[i]
    containsneg1 <- grepl("- 1", slope_term, fixed = TRUE)
    containspls0 <- grepl("+ 0", slope_term, fixed = TRUE)

    if(!is.na(slope_term) == TRUE) {
      if(slope_term == "1") {
        prop_out$termtype[i] <- "int"
      } else if (containsneg1 == TRUE || containspls0 == TRUE) {
        prop_out$termtype[i] <- "slope"
      } else {
        prop_out$termtype[i] <- "int_and_slope"
      }
    }
  }

  # Define fixed and random terms
  prop_out$vartype <- as.character(NA)
  for(i in 1:nrow(prop_out)) {
    trm<-prop_out$terms[i]

    if(trm %in% fixed) {
      prop_out$vartype[i]<-"fixed"
    } else if(trm =="(Intercept)") {
      prop_out$vartype[i]<-"fixed"
    } else {
      prop_out$vartype[i]<-"random"
    }
  }

  # Identify the variable classes in the data frame
  prop_out$varclass <- as.character(NA)
  for(i in 1:nrow(prop_out)) {
    trm <- prop_out$terms[i]

    #exclude interactions since we just want the main data frame variables
    len <- length(strsplit(trm, ":")[[1]])

    if(trm != "(Intercept)" && trm != "Residuals" && len == 1){
      prop_out$varclass[i] <- class(model$frame[ ,trm])
    }
  }

  return(prop_out)
}
