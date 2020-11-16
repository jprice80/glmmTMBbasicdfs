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

  # Determine if we have a random int, slope, or int_and_slope
  prop_out$termtype <- as.character(NA)
  for(i in 1:nrow(prop_out)) {

    slope_term <- prop_out$slopeterm[i]

    # Identify the slope class
    if(slope_term != "1" && is.na(slope_term) == FALSE){

      containsneg1 <- grepl("- 1", slope_term, fixed = TRUE)
      containspls0 <- grepl("+ 0", slope_term, fixed = TRUE)

      if(containsneg1 == TRUE) {
        eval_slope <- trimws(strsplit(as.character(slope_term), "\\- ")[[1]][1])
      } else if(containspls0 == TRUE) {
        eval_slope <- trimws(strsplit(as.character(slope_term), "\\+ ")[[1]][1])
      } else {
        rslope_int <- TRUE
      }

      if(containsneg1 == TRUE || containspls0 == TRUE || rslope_int == TRUE) {
        fixed <- prop_out[prop_out$vartype=="fixed", ]
        slope_class <- fixed[which(fixed$terms == eval_slope), "varclass"]
      } else {
        slope_class <- "factor"
      }
    }


    if(!is.na(slope_term) == TRUE) {
      if(slope_term == "1" || slope_class == "factor") {
        prop_out$termtype[i] <- "int"
      } else if (containsneg1 == TRUE || containspls0 == TRUE) {
        prop_out$termtype[i] <- "slope"
      } else {
        prop_out$termtype[i] <- "int_and_slope"
      }
    }
  }

  return(prop_out)
}


#====================================== Best Rules ===================================================


#' Identify the Most Appropriate Rules to Apply to each Fixed Effect Term
#'
#' @param model
#'
#' @return a \code{\link{data.frame}}
#' @keywords internal
best_rules <- function(model){
  mp<-model_properties(model)

  fixed <- mp[mp$vartype == "fixed",]
  random <- mp[mp$vartype == "random",]

  # replace p0 n1 terms
  for(i in 1:nrow(random)) {

    termtype <- random$termtype[i]

    if(termtype == "slope"){

      slope_term <- random$slopeterm[i]

      containsneg1 <- grepl("- 1", slope_term, fixed = TRUE)
      containspls0 <- grepl("+ 0", slope_term, fixed = TRUE)

      if(containsneg1 == TRUE) {
        random$slopeterm[i] <- trimws(strsplit(as.character(slope_term), "\\- ")[[1]][1])
      } else if(containspls0 == TRUE) {
        random$slopeterm[i] <- trimws(strsplit(as.character(slope_term), "\\+ ")[[1]][1])
      }
    }
  }

  # Determine the best rules to apply to each fixed effect
  slope_term_vec <- random$slopeterm
  rules_out <- data.frame()
  for(i in 1:nrow(fixed)) {

    term <- fixed$terms[i]
    pos <-which(slope_term_vec %in% term)

    if(identical(pos, integer(0))){
      rules <- "int"
      slope_term <- "1"
      intercept_term <- as.character(NA)
      covar <- as.character(NA)
    } else {
      rules <- random$termtype[pos]
      slope_term <- random$slopeterm[pos]
      intercept_term <- random$terms[pos]
      covar <- random$covar[pos]
    }

    currow <- data.frame(terms = term, slopeterm=slope_term, interceptterm=intercept_term, covar, rules)
    rules_out <- rbind(rules_out, currow)
  }

  return(rules_out)
}

