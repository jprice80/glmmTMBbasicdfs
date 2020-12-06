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
      cls <- class(model$frame[ ,trm])
      cls <- cls[length(class(model$frame[ ,trm]))]
      prop_out$varclass[i] <- cls
    } else {
      prop_out$varclass[i] <- as.character(NA)
    }
  }

  # Determine if we have a random int, slope, or int_and_slope
  prop_out$termtype <- as.character(NA)
  for(i in 1:nrow(prop_out)) {

    slope_term <- prop_out$slopeterm[i]
    covar <- prop_out$covar[i]

    # Identify the slope class
    if(slope_term != "1" && is.na(slope_term) == FALSE){

      containsneg1 <- grepl("- 1", slope_term, fixed = TRUE)
      containspls0 <- grepl("+ 0", slope_term, fixed = TRUE)

      eval_slope <- as.character(NA)
      if(containsneg1 == TRUE) {
        eval_slope <- trimws(strsplit(as.character(slope_term), "\\- ")[[1]][1])
      } else if(containspls0 == TRUE) {
        eval_slope <- trimws(strsplit(as.character(slope_term), "\\+ ")[[1]][1])
      } else {
        eval_slope <- as.character(slope_term)
        rslope_int <- TRUE
      }

      if(containsneg1 == TRUE || containspls0 == TRUE || rslope_int == TRUE) {
        fixed <- prop_out[prop_out$vartype=="fixed", ]
        slope_class <- fixed[which(fixed$terms == eval_slope), "varclass"]
      } else {
        slope_class <- "factor"
      }

      # Catch any slope effects that are not specified as fixed effects in prop_out
      if(identical(slope_class, character(0))){
        slope_class <- class(model$frame[ ,eval_slope])
      }
    }




    #Need to fix this for categorical random slope

    if(is.na(slope_term) == FALSE) {
      #if(slope_term == "1" || slope_class == "factor")
      if(slope_term == "1") {
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


#====================================== Formula Writer ===================================================

#' AOV slope and slope/intercept models equation developer
#'
#' @param model
#'
#' @return a \code{\link{formula}}
#' @keywords internal
aov_formula_writer <- function(model){
  form <- model$modelInfo$allForm$combForm
  fixed <- unique(attr(terms(model), "term.labels"))
  random <- unique(model$modelInfo$grpVar)
  y_name<- names(model$modelInfo$respCol)
  N <- model$modelInfo$nobs
  intercept<-attr(model$modelInfo$terms$cond$fixed, "intercept")
  num_random_eq <- length(model$modelInfo$reStruc$condReStruc)

  #TODO: check no intercept models to make sure this works
  if(intercept == 1){
    fixed_eq <- paste(fixed, collapse = "+")
  } else {
    fixed_eq <- paste("-1 + ", paste(fixed, collapse = "+"), sep = "+")
  }



  # Random Effects equation developer
  random_vec <- c()
  for(i in 1:num_random_eq) {

    eq <- unlist(names(model$modelInfo$reStruc$condReStruc)[i])

    # Replace -1 with +0
    containsneg1 <- grepl("- 1", eq, fixed = TRUE)
    if(containsneg1 == TRUE){
      eq <- gsub("- 1", "+ 0", eq)
    }

    terms <- trimws(strsplit(as.character(eq), "\\+")[[1]])
    # Identify the type of random effects we have
    if(trimws(strsplit(terms, "\\|")[[1]])[[1]] == "1") {
      model_type <- "int"
    } else if(any(grepl("0", terms) == TRUE) || any(grepl("-1", terms) == TRUE)) {
      model_type <- "slope"
    } else {
      model_type <- "int_slope"
    }


    if(model_type == "int"){
      rterm <- trimws(strsplit(terms, "\\|")[[1]])[[2]]
      random_vec[length(random_vec)+1] <- rterm
    } else if(model_type == "slope"){

      # Identify the right hand side of the bar
      # This is the appropriate intercept term
      for(j in 1:length(terms)){
        is_bar <- grep("\\|", terms[j])

        if(length(is_bar) > 0) {
          rhs <-trimws(strsplit(terms[j], "\\|")[[1]])[2]
        }
      }

      # Identify the left hand side of the bar and cross each term with the right hand side
      # This is the random slope side
      for(j in 1:length(terms)) {
        lhs <- trimws(strsplit(terms[j], "\\|")[[1]])[1]
        is_bar <- grep("\\|", terms[j])

        if(length(lhs) != 0 && lhs != "0") {
          rterm <- paste(lhs,rhs,sep=":")
          random_vec[length(random_vec)+1] <- rterm
        }
      }

    } else if (model_type == "int_slope") {

      # Identify the right hand side of the bar
      # This is the appropriate intercept term
      for(j in 1:length(terms)){
        is_bar <- grep("\\|", terms[j])

        if(length(is_bar) > 0) {
          rhs <-trimws(strsplit(terms[j], "\\|")[[1]])[2]

          # Place in the rhs for the int term
          random_vec[length(random_vec)+1] <- rhs
        }
      }

      # Identify the left hand side of the bar and cross each term with the right hand side
      # This is the random slope side
      for(j in 1:length(terms)) {
        lhs <- trimws(strsplit(terms[j], "\\|")[[1]])[1]
        is_bar <- grep("\\|", terms[j])

        if(length(lhs) != 0) {
          rterm <- paste(lhs,rhs,sep=":")
          random_vec[length(random_vec)+1] <- rterm
        }
      }
    }
  }

  random_eq <- paste(random_vec, collapse="+")
  form <- formula(paste(y_name, "~", paste(fixed_eq, random_eq, sep="+")))

  return(form)
}


#====================================== Exclusive Terms ===================================================


#' Identify Exclusive terms used to identify levels of random slope and intercept variables
#'
#' @param model
#' @param term
#'
#' @return a \code{\link{data.frame}}
#' @keywords internal
exclusive_terms<-function(model, term){
  y_name<- names(model$modelInfo$respCol)
  br <- best_rules(model)
  temp <- suppressMessages(br %>% filter(slopeterm ==  !!sym(term), rules == "slope")) %>% arrange(interceptterm)

  temp$length <- NA
  for(i in 1:nrow(temp)){
    temp$length[i] <- length(strsplit(temp$interceptterm[i], ":")[[1]])
  }

  temp <- temp %>% arrange(length)

  if(nrow(br) > 1){
    exclusive_terms <- c()
    for(j in 2:nrow(temp)){
      trm1 <- temp$interceptterm[j-1]
      trm2 <- temp$interceptterm[j]

      out <- strsplit(as.character(trm2), "\\:")[[1]] %in% strsplit(as.character(trm1), "\\:")[[1]]

      individual_terms <- strsplit(as.character(trm2), "\\:")[[1]]
      individual_terms <- individual_terms[out==FALSE]

      if(length(individual_terms) > 1){
        individual_terms <- paste0(individual_terms,":")
      }

      exclusive_terms[length(exclusive_terms)+1] <- individual_terms
    }
  }

  temp$exclusive_terms <- c(as.character(NA), exclusive_terms)

  return(temp)
}


