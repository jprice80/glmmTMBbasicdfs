#' Performs basic ANOVA degrees of freedom calculations for input into downstream
#' functions.
#'
#' @param model a \code{\link{glmmTMB}} model object
#'
#' @return a \code{\link{data.frame}}
#' @keywords internal
base_aov_dfs<-function(model){

  # Model Matrix
  form <- model$modelInfo$allForm$combForm
  fixed <- unique(attr(terms(model), "term.labels"))
  random <- unique(model$modelInfo$grpVar)
  y <- length(model$frame[,model$modelInfo$respCol])
  y_name<- names(model$modelInfo$respCol)
  N <- model$modelInfo$nobs
  intercept<-attr(model$modelInfo$terms$cond$fixed, "intercept")

  # Define term numbers
  #vecterms <- c(fixed,random)
  #newform <- as.formula(paste(y_name,"~ ",paste(vecterms, collapse="+")))
  fullterms <- unique(attr(terms(form), "term.labels"))


  # Define model matrix
  mf1 <- model.frame(formula = form, data=model$frame)
  mm1 <- model.matrix(form,mf1)


  # Establish output dataframe and count the number of terms
  if(intercept==1){
    basic_dfs_out<-data.frame(terms=c("(Intercept)", fullterms))
  } else {
    basic_dfs_out<-data.frame(terms=fullterms)
  }

  basic_dfs_out$basic_dfs_effectnum<-seq(from = 0, to = (nrow(basic_dfs_out)-1))

  # QR decomposition information
  my_qr <- qr(mm1)
  pivot <- my_qr$pivot
  rank <- my_qr$rank

  # Identify the appropriate column indices and number of terms
  full_matrix_columns_index <- attr(mm1, "assign")
  full_matrix_names <- dimnames(mm1)[[2]]

  # Extract column terms corresponding to appropriate pivots to generate dfs
  valid_matrix_columns_index <- full_matrix_columns_index[pivot[1L:rank]]
  unique_index <- unique(valid_matrix_columns_index)
  nterms <- length(unique_index)

  if(nrow(basic_dfs_out) != nterms){
    warning("Model is overparameterized resulting. Computation continues in spite of this warning.
            Please consider respecifying the model. Interpret this model with extreme caution.")
  }

  # Drop any terms from ANOVA table not existing in pivots
  basic_dfs_out <- basic_dfs_out[match(unique_index, basic_dfs_out$basic_dfs_effectnum),]
  basic_dfs_out$basic_dfs_effectnum <- NULL


  # effects calc
  # Q <- qr.Q(my_qr)
  # effects <- t(Q) %*% model$frame[,1]
  # valid_effect_names <- full_matrix_names[pivot[1L:rank]]
  #
  # if (!is.null(effects)){
  #   effects <- as.matrix(effects)[seq_along(asgn), , drop = FALSE]
  #   effects <- data.frame(valid_effect_names, groups=full_matrix_columns_index[pivot[1L:rank]], effects)
  # }

  # calculate the appropriate basic aov degrees of freedom
  df<-c()
  for (i in seq(nterms)) {

    #identify the parameter estimates (effects) that this iteration corresponds to
    ai <- (valid_matrix_columns_index == unique_index[i])

    # Sum the total number of effects for each term
    df <- c(df, sum(ai))
  }

  basic_dfs_out<-cbind(basic_dfs_out, df)

  # Calculate the residual dfs
  if(intercept==1){
    resids <- N - sum(df, na.rm = TRUE)
  } else {
    resids <- N - sum(df[-1], na.rm = TRUE)
  }

  resids<-data.frame(terms="Residuals", df=resids)
  basic_dfs_out<-rbind(basic_dfs_out, resids)

  # Define fixed and random terms
  basic_dfs_out$vartype<-as.character(NA)

  for(i in 1:nrow(basic_dfs_out)){
    trm<-basic_dfs_out$terms[i]

    if(trm %in% fixed){
      basic_dfs_out$vartype[i]<-"fixed"
    } else if(trm =="(Intercept)") {
      basic_dfs_out$vartype[i]<-"fixed"
    } else {
      basic_dfs_out$vartype[i]<-"random"
    }
  }

  # Define whether there are associated random slope terms
  random_terms <- basic_dfs_out[which(basic_dfs_out$vartype == "random"),]
  random_terms <- random_terms[which(random_terms$terms != "Residuals"),]

  for(i in 1:nrow(random_terms)) {
    random_effect <- model$modelInfo$reStruc$condReStruc[i]
    covar <- names(model$modelInfo$reStruc$condReStruc[[i]]$blockCode)
    random_effect_name <- names(random_effect[1])
    lhs <- trimws(strsplit(as.character(random_effect_name), "\\| ")[[1]][1])
    random_terms$slopeterm[i] <- lhs
    random_terms$covar[i] <- covar
  }

  basic_dfs_out <- left_join(basic_dfs_out, random_terms, by=c("terms", "df", "vartype"))

  # Determine if we have a random int, slope, or int_and_slope
  basic_dfs_out$termtype <- as.character(NA)
  for(i in 1:nrow(basic_dfs_out)){

    slope_term <- basic_dfs_out$slopeterm[i]
    containsneg1 <- grepl("- 1", slope_term, fixed = TRUE)
    containspls0 <- grepl("+ 0", slope_term, fixed = TRUE)

    if(!is.na(slope_term) == TRUE){
      if(slope_term == "1"){
        basic_dfs_out$termtype[i] <- "int"
      } else if (containsneg1 == TRUE || containspls0 == TRUE){
        basic_dfs_out$termtype[i] <- "slope"
      } else {
        basic_dfs_out$termtype[i] <- "int_and_slope"
      }
    }
  }

  return(basic_dfs_out)
}


