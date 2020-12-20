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
  fullterms <- unique(attr(terms(form), "term.labels"))

  # Define model matrix
  mf1 <- model.frame(formula = form, data = model$frame)
  mm1 <- model.matrix(form, mf1)


  # Establish output dataframe and count the number of terms
  if(intercept==1){
    basic_dfs_out<-data.frame(terms=c("(Intercept)", fullterms))
  } else {
    basic_dfs_out<-data.frame(terms=fullterms)
  }

  basic_dfs_out$basic_dfs_effectnum<-seq(from = 0, to = (nrow(basic_dfs_out)-1))
  basic_dfs_out$df<-NA

  # QR decomposition information
  my_qr <- qr(mm1)
  pivot <- my_qr$pivot
  rank <- my_qr$rank

  # Identify the appropriate column indices and number of terms
  full_matrix_columns_index <- attr(mm1, "assign")
  full_matrix_names <- dimnames(mm1)[[2]]

  # Extract column terms corresponding to appropriate pivots to generate dfs
  valid_matrix_columns_index <- full_matrix_columns_index[pivot[1L:rank]]

  # calculate the appropriate basic aov degrees of freedom
  for (i in 1:nrow(basic_dfs_out)) {

    effnum <- basic_dfs_out$basic_dfs_effectnum[i]
    effname <- basic_dfs_out$terms[i]

    #identify the parameter estimates (effects) that this iteration corresponds to
    ai <- (valid_matrix_columns_index == effnum)

    # Sum the total number of effects for each term
    df <- sum(ai)

    if(!is.na(df) && df != 0){
      basic_dfs_out[which(basic_dfs_out$terms == effname),"df"] <- df
    } else {
      basic_dfs_out[which(basic_dfs_out$terms == effname),"df"] <- NA
    }
  }

  basic_dfs_out$basic_dfs_effectnum <- NULL

  #resid_val <- N - sum(basic_dfs_out$df, na.rm = TRUE)
  resid_val <- resid_calc(model)
  resids <- data.frame(terms="Residuals", df=resid_val)
  basic_dfs_out <- rbind(basic_dfs_out, resids)

  if(any(is.na(basic_dfs_out$df))){
    warning("Unable to determine fully determine degrees of freedom for this model using the specified method.
            Computation continues in spite of this warning. Interpret this model with caution.")
  }

  # Define fixed and random terms
  basic_dfs_out$vartype <- as.character(NA)
  for(i in 1:nrow(basic_dfs_out)) {
    trm<-basic_dfs_out$terms[i]

    if(trm %in% fixed) {
      basic_dfs_out$vartype[i]<-"fixed"
    } else if(trm =="(Intercept)") {
      basic_dfs_out$vartype[i]<-"fixed"
    } else {
      basic_dfs_out$vartype[i]<-"random"
    }
  }

  return(basic_dfs_out)
}


#====================================== Formula and Data basic dfs ===================================================


#' Performs basic ANOVA degrees of freedom calculations for input into downstream
#' functions with a formula and data supplied separately. This useful for random
#' slope / slope & int models
#'
#' @param formula \code{\link{formula}}
#' @param data a \code{\link{data.frame}}
#'
#' @return a \code{\link{data.frame}}
#' @keywords internal
data_aov_dfs <- function(formula, data){

  N <- nrow(data)

  # Model Matrix
  form <- formula

  # Define term numbers
  fullterms <- unique(attr(terms(form), "term.labels"))

  # Define model matrix
  mf1 <- model.frame(formula = form, data = data)
  mm1 <- model.matrix(form, mf1)

  intercept <- attr(formula, "Intercept")

  # Establish output dataframe and count the number of terms
  if(all(mm1[,1] == 1)){
    basic_dfs_out<-data.frame(terms=c("(Intercept)", fullterms))
  } else {
    basic_dfs_out<-data.frame(terms=fullterms)
  }

  basic_dfs_out$basic_dfs_effectnum<-seq(from = 0, to = (nrow(basic_dfs_out)-1))
  basic_dfs_out$df<-NA

  # QR decomposition information
  my_qr <- qr(mm1)
  pivot <- my_qr$pivot
  rank <- my_qr$rank

  # Identify the appropriate column indices and number of terms
  full_matrix_columns_index <- attr(mm1, "assign")
  full_matrix_names <- dimnames(mm1)[[2]]

  # Extract column terms corresponding to appropriate pivots to generate dfs
  valid_matrix_columns_index <- full_matrix_columns_index[pivot[1L:rank]]

  # calculate the appropriate basic aov degrees of freedom
  for (i in 1:nrow(basic_dfs_out)) {

    effnum <- basic_dfs_out$basic_dfs_effectnum[i]
    effname <- basic_dfs_out$terms[i]

    #identify the parameter estimates (effects) that this iteration corresponds to
    ai <- (valid_matrix_columns_index == effnum)

    # Sum the total number of effects for each term
    df <- sum(ai)

    if(!is.na(df) && df != 0){
      basic_dfs_out[which(basic_dfs_out$terms == effname),"df"] <- df
    } else {
      basic_dfs_out[which(basic_dfs_out$terms == effname),"df"] <- NA
    }
  }

  basic_dfs_out$basic_dfs_effectnum <- NULL

  resid_val <- N - sum(basic_dfs_out$df, na.rm = TRUE)
  #resid_val <- resid_calc(model)
  resids <- data.frame(terms="Residuals", df=resid_val)
  basic_dfs_out <- rbind(basic_dfs_out, resids)

  if(any(is.na(basic_dfs_out$df))){
    warning("Unable to determine fully determine degrees of freedom for this model using the specified method.")
  }

  return(basic_dfs_out)
}


#====================================== single random slope ===================================================

individual_rslope <- function(lhs, rhs, data) {

  form <- paste("~", paste(rhs, paste(lhs, rhs, sep = ":"), sep = "+"))
  form <- as.formula(form)

  # Define term numbers
  fullterms <- unique(attr(terms(form), "term.labels"))

  mf1 <- model.frame(formula = form, data = data)
  mm1 <- model.matrix(form, mf1)

  intercept <- attr(formula, "Intercept")

  # Establish output dataframe and count the number of terms
  if(all(mm1[,1] == 1)){
    basic_dfs_out<-data.frame(terms=c("(Intercept)", fullterms))
  } else {
    basic_dfs_out<-data.frame(terms=fullterms)
  }

  basic_dfs_out$basic_dfs_effectnum<-seq(from = 0, to = (nrow(basic_dfs_out)-1))
  basic_dfs_out$df<-NA

  # QR decomposition information
  my_qr <- qr(mm1)
  pivot <- my_qr$pivot
  rank <- my_qr$rank

  # Identify the appropriate column indices and number of terms
  full_matrix_columns_index <- attr(mm1, "assign")
  full_matrix_names <- dimnames(mm1)[[2]]

  # Extract column terms corresponding to appropriate pivots to generate dfs
  valid_matrix_columns_index <- full_matrix_columns_index[pivot[1L:rank]]

  # calculate the appropriate basic aov degrees of freedom
  for (i in 1:nrow(basic_dfs_out)) {

    effnum <- basic_dfs_out$basic_dfs_effectnum[i]
    effname <- basic_dfs_out$terms[i]

    #identify the parameter estimates (effects) that this iteration corresponds to
    ai <- (valid_matrix_columns_index == effnum)

    # Sum the total number of effects for each term
    df <- sum(ai)

    if(!is.na(df) && df != 0){
      basic_dfs_out[which(basic_dfs_out$terms == effname),"df"] <- df
    } else {
      basic_dfs_out[which(basic_dfs_out$terms == effname),"df"] <- 0
      warning("Unable to determine fully determine degrees of freedom for this model using the specified method.")
    }
  }

  basic_dfs_out$basic_dfs_effectnum <- NULL

  df <- basic_dfs_out[nrow(basic_dfs_out), "df"]

  return(df)
}


#====================================== single random intercept and slope ===================================================

individual_rint_rslope <- function(lhs, rhs, data) {

  form <- paste("~", paste(lhs, rhs, paste(lhs, rhs, sep = ":"), sep = "+"))
  form <- as.formula(form)

  # Define term numbers
  fullterms <- unique(attr(terms(form), "term.labels"))

  mf1 <- model.frame(formula = form, data = data)
  mm1 <- model.matrix(form, mf1)

  intercept <- attr(formula, "Intercept")

  # Establish output dataframe and count the number of terms
  if(all(mm1[,1] == 1)){
    basic_dfs_out<-data.frame(terms=c("(Intercept)", fullterms))
  } else {
    basic_dfs_out<-data.frame(terms=fullterms)
  }

  basic_dfs_out$basic_dfs_effectnum<-seq(from = 0, to = (nrow(basic_dfs_out)-1))
  basic_dfs_out$df<-NA

  # QR decomposition information
  my_qr <- qr(mm1)
  pivot <- my_qr$pivot
  rank <- my_qr$rank

  # Identify the appropriate column indices and number of terms
  full_matrix_columns_index <- attr(mm1, "assign")
  full_matrix_names <- dimnames(mm1)[[2]]

  # Extract column terms corresponding to appropriate pivots to generate dfs
  valid_matrix_columns_index <- full_matrix_columns_index[pivot[1L:rank]]

  # calculate the appropriate basic aov degrees of freedom
  for (i in 1:nrow(basic_dfs_out)) {

    effnum <- basic_dfs_out$basic_dfs_effectnum[i]
    effname <- basic_dfs_out$terms[i]

    #identify the parameter estimates (effects) that this iteration corresponds to
    ai <- (valid_matrix_columns_index == effnum)

    # Sum the total number of effects for each term
    df <- sum(ai)

    if(!is.na(df) && df != 0){
      basic_dfs_out[which(basic_dfs_out$terms == effname),"df"] <- df
    } else {
      basic_dfs_out[which(basic_dfs_out$terms == effname),"df"] <- 0
      warning("Unable to determine fully determine degrees of freedom for this model using the specified method.")
    }
  }

  basic_dfs_out$basic_dfs_effectnum <- NULL

  df <- basic_dfs_out[nrow(basic_dfs_out), "df"]

  return(df)
}

# residual calc
resid_calc <- function(model){

  resid_form <- aov_formula_writer(model)
  df_resid <- data_aov_dfs(resid_form, data = model$frame)
  resid <- df_resid$df[nrow(df)]
  return(resid)
}


