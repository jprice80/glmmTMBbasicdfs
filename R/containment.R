#' Containment degrees of freedom calculation
#'
#' @param a \code{\link{glmmTMB}} model object
#' @param type type type of test, \code{II, III, 2, 3}. Roman numerals are equivalent
#'     to the corresponding Arabic numerals.
#'
#' @return a \code{data.frame}
#'
#' @keywords internal
containment_aov <- function(model = model, type = type){
  y_name<- names(model$modelInfo$respCol)
  dc <- dataClasses(model)
  TMBaov <- suppressPackageStartupMessages(car::Anova(model, type=type))

  # Pull the DFs associated with each term
  basic_aov_dfs <- base_aov_dfs(model)

  # Setup the final output
  if(row.names(TMBaov)[1] == "(Intercept)") {
    df_output <- basic_aov_dfs[basic_aov_dfs$vartype=="fixed", ]
  } else {
    df_output <- basic_aov_dfs[basic_aov_dfs$vartype=="fixed", ][-1, ]
  }

  df_output$vartype <- NULL
  df_output$denDf <- NA

  # Identify rules for each term
  rules <- best_rules(model)
  basic_aov_dfs <- left_join(basic_aov_dfs, rules, by="terms")

  # Correct the terms if no intercept is estimated
  if(row.names(TMBaov)[1] == "(Intercept)") {
    fixed <- basic_aov_dfs[basic_aov_dfs$vartype=="fixed", ]
    fixed$rules[1] <- "int"
    fixed$slopeterm[1] <- 1
  } else {
    fixed <- basic_aov_dfs[basic_aov_dfs$vartype=="fixed", ][-1, ]
  }

  # Identify the random terms
  random <- basic_aov_dfs[basic_aov_dfs$vartype=="random", ]

  #sort the data by size to make sure the right df is chosen
  random<-random[order(random$df),]

  # Apply the inner-outer DFs for random intercept terms
  fixed$denDf <- NA
  fixed$vartype <- NULL
  fixed$covar <- NULL
  for(i in 1:nrow(fixed)){
    ftrm <- strsplit(fixed$terms[i], ":")[[1]]
    datacls <- dc[match(ftrm, dc$terms), "class"]
    rule <- fixed$rules[i]

    #=================================================== Random Intercept Part ============================================================

    if(rule == "int"){
      if(!is.na(datacls[1])){
        for(j in 1:nrow(random)){
          rtrm<-strsplit(random$terms[j], ":")[[1]]

          if(rtrm[1] != "Residuals"){

            if(all(datacls %in% "factor")==TRUE){

              # categorical variable containment

              if(all(ftrm %in% rtrm)==TRUE){
                fixed$denDf[i] <- random$df[j]
                break;
              } else {
                fixed$denDf[i] <- random$df[nrow(random)]
              }

            } else {
              fixed$denDf[i] <- random$df[nrow(random)]
            }
          } else {
            fixed$denDf[i] <- random$df[nrow(random)]
          }
        }
      } else {
        fixed$denDf[i] <- random$df[nrow(random)]
      }
    } else {
      fixed$denDf[i] <- NA
    }
  }

  # Place DFs in the output data frame
  for(i in 1:nrow(df_output)){
    trm <- df_output[i,1]
    pos <- match(trm, fixed$terms)
    denDf <- fixed$denDf[pos]

    df_output$denDf[i] <- denDf
  }


  #=================================================== Random Slope/Int and Slope Part ============================================================


  if(("slope" %in% rules$rules) || ("int_and_slope" %in% rules$rules)){
    form <- aov_formula_writer(model)
    df_output1 <- data_aov_dfs(form, data = model$frame)

    fullterms <- c(unique(attr(terms(form), "term.labels")), "Residuals")
    fullterms <- data.frame(terms = fullterms)

    df_output2 <- left_join(fullterms, df_output1, by="terms")

    # What to do with NA returned by aov?
    # Manually compute dfs?
    if(any(is.na(df_output2$df)) == TRUE) {
      na_terms <- df_output2[is.na(df_output2$df), "terms"]

      rules$both <- paste(rules$slopeterm, rules$interceptterm, sep = ":")
      manual_terms <- rules %>% select(slopeterm, interceptterm, rules, both) %>% filter(both == na_terms)

      df_vec = c()
      for(i in 1:nrow(manual_terms)){
        lhs <- manual_terms$slopeterm[i]
        rhs <- manual_terms$interceptterm[i]
        rule <- manual_terms$rules[i]

        if(rule == "int_and_slope") {
          df_vec[length(df_vec) + 1] <- individual_rint_rslope(lhs, rhs, data = model$frame)
        } else if (rule == "slope") {
          df_vec[length(df_vec) + 1] <- individual_rslope(lhs, rhs, data = model$frame)
        }
      }

      # Insert the minimum df value into the appropriate space
      df_output2[which(df_output2$terms == na_terms), "df"] <- min(df_vec, na.rm = TRUE)

    }
  }

  # Next find the correct DF for each fixed effect term
  rules_ready <- data.frame(terms = row.names(TMBaov))
  rules2 <- rules %>% select(terms, rules) %>% distinct()
  rules2 <- left_join(rules_ready, rules2, by="terms")

  if(is.na(rules2$rules[1])){
    rules2$rules[1] <-"int"
  }

  for(i in 1:nrow(rules2)){
    rule <- rules2$rules[i]
    term <- rules2$terms[i]

    # Find the current df as calculated by the random intercept computation
    # Consider this a possible df value
    # May need to be changed later to not consider this
    curdf <- df_output$denDf[i]

    if(rule != "int"){

      #start the vector with the random intercept value
      random_vec <- c(curdf)
      for(j in 1:nrow(df_output2)){
        interact <- grepl("\\:", df_output2$terms[j])
        if(interact == TRUE){
          inter_vec <- strsplit(df_output2$terms[j], ":")[[1]]

          if(term %in% inter_vec){
            both <- rules %>% select(slopeterm, interceptterm) %>% filter(slopeterm == term)

            for(k in 1:nrow(both)){
              both_vec <- strsplit(both$interceptterm[k], ":")[[1]]
              if(all(both_vec %in% inter_vec)) {
                random_vec[length(random_vec) + 1] <- df_output2$df[j]
              }
            }
          }
        }
      }

      df <- min(random_vec, na.rm = TRUE)

      # If random slope & slope_int terms were considered pick the minimum df
      if(i > 1) {
        term2 <- rules2$terms[i-1]
        if(term2 == term){
          df2 <- df_output$denDf[i-1]
          if(df < df2){
            df_output$denDf[i] <- df
          }
        } else {
          df_output$denDf[i] <- df
        }
      } else {
        df_output$denDf[i] <- df
      }
    }
  }

  #Complete output
  chisq <- as.vector(TMBaov$Chisq)
  nDF <- as.vector(TMBaov$Df)
  Fval <- chisq/nDF
  dDF <- df_output$denDf
  Pval <- pf(Fval, nDF, dDF, lower.tail = FALSE)

  aod <- data.frame(numDF = nDF, denDF = dDF, Fvalue = round(Fval, 2), pvalue = round(Pval, 4))
  row.names(aod) <- df_output$terms
  class(aod) <- c("bdf", "containment", "data.frame")

  if (type == 3) {
    attr(aod, "heading") <-  paste("Analysis of Deviance Table (Type III F-tests)", "\n\nResponse: ", y_name)
  } else if (type == 2){
    attr(aod, "heading") <-  paste("Analysis of Deviance Table (Type II F-tests)", "\n\nResponse: ", y_name)
  }

  return(aod)
}


#====================================== Summary ===================================================


#' Containment Parameter Estimates Table Calculations
#'
#' @param model
#'
#' @param model a \code{\link{glmmTMB}} model object
#'
#' @return a \code{data.frame}
containment_summary <- function(model = model) {

  # Set up basic summary table information
  basic_aov_dfs <- containment_aov(model)
  param_est <- as.data.frame(summary(model)$coef$cond)
  param_names <- row.names(param_est)
  est <- param_est$Estimate
  se <- param_est$`Std. Error`
  tval <- param_est$`z value`

  newParmEst <- data.frame(Estimate=est, Std.Error=se, df=NA, t.value=tval, Probt=NA)
  row.names(newParmEst) <- param_names

  for(i in 1:nrow(basic_aov_dfs)) {
    #query the anova table term
    aovtrm <- strsplit(row.names(basic_aov_dfs)[i], ":")[[1]]
    denDF <- basic_aov_dfs$df[i]

    for(j in 1:nrow(newParmEst)) {

      #query the summary table terms to compare to the anova term
      param_names <- row.names(newParmEst)[j]
      matchvec <- c()
      param_names_vec <- strsplit(param_names, ":")[[1]]

      if(length(param_names_vec) == length(aovtrm)) {

        # split up interactions and get a vector
        for(k in 1:length(param_names_vec)) {

          # Check for the correctness of each term in the vector
          amatch<-agrep(aovtrm[k], param_names_vec[k])

          if(length(amatch) == 0) {
            amatch <- 0
          }

          matchvec[length(matchvec)+1] <- amatch
        }

        if(all(matchvec == 1)) {
          newParmEst$df[j] <- denDF
          newParmEst$Probt[j] <- 2*pt(-abs(newParmEst$t.value[j]), df = denDF)
        }
      }
    }
  }

  #Assign the Intercept to be tested against the residual dfs
  if(row.names(newParmEst)[1] == "(Intercept)"){
    denDF<-basic_aov_dfs[nrow(basic_aov_dfs), "df"]

    newParmEst$df[1] <- denDF
    newParmEst$Probt[1] <- 2*pt(-abs(newParmEst$t.value[1]), df = denDF)
  }

  class(newParmEst) <- c("bdf", "nlme", "data.frame")
  y_name<- names(model$modelInfo$respCol)

  if (type == 3) {
    attr(aod, "heading") <-  paste("Conditional model:", "\n\nResponse: ", y_name)
  } else if (type == 2){
    attr(aod, "heading") <-  paste("Conditional model:", "\n\nResponse: ", y_name)
  }

  return(newParmEst)
}

