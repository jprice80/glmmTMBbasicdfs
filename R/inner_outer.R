#' Inner-Outer ANOVA and Degrees of Freedom Calculations
#'
#' @param model a \code{\link{glmmTMB}} model object
#' @param type type of test, \code{II, III, 2, 3}. Roman numerals are equivalent
#'     to the corresponding Arabic numerals.
#'
#' @return a \code{data.frame}
inner_outer_aov <- function(model = model, type = type){
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

            if(length(rtrm) == 1 && length(ftrm) == 1){
              temp<-suppressMessages(model$frame %>% select(rtrm, ftrm) %>% group_by(!!sym(rtrm)) %>% summarise(count = n_distinct(!!sym(ftrm), na.rm = TRUE)))
            } else if(length(ftrm) > 1) {
              temp<-suppressMessages(model$frame %>% select(rtrm, ftrm) %>% group_by_at(rtrm) %>% summarise(count = n_distinct(ftrm, na.rm = TRUE)))
            } else {
              temp<-suppressMessages(model$frame %>% select(rtrm, ftrm) %>% group_by_at(rtrm) %>% summarise(count = n_distinct(!!sym(ftrm), na.rm = TRUE)))
            }

            # Check to see if we have all 1s
            # If we do, we know we have found the lowest level

            if(all(temp$count == 1)){
              fixed$denDf[i] <- random$df[j]
              break;
            } else {
              fixed$denDf[i] <- random[nrow(random),2]
            }
          } else {
            fixed$denDf[i] <- random[nrow(random),2]
          }
        }
      } else {
        fixed$denDf[i] <- random[nrow(random),2]
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




  #=================================================== Random Slope Part ============================================================


  if(("slope" %in% rules$rules) == TRUE){
    slope_terms <- unique(fixed$slopeterm)
    slope_terms <- slope_terms[!slope_terms=="1"]
    num_unique_slope_terms <- length(slope_terms)

    allterms <- unique(fixed$terms)

    if("(Intercept)" %in% allterms){
      allterms <- allterms[2:length(allterms)]
      fixedform <- paste(y_name, "~", paste(allterms, collapse="+"))
    } else {
      allterms <- allterms[2:length(allterms)]
      fixedform <- paste(y_name, "~ -1 +", paste(allterms, collapse="+"))
    }

    dfvec <- c()
    for(i in 1:num_unique_slope_terms) {
      slope_trm <- slope_terms[i]
      trms <- fixed[which(fixed$slopeterm == slope_trm),]

      # Try each possible intercept term specified for each slope level
      for(j in 1:nrow(trms)) {

        # Check to see if individual random terms exist in the current random intercept interaction term
        intercept_terms <- strsplit(as.character(trms$interceptterm[j]), "\\:")[[1]]
        current_random_terms <- c()
        for(k in 1:length(random$terms)) {
          test_terms <- strsplit(as.character(random$terms[k]), "\\:")[[1]]
          in_it <- all(test_terms %in% intercept_terms)

          if(in_it == TRUE){
            current_random_terms[length(current_random_terms) + 1] <- random$terms[k]
          }
        }

        # Drop the associated random intercept term since we are only performing a random slope model
        # Remove this for random slope and intercept models
        temp_term <- strsplit(current_random_terms[length(current_random_terms)], "\\:")[[1]]
        if(all(temp_term %in% intercept_terms)){
          current_random_terms <- current_random_terms[-length(current_random_terms)]
        }

        randomform <- paste(current_random_terms, collapse = "+")
        slopeterm <- paste(slope_trm, trms$interceptterm[j], sep = ":")
        slopeform <- as.formula(paste(fixedform, randomform, slopeterm, sep="+"))
        slope_vec <- strsplit(slopeterm, "\\:")[[1]]

        slope_aov_dfs <- data_aov_dfs(formula = slopeform, data = model$frame)

        # Find the appropriate DFs
        for(k in 1:nrow(slope_aov_dfs)) {
          temp_term <- slope_aov_dfs$terms[k]
          temp_vec <- strsplit(temp_term, "\\:")[[1]]

          if(all(slope_vec %in% temp_vec)== TRUE) {
            dfs <- slope_aov_dfs$df[k]
            dfvec[length(dfvec) + 1] <- dfs
            break;
          }
        }
      }

      df <- min(dfvec, na.rm = TRUE)
      df_output[which(df_output$terms == slope_trm), "denDf"] <- df
    }
  }


  #=================================================== Random Slope and Intercept Part ============================================================


  if(("int_and_slope" %in% rules$rules) == TRUE){
    slope_terms <- unique(fixed$slopeterm)
    slope_terms <- slope_terms[!slope_terms=="1"]
    num_unique_slope_terms <- length(slope_terms)

    allterms <- unique(fixed$terms)

    if("(Intercept)" %in% allterms){
      allterms <- allterms[2:length(allterms)]
      fixedform <- paste(y_name, "~", paste(allterms, collapse="+"))
    } else {
      allterms <- allterms[2:length(allterms)]
      fixedform <- paste(y_name, "~ -1 +", paste(allterms, collapse="+"))
    }

    dfvec <- c()
    for(i in 1:num_unique_slope_terms) {
      slope_trm <- slope_terms[i]
      trms <- fixed[which(fixed$slopeterm == slope_trm),]

      # Try each possible intercept term specified for each slope level
      for(j in 1:nrow(trms)) {

        # Check to see if individual random terms exist in the current random intercept interaction term
        intercept_terms <- strsplit(as.character(trms$interceptterm[j]), "\\:")[[1]]
        current_random_terms <- c()
        for(k in 1:length(random$terms)) {
          test_terms <- strsplit(as.character(random$terms[k]), "\\:")[[1]]
          in_it <- all(test_terms %in% intercept_terms)

          if(in_it == TRUE){
            current_random_terms[length(current_random_terms) + 1] <- random$terms[k]
          }
        }

        randomform <- paste(current_random_terms, collapse = "+")
        slopeterm <- paste(slope_trm, trms$interceptterm[j], sep = ":")
        slopeform <- as.formula(paste(fixedform, randomform, slopeterm, sep="+"))
        slope_vec <- strsplit(slopeterm, "\\:")[[1]]

        slope_aov_dfs <- data_aov_dfs(formula = slopeform, data = model$frame)

        # Find the appropriate DFs
        for(k in 1:nrow(slope_aov_dfs)) {
          temp_term <- slope_aov_dfs$terms[k]
          temp_vec <- strsplit(temp_term, "\\:")[[1]]

          if(all(slope_vec %in% temp_vec)== TRUE) {
            dfs <- slope_aov_dfs$df[k]
            dfvec[length(dfvec) + 1] <- dfs
            break;
          }
        }
      }

      df <- min(dfvec, na.rm = TRUE)
      df_output[which(df_output$terms == slope_trm), "denDf"] <- df
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
  class(aod) <- c("bdf", "inner_outer", "data.frame")

  if (type == 3) {
    attr(aod, "heading") <-  paste("Analysis of Deviance Table (Type III F-tests)", "\n\nResponse: ", y_name)
  } else if (type == 2){
    attr(aod, "heading") <-  paste("Analysis of Deviance Table (Type II F-tests)", "\n\nResponse: ", y_name)
  }

  return(aod)
}


#====================================== Summary ===================================================

#' Inner-Outer Parameter Estimates Table Calculations
#'
#' @param model a \code{\link{glmmTMB}} model object
#'
#' @return a \code{data.frame}
inner_outer_summary <- function(model = model) {

  # Set up basic summary table information
  basic_aov_dfs <- inner_outer_aov(model)
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
    denDF <- basic_aov_dfs$denDF[i]

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

  return(newParmEst)
}
