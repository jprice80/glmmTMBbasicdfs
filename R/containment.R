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

  #Next find the correct DF for each fixed effect term
  rules_ready <- data.frame(terms = row.names(TMBaov))
  rules2 <- rules %>% select(terms, rules) %>% distinct()
  rules2 <- left_join(rules_ready, rules2, by="terms")

  if(is.na(rules2$rules[1])){
    rules2$rules[1] <-"int"
  }

  for(i in 1:nrow(rules2)){
    rule <- rules2$rules[i]
    term <- rules2$terms[i]

    if(rule != "int"){
      random_vec <- c()
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
