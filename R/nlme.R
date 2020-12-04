#' nlme ANOVA and Degrees of Freedom Calculations
#'
#' @param model a \code{\link{glmmTMB}} model object
#' @param type type of test, \code{II, III, 2, 3}. Roman numerals are equivalent
#'     to the corresponding Arabic numerals.
#'
#' @return a \code{data.frame}
#' @keywords internal
nlme_aov <- function(model = model, type = type){

  dc <- dataClasses(model)

  TMBaov <- suppressPackageStartupMessages(car::Anova(model, type=type))

  # Pull the DFs associated with each term
  basic_aov_dfs <- base_aov_dfs(model)

  # Correct the terms if no intercept is estimated
  if(row.names(TMBaov)[1] == "(Intercept)") {
    fixed <- basic_aov_dfs[basic_aov_dfs$vartype=="fixed", ]
  } else {
    fixed <- basic_aov_dfs[basic_aov_dfs$vartype=="fixed", ][-1, ]
  }

  random <- basic_aov_dfs[basic_aov_dfs$vartype=="random", ]

  #sort the data by size to make sure the right df is chosen
  random<-random[order(random$df),]

  # Apply the nlme DFs
  fixed$vartype <- NULL
  fixed$denDf <- NA
  for(i in 1:nrow(fixed)){
    ftrm <- strsplit(fixed$terms[i], ":")[[1]]
    datacls <- dc[match(ftrm, dc$terms), "class"]

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
  }


  #Complete output
  chisq <- as.vector(TMBaov$Chisq)
  nDF <- as.vector(TMBaov$Df)
  Fval <- chisq/nDF
  dDF <- fixed$denDf
  Pval <- pf(Fval, nDF, dDF, lower.tail = FALSE)

  aod <- data.frame(numDF = nDF, denDF = dDF, Fvalue = round(Fval, 2), pvalue = round(Pval, 4))
  row.names(aod) <- fixed$terms
  class(aod) <- c("bdf", "nlme", "data.frame")
  y_name<- names(model$modelInfo$respCol)

  if (type == 3) {
    attr(aod, "heading") <-  paste("Analysis of Deviance Table (Type III F-tests)", "\n\nResponse: ", y_name)
  } else if (type == 2){
    attr(aod, "heading") <-  paste("Analysis of Deviance Table (Type II F-tests)", "\n\nResponse: ", y_name)
  }

  return(aod)
}


#====================================== Summary ===================================================

#' nlme Parameter Estimates Table Calculations
#'
#' @param model a \code{\link{glmmTMB}} model object
#'
#' @return a \code{data.frame}
nlme_summary <- function(model = model) {

  # Set up basic summary table information
  basic_aov_dfs <- nlme_aov(model)
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
