#' Inner-Outer degrees of freedom
#'
#' @param model a \code{\link{glmmTMB}} model object
#' @param type type of test, \code{II, III, 2, 3}. Roman numerals are equivalent
#'     to the corresponding Arabic numerals.
#'
#' @return a \code{data.frame}
inner_outer<-function(model = model, type = type){

  dc <- dataClasses(model)

  TMBaov <- car::Anova(model, type=type)

  #Pull the DFs associated with each term
  basic_aov_dfs <- base_aov_dfs(model)

  fixed <- basic_aov_dfs[basic_aov_dfs$vartype=="fixed", ]
  random <- basic_aov_dfs[basic_aov_dfs$vartype=="random", ]

  #sort the data by size to make sure the right df is chosen
  random <- random[order(random$df), ]

  #Apply the inner-outer DFs
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
  class(aod) <- c("bdf", "inner_outer", "data.frame")
  attr(aod, "heading") <-  paste("Analysis of Variance Table")

  return(aod)
}


#====================================== Summary ===================================================


