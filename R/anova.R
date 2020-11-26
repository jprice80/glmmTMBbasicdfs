#' ANOVA Tables using the F distribution for glmmTMB models
#'
#' @description The \code{anova} function will perform an analysis of variance on glmmTMB
#'     objects using the F distribution. Denominator degrees of freedom calculations are
#'     performed using either the \code{\link{nlme}}, \code{inner-outer}, or \code{containment}
#'     methods.
#'
#' @param model a \code{\link{glmmTMB}} model object
#' @param type type of test, \code{II, III, 2, 3}. Roman numerals are equivalent
#'     to the corresponding Arabic numerals.
#' @param ddf the method for computing the denominator degrees of freedom and F-statistics.
#'     ddf = \code{\link{nlme}} corresponds to inner-outer denominator degrees of freedom
#'     rules used by \code{\link{nlme}}. ddf = \code{inner-outer} primarily corresponds
#'     to classic \code{\link{nlme}} denominator degrees of freedom calculations but
#'     makes additional corrections for random slope models. \code{containment} corresponds
#'     to the SAS default degrees of freedom calculations for generalized linear mixed models.
#'     Default is \code{inner-outer}.
#' @param test.statistic designates whether the \code{F} or \code{Chisq} distribution
#'     is used to compute the ANOVA tables.
#' @param contr_sum a \code{boolean} indicating whether to apply 'sum to zero'
#'     contrasts when performing the ANOVA calculation. Default is \code{TRUE}.
#'
#' @details Insert details here later
#' @note Insert notes here later
#' @references Insert references here later
#'
#' @return a \code{data.frame}
#' @export
#'
#' @examples Insert examples here later
anova <- function(model, type = 3, ddf = "inner-outer", test.statistic="F", contr_sum = TRUE){

  if (class(model) != "glmmTMB") {
    stop ("Only glmmTMB models are supported")
  }

  #TODO: check to see if this makes a difference
  if (contr_sum == TRUE){
    current_contrast_settings <- options("contrasts")
    options(contrasts = c("contr.sum", "contr.poly"))
    on.exit(options(current_contrast_settings))
  }

  if (type == "III" || type == 3) {
    type = 3
  } else if (type == "II" || type == 2) {
    type = 2
  } else {
    stop ("Specified type not supported at this time")
  }

  if(test.statistic == "F"){
    if(ddf == "nlme") {
      aov_out <- nlme_aov(model, type)
    } else if (ddf == "inner-outer") {
      aov_out <- inner_outer_aov(model, type)
    } else if (ddf == "containment") {
      aov_out <- containment_aov(model, type)
    } else {
      stop ("Only nlme, inner-outter, and containment methods are supported at this time")
    }
  } else if(test.statistic == "Chisq") {
    aov_out <- suppressForeignCheck(glmmTMB:::Anova.glmmTMB(model, type = type))
  } else {
    cat("Only F and Chisq test statistics are supported at this time")
  }

  return(aov_out)
}
