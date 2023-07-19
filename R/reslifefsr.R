
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'
#
#
#
#Written by: Andrew Crawford and Zekai Wang


#' Calculate Residual Life Values Using a flexsurvreg() Object
#' @description Calculate residual life values using a flexsurvreg() object. Contains an option to supply new data
#' and returns the output as a vector.
#' @param obj Name of a flexsurvreg() object from which data is extracted.
#' @param life Value at which the user wants to calculate residual life. Given as a scalar.
#' @param p percentile, default is .5.
#' @param type can be 'mean', 'median', 'percentile', or 'all'. Default is
#' 'mean'.
#' @param newdata a data frame containing new data values to calculate residual
#' life for. Default is a blank data frame.
#'
#' @return A vector of residual life values
#' @export
#'
#' @examples
#' library(flexsurv)
#' fitg <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ 1,
#' data = ovarian, dist="gengamma")
#' reslifefsr(obj = fitg, life = 6, p= .75, type= 'all')
#'
#' fitg2 <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ age,
#' data = ovarian, dist="gengamma")
#' df_new = data.frame(age = 12)
#' reslifefsr(obj = fitg2, life = 3, type = 'median', newdata = df_new)
reslifefsr <- function(obj, life, p=.5, type = 'mean', newdata = data.frame()) {
  stopifnot(class(newdata)=='data.frame')
  if (obj$dlist$name == 'gamma'){
    return(gamma.rl(obj, life, p, type, newdata))
  }
  if (obj$dlist$name == 'gompertz'){
    return(gompertz.rl(obj, life, p, type, newdata))
  }
  if (obj$dlist$name == 'llogis'){
    return(llogis.rl(obj, life, p, type, newdata))
  }
  if (obj$dlist$name == 'exp'){
    return(exp.rl(obj, life, p, type, newdata))
  }
  if (obj$dlist$name == 'lnorm'){
    return(lnorm.rl(obj, life, p, type, newdata))
  }
  if (obj$dlist$name == 'weibull.quiet'){
    return(weibull_mlr(obj, life, p, type, newdata))
  }
  if (obj$dlist$name == 'gengamma.orig'){
    return(gengamma.orig(obj, life, p, type, newdata))
  }
  if (obj$dlist$name == 'gengamma'){
    return(gengamma(obj, life, p, type, newdata))
  }
  if (obj$dlist$name == 'genf'){
    return(genF(obj, life, p, type, newdata))
  }
  if (obj$dlist$name == 'genf.orig'){
    return(genF.orig(obj, life, p, type, newdata))
  }


  #add in other if statements
  else{
    return('Invalid Distribution!')
  }
}




