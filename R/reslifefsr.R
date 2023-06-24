
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


#' Calculate Residual Life Values Using 'flexsurvreg'
#'
#' @param data Name of a 'flexsurvreg' object from which data is extracted
#' @param life Value of the 'given' used to calculate residual life
#' @param p percentile, default is .5
#' @param type can be 'mean', 'median', 'percentile', or 'all'. Default is
#' 'mean'
#' @param newdata a data frame containing new data values to calculate residual
#' life for. Default is a blank data frame.
#'
#' @return Residual life values
#' @export
#'
#' @examples
#' reslifefsr(flexsurvreg, 6, .75, 'all')
#' reslifefsr(flexsurvreg, 3, type = 'median', newdata = df_new)
reslifefsr <- function(data, life, p=.5, type = 'mean', newdata = data.frame()) {
  stopifnot(class(newdata)=='data.frame')
  if (data$dlist$name == 'gamma'){
    return(gamma.rl(data, life, p, type, newdata))
  }
  if (data$dlist$name == 'gompertz'){
    return(gompertz.rl(data, life, p, type, newdata))
  }
  if (data$dlist$name == 'llogis'){
    return(llogis.rl(data, life, p, type, newdata))
  }
  if (data$dlist$name == 'exp'){
    return(exp.rl(data, life, p, type, newdata))
  }
  if (data$dlist$name == 'lnorm'){
    return(lnorm.rl(data, life, p, type, newdata))
  }
  if (data$dlist$name == 'weibull.quiet'){
    return(weibull_mlr(data, life, p, type, newdata))
  }
  if (data$dlist$name == 'gengamma.orig'){
    return(gengamma.orig(data, life, p, type, newdata))
  }
  if (data$dlist$name == 'gengamma'){
    return(gengamma(data, life, p, type, newdata))
  }


  #add in other if statements
  else{
    return('Invalid Distribution!')
  }
}



