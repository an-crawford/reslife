
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


#' Calculate Residual Life Values
#'
#' @param data Name of a 'flexsurvreg' object from which data is extracted
#' @param life Value of the 'given' used to calculate residual life
#' @param p percentile, default is .5
#' @param type can be 'mean', 'median', 'percentile', or 'all'. Default is 'mean'
#'
#' @return Residual life values
#' @export
#'
#' @examples
#' residLife(flexsurvreg, 6, .75, 'all')
#' residLife(flexsurvreg, 3, type = 'median')
residLife <- function(data, life, p=.5, type = 'mean', newdata = data.frame()) {
  #stopifnot(class(newdata)=='data.frame')
  if (data$dlist$name == 'gamma'){
    return(gamma.rl(data, life, p, type))
  }
  if (data$dlist$name == 'gompertz'){
    return(gompertz.rl(data, life, p, type))
  }
  if (data$dlist$name == 'llogis'){
    return(llogis.rl(data, life, p, type))
  }
  if (data$dlist$name == 'exp'){
    return(exp.rl(data, life, p, type))
  }
  if (data$dlist$name == 'lnorm'){
    return(lnorm.rl(data, life, p, type))
  }
  if (data$dlist$name == 'weibull.quiet'){
    return(weibull_mlr(data, life, p, type))
  }


  #add in other if statements
  else{
    return('Invalid Distribution!')
  }
}


library(flexsurv)
fsr = flexsurvreg(formula = Surv(recyrs, censrec) ~ group, data = bc,dist = "weibull")
newdat = t(c(censrec = 0, group = 'Poor'))
nd = as.data.frame(newdat)
#nd$V3 = as.factor(nd$V3)
predict(fsr, newdata = nd)


nd
newdat

residLife(fsr, 4, newdata= newdat)
class(nd)


n = c(10, 20, 30)


