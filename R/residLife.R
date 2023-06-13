
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
#' @param type can be 'mean', 'median', 'percentile', or 'all'. Default is
#' 'mean'
#' @param newdata a data frame containing new data values to calculate residual
#' life for. Default is a blank data frame.
#'
#' @return Residual life values
#' @export
#'
#' @examples
#' residLife(flexsurvreg, 6, .75, 'all')
#' residLife(flexsurvreg, 3, type = 'median')
residLife <- function(data, life, p=.5, type = 'mean', newdata = data.frame()) {
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


  #add in other if statements
  else{
    return('Invalid Distribution!')
  }
}



devtools::load_all()



library(flexsurv)
group = c("Medium", 'Good', "Poor")
age = c(43, 35, 39)
newdata = data.frame(group, age)
newd = data.frame(group)
fitg <- flexsurvreg(formula = Surv(futime, fustat) ~ age, data = ovarian, dist = "weibull")
fsr = flexsurvreg(formula = Surv(recyrs, censrec) ~ group, data = bc,dist = "weibull")
fsr1 = flexsurvreg(formula = Surv(recyrs, censrec) ~ 1, data = bc,dist = "weibull")
newbc <- bc
newbc$age <- rnorm(dim(bc)[1], mean = 65-scale(newbc$recyrs, scale=FALSE),sd = 5)
fsr2 =  flexsurvreg(Surv(recyrs, censrec) ~ group+age, data=newbc, dist = 'weibull')

residLife(fsr, 4)
residLife(fsr1, 4)
residLife(fsr, 4,p = .6, type = 'all')
residLife(fsr1, 4,p = .6, type = 'all')
residLife(fsr, 4,p = .6, type = 'all', newdata= newdata) #errors because input data does not match flexsurv
residLife(fsr, 4,p = .6, type = 'all', newdata= newd)
residLife(fsr1, 4,p = .6, type = 'all', newdata= newd)
residLife(fsr2, 4,p = .6, type = 'all', newdata= newdata)

