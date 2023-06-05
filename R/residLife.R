
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

residLife <- function(data, life, p=.5, type = 'mean') {
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
fsr = flexsurvreg(Surv(recyrs, censrec)~ group, data = bc, dist = 'weibull')
residLife(fsr, 4)
fsr$dlist$name
