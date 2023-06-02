
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

residLife <- function(data, x, p=.5, type = 'mean') {
  if (data$dlist$name == 'gamma'){
    return(gamma.rl(data, x, p, type))
  }
  if (data$dlist$name == 'gompertz'){
    return(gompertz.rl(data, x, p, type))
  }
  if (data$dlist$name == 'llogis'){
    return(llogis.rl(data, x, p, type))
  }
  if (data$dlist$name == 'exp'){
    return(exp.rl(data, x, p, type))
  }
  if (data$dlist$name == 'lnorm'){
    return(lnorm.rl(data, x, p, type))
  }


  #add in other if statements
  else{
    return('Invalid Distribution!')
  }
}
fsr = flexsurvreg(Surv(recyrs, censrec)~ 1, data = bc, dist = 'exp')
residLife(fsr, 4)
