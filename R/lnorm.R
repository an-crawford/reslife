data("bc")
fsr = flexsurvreg(Surv(recyrs, censrec)~ 1, data = bc, dist = 'lnorm')
length(fsr$coefficients)
fsr$res[1][1]
fsr
