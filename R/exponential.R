data("bc")
fsr = flexsurvreg(Surv(recyrs, censrec)~ 1, data = bc, dist = 'exp')
length(fsr$coefficients)
fsr$res[1][1]
fsr
