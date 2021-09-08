library(pwr)

qf(.5, 10, 20, lower = F)
pf(.5, 10, 20, lower = F)

pwr.t.test(100, .1, .05, type = 'two.sample')
pwr.t.test(100, .1, .05, type = 'one.sample')
pwr.chisq.test(.1, 10, 10, .05)
pwr.f2.test(10, 20, .1, .05)
pwr.anova.test(2, 3, .1, .05)

ES.h(.1, .2)
ES.w1(c(.5, .5), c(.6, .4))
P <- matrix(c(.4, .2, .3, .1), nrow = 2)
ES.w2(P)